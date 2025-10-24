# ============================================================
#  Tang 风格癌症共病网络分析管线（tidy 封装版）
# ============================================================

library(tidyverse)
library(igraph)
library(tidygraph)
library(ggraph)
library(patchwork)
library(ggrepel)
library(RColorBrewer)

# ------------------------------------------------------------
# 0️⃣ 基础设置：颜色 + 系统分类
# ------------------------------------------------------------
icd_system_levels <- c(
  "心血管系统", "呼吸系统", "消化系统", "内分泌代谢",
  "泌尿生殖系统", "肌肉骨骼系统", "神经系统",
  "恶性肿瘤", "精神行为障碍", "其他系统"
)
icd_system_colors <- c(
  "心血管系统"     = "#E41A1C",
  "呼吸系统"       = "#377EB8",
  "消化系统"       = "#984EA3",
  "内分泌代谢"     = "#4DAF4A",
  "泌尿生殖系统"   = "#A65628",
  "肌肉骨骼系统"   = "#FF7F00",
  "神经系统"       = "#999999",
  "恶性肿瘤"       = "#E6AB02",
  "精神行为障碍"   = "#66C2A5",
  "其他系统"       = "#B3B3B3"
)

# ------------------------------------------------------------
# 1️⃣ 网络构建函数（已优化）
# ------------------------------------------------------------
build_stage_graph <- function(df_stage,
                              id_col = "hadm_id",
                              node_min_prev = 0.01,
                              node_min_n = 50,
                              edge_min_n11 = 20,
                              edge_min_jaccard = 0.05,
                              edge_fdr_alpha = 0.05) {
  id_col <- rlang::sym(id_col) 
  N_stage <- df_stage %>% distinct(!!id_col) %>% nrow()
  node_stats <- df_stage %>%
    distinct(!!id_col, icd_code_norm) %>%
    count(icd_code_norm, name = "n_pat") %>%
    mutate(prev = n_pat / N_stage)
  keep_nodes <- node_stats %>%
    filter(prev >= node_min_prev, n_pat >= node_min_n) %>%
    pull(icd_code_norm)
  df_stage_f <- df_stage %>% filter(icd_code_norm %in% keep_nodes)
  if (nrow(df_stage_f) == 0) return(make_empty_graph())
  
  pairs <- df_stage_f %>%
    distinct(!!id_col, icd_code_norm) %>%
    group_by(!!id_col) %>%
    reframe(icds = list(sort(unique(icd_code_norm)))) %>%
    mutate(pairs = map(icds, ~ {
      v <- .x
      if (length(v) < 2) return(NULL)
      as.data.frame(t(combn(v, 2)), stringsAsFactors = FALSE) %>%
        setNames(c("a","b"))
    })) %>%
    select(-icds) %>%
    unnest(pairs)
  if (is.null(pairs) || nrow(pairs) == 0) return(make_empty_graph())
  
  e_counts <- pairs %>% count(a, b, name = "n11")
  marg <- df_stage_f %>% distinct(!!id_col, icd_code_norm) %>%
    count(icd_code_norm, name = "n1dot")
  e_stats <- e_counts %>%
    left_join(marg %>% rename(a = icd_code_norm, n1dot_a = n1dot), by = "a") %>%
    left_join(marg %>% rename(b = icd_code_norm, n1dot_b = n1dot), by = "b") %>%
    mutate(
      n10 = n1dot_a - n11,
      n01 = n1dot_b - n11,
      n00 = N_stage - n11 - n10 - n01,
      jaccard = n11 / (n11 + n10 + n01),
      lift = (n11 / N_stage) / ((n1dot_a / N_stage) * (n1dot_b / N_stage)),
      phi = ((as.numeric(n11) * as.numeric(n00)) - (as.numeric(n10) * as.numeric(n01))) /
        sqrt((as.numeric(n11 + n10)) *
               (as.numeric(n11 + n01)) *
               (as.numeric(n00 + n10)) *
               (as.numeric(n00 + n01)))
    ) %>%
    mutate(p_fisher = pmap_dbl(list(n11, n10, n01, n00),
                               ~ fisher.test(matrix(c(..1, ..2, ..3, ..4), nrow = 2))$p.value),
           q_fdr = p.adjust(p_fisher, method = "BH")) %>%
    filter(n11 >= edge_min_n11,
           jaccard >= edge_min_jaccard,
           q_fdr < edge_fdr_alpha)
  if (nrow(e_stats) == 0) return(make_empty_graph())
  
  g <- graph_from_data_frame(
    d = e_stats %>%
      transmute(from = a, to = b, weight = n11, jaccard, lift, phi, q_fdr),
    vertices = node_stats %>%
      filter(icd_code_norm %in% unique(c(e_stats$a, e_stats$b))) %>%
      transmute(name = icd_code_norm, n_pat, prev),
    directed = FALSE
  )
  simplify(g)
}

# ------------------------------------------------------------
# 2️⃣ Tang 风格三阶段绘图封装
# ------------------------------------------------------------
plot_tang_networks <- function(g_pre, g_onset, g_post, title_suffix = "") {
  p1 <- make_core_plot_tang_final_v2(g_pre,   paste0("Pre-Cancer 共病网络", title_suffix))
  p2 <- make_core_plot_tang_final_v2(g_onset, paste0("Onset 共病网络", title_suffix))
  p3 <- make_core_plot_tang_final_v2(g_post,  paste0("Post-Cancer 共病网络", title_suffix))
  (p1 + p2 + p3) + plot_layout(guides = "collect") &
    theme(
      legend.position = "bottom",
      legend.title = element_text(family = "STKaiti", size = 11, face = "bold"),
      legend.text  = element_text(family = "STKaiti", size = 9)
    )
}

# ------------------------------------------------------------
# 3️⃣ 全局指标汇总与可视化
# ------------------------------------------------------------
summarize_network_metrics <- function(g_list, labels = c("Pre", "Onset", "Post")) {
  metrics <- map2_dfr(g_list, labels, function(g, stage) {
    tibble(
      stage = stage,
      n_nodes = vcount(g),
      n_edges = ecount(g),
      density = edge_density(g),
      avg_degree = mean(degree(g)),
      transitivity = transitivity(g, type = "globalundirected", isolates = "zero")
    )
  })
  
  p <- metrics %>%
    pivot_longer(cols = -stage, names_to = "metric", values_to = "value") %>%
    mutate(stage = factor(stage, levels = labels)) %>%
    ggplot(aes(stage, value, group = 1, color = metric)) +
    geom_line() + geom_point(size = 3) +
    facet_wrap(~ metric, scales = "free_y") +
    theme_minimal(base_family = "STKaiti") +
    labs(title = "三阶段网络全局拓扑指标对比", x = "阶段", y = "值")
  
  list(summary = metrics, plot = p)
}

# ------------------------------------------------------------
# 4️⃣ 系统构成比例分析
# ------------------------------------------------------------
analyze_system_composition <- function(g_list, labels = c("Pre", "Onset", "Post")) {
  get_sys <- function(g, stage) {
    g %>% as_tbl_graph() %>% activate(nodes) %>%
      mutate(system = classify_icd_system(name)) %>%
      as_tibble() %>% count(system) %>%
      mutate(stage = stage, prop = n / sum(n))
  }
  df_sys <- map2_dfr(g_list, labels, get_sys)
  
  p <- ggplot(df_sys, aes(x = reorder(system, prop), y = prop, fill = stage)) +
    geom_col(position = "dodge") + coord_flip() +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal(base_family = "STKaiti") +
    ggtitle("各阶段共病网络的系统节点构成")
  
  list(data = df_sys, plot = p)
}

# ------------------------------------------------------------
# 5️⃣ 节点度变化趋势（复用前面函数）
# ------------------------------------------------------------
analyze_degree_trend <- function(g_pre, g_onset, g_post) {
  res <- plot_degree_trend(g_pre, g_onset, g_post)
  list(data = res$change_table, plot = res$plot)
}

# ------------------------------------------------------------
# 6️⃣ 主函数：一键运行完整管线
# ------------------------------------------------------------
run_full_pipeline <- function(df, id_col = "hadm_id",
                              exclude_mgmt = TRUE,
                              node_min_prev = 0.01,
                              node_min_n = 50,
                              edge_min_n11 = 20,
                              edge_min_jaccard = 0.05,
                              edge_fdr_alpha = 0.05) {
  
  # 筛选管理性诊断
  if (exclude_mgmt) {
    df <- df %>%
      filter(!grepl("^(Z|V|W|X|Y|E|R|S|T)", toupper(icd_code_norm)))
  }
  
  # 构建三阶段网络
  g_pre   <- build_stage_graph(filter(df, cancer_stage=="pre-cancer"),   id_col)
  g_onset <- build_stage_graph(filter(df, cancer_stage=="onset"),        id_col)
  g_post  <- build_stage_graph(filter(df, cancer_stage=="post-cancer"),  id_col)
  
  # 核心输出
  networks_plot <- plot_tang_networks(g_pre, g_onset, g_post)
  metrics <- summarize_network_metrics(list(g_pre, g_onset, g_post))
  systems <- analyze_system_composition(list(g_pre, g_onset, g_post))
  trend <- analyze_degree_trend(g_pre, g_onset, g_post)
  
  list(
    graphs = list(pre = g_pre, onset = g_onset, post = g_post),
    plots = list(networks = networks_plot, metrics = metrics$plot,
                 systems = systems$plot, trend = trend$plot),
    summaries = list(metrics = metrics$summary,
                     system = systems$data,
                     degree_change = trend$data)
  )
}


# 主分析（排除管理诊断）
res_main <- run_full_pipeline(diag_groups_clean, id_col = "hadm_id", exclude_mgmt = TRUE)

# 包含管理/外因诊断的敏感性分析
res_mgmt <- run_full_pipeline(diag_groups_clean, id_col = "hadm_id", exclude_mgmt = FALSE)

# === 管理 / 外因诊断过滤函数 ===
is_management_or_external <- function(code){
  code <- toupper(code)
  grepl("^(Z|V|W|X|Y|E)", code)
}
is_symptom_or_injury <- function(code){
  code <- toupper(code)
  grepl("^(R|S|T)", code)
}

# === 提取管理外因诊断子集 ===
diag_groups_mgmt_only <- diag_groups_clean %>%
  filter(is_management_or_external(icd_code_norm) | is_symptom_or_injury(icd_code_norm))

res_mgmt_only <- run_full_pipeline(
  diag_groups_mgmt_only,
  id_col = "hadm_id",
  exclude_mgmt = FALSE  #  此处必须为 FALSE，否则会被再次过滤掉
)

# 输出主结果图
res_main$plots$networks      # Tang-style 三阶段网络
res_mgmt$plots$networks
res_mgmt_only$plots$networks
res_main$plots$metrics       # 全局指标趋势
res_main$plots$systems       # 系统节点构成柱状图
res_main$plots$trend         # 节点度变化趋势图

# 查看全局指标表
res_main$summaries$metrics

# 查看系统构成表
res_main$summaries$system

# 查看节点度变化表
head(res_main$summaries$degree_change, 10)


