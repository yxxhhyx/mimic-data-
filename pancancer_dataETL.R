library(DBI)
library(dplyr)
library(tidyverse)
library(dbplyr)
library(comorbidity)
library(icdcoder)
library(data.table)

con <- dbConnect(
  RPostgres::Postgres(),
  dbname = "mimiciv",
  host = "localhost",
  port = 5432,
  user = "postgres",
  password = "7929389"
)


cancer <- dbGetQuery(con, "
SELECT *
FROM mimiciv_derived.all_patient_features_v4
WHERE cancer_stage IN ('onset', 'post-cancer', 'pre-cancer');
")

dod <- dbGetQuery(con, "SELECT subject_id, dod FROM mimiciv_hosp.patients;")

# 2. 合并进 cancer 表
cancer_full <- cancer |>
  dplyr::left_join(dod, by = "subject_id") |>
  dplyr::mutate(
    death_date = as.POSIXct(dod),
    died_within_30d = ifelse(!is.na(death_date) & 
                               death_date <= dischtime + lubridate::days(30), 1, 0),
    died_within_90d = ifelse(!is.na(death_date) & 
                               death_date <= dischtime + lubridate::days(90), 1, 0),
    died_after_discharge = ifelse(hospital_expire_flag == 0 &
                                    !is.na(death_date) & death_date > dischtime, 1, 0)
  )


cancer_attr <- cancer_full %>%
  group_by(subject_id) %>%
  arrange(admittime, .by_group = TRUE) %>%
  mutate(
    # 唯一区间归属：dod 落在 [dischtime, next_admittime) 或院内死亡
    died_in_hosp_attr  = ifelse(!is.na(death_date) & hospital_expire_flag == 1 & death_date <= dischtime, 1, 0),
    died_postdisch_attr = ifelse(!is.na(death_date) & death_date > dischtime &
                                   (is.na(lead(admittime)) | death_date <= lead(admittime)),
                                 1, 0),
    died_attr_flag = ifelse(died_in_hosp_attr == 1 | died_postdisch_attr == 1, 1, 0),
    
    # 仅在归属行上计算 30/90 天死亡，其余一律 0
    died_within_30d_attr = ifelse(died_attr_flag == 1 & !is.na(death_date) &
                                    death_date <= dischtime + days(30), 1, 0),
    died_within_90d_attr = ifelse(died_attr_flag == 1 & !is.na(death_date) &
                                    death_date <= dischtime + days(90), 1, 0)
  ) %>%
  ungroup()


cancer_dig <- dbGetQuery(con, "
  SELECT 
      d.subject_id,
      d.hadm_id,
      d.seq_num,
      d.icd_code,
      d.icd_version,
      dd.long_title
  FROM mimiciv_hosp.diagnoses_icd AS d
  LEFT JOIN mimiciv_hosp.d_icd_diagnoses AS dd
    ON d.icd_code = dd.icd_code
   AND d.icd_version = dd.icd_version
  WHERE d.hadm_id IN (
      SELECT hadm_id
      FROM mimiciv_derived.all_patient_features_v4
      WHERE cancer_stage IN ('onset', 'post-cancer', 'pre-cancer')
  );
")

diagnoses_clean_cancer <- cancer_dig |>
  mutate(icd_code = str_trim(icd_code),  # 去除首尾空格
         icd_code = gsub("[^0-9A-Za-z.]", "", icd_code))  # 只保留字母数字和点号

ICD_map_cancer <- diagnoses_clean_cancer |>
  filter(icd_version == 9) |>
  pull(icd_code) |>
  convICD("icd9") |>
  as_tibble() |>
  distinct(icd9, .keep_all = TRUE) |>
  
  # === 已有的手动映射修正 ===
  mutate(icd10 = case_when(
    icd9 == "0414" ~ "B9620",   # Enterococcus infection
    icd9 == "1733" ~ "C4490",   # Skin, unspecified site (default trunk)
    icd9 == "1732" ~ "C4490",   # Skin, unspecified site
    icd9 == "7863" ~ "R042",    # Hemoptysis
    icd9 == "2740" ~ "M109",    # Gout, unspecified
    icd9 == "2750" ~ "E8319",   # Disorder of copper metabolism
    icd9 == "7876" ~ "R197",    # Abnormal bowel sounds
    icd9 == "4881" ~ "J09X1",   # Influenza A/H1N1
    icd9 == "4251" ~ "I422",    # Other cardiomyopathy
    icd9 == "5185" ~ "J951",    # Postprocedural respiratory failure
    
    # === 新增48项 + 修正 ===
    icd9 == "1736" ~ "C4460",   # Malignant neoplasm skin, lower limb
    icd9 == "5997" ~ "R31",     # Hematuria
    icd9 == "2874" ~ "D696",    # Secondary thrombocytopenia
    icd9 == "2849" ~ "D619",    # Pancytopenia, unspecified
    icd9 == "4539" ~ "I8290",   # Venous thrombosis, unspecified vein
    icd9 == "2765" ~ "E87",     # Fluid overload
    icd9 == "5698" ~ "K928",    # Digestive complications NEC
    icd9 == "V8554" ~ "Z6841",  # BMI 40+
    icd9 == "7845" ~ "R471",    # Speech disturbance NEC
    icd9 == "6119" ~ "N649",    # Breast disorders NEC
    icd9 == "9998" ~ "T8089XA", # Transfusion reaction
    icd9 == "7951" ~ "R76",     # TB skin test reaction
    icd9 == "4446" ~ "I740",    # Abdominal aortic embolism
    icd9 == "7995" ~ "R458",    # Nervousness / emotional disturbance
    icd9 == "7470" ~ "Q258",    # Pulmonary artery anomaly
    icd9 == "3489" ~ "G939",    # Brain disorders NEC
    icd9 == "V1389" ~ "Z87898", # Personal history of other diseases
    icd9 == "5162" ~ "J841",    # Idiopathic fibrosing alveolitis
    icd9 == "V1089" ~ "Z859",   # Personal history of malignancy, unspecified
    icd9 == "3128" ~ "F6989",   # Behavioral problems NEC
    icd9 == "2794" ~ "D89",     # Autoimmune disease NEC
    icd9 == "7888" ~ "R39",     # Urinary symptoms
    icd9 == "7299" ~ "M7989",   # Soft tissue disorders NEC
    icd9 == "9690" ~ "T4396XA", # Antidepressant poisoning
    icd9 == "2882" ~ "D729",    # Disease of white blood cells
    icd9 == "V2502" ~ "Z30430", # Insertion of IUD
    icd9 == "78060" ~ "R509",   # Fever, unspecified
    icd9 == "1735" ~ "C4460",   # Skin, upper limb
    icd9 == "9971" ~ "J9581",   # Surgical complication, respiratory
    icd9 == "V451" ~ "Z992",    # Dependence on renal dialysis
    icd9 == "1734" ~ "C4439",   # Skin, neck/scalp
    icd9 == "5968" ~ "N329",    # Bladder disorder NEC
    icd9 == "2866" ~ "D689",    # Hemorrhagic disorder
    icd9 == "5128" ~ "J9312",   # Other spontaneous pneumothorax
    icd9 == "9995" ~ "T8069XA", # Serum reaction NEC
    icd9 == "7938" ~ "R918",    # Abnormal finding on lung field imaging
    icd9 == "E9271" ~ "X5089XA",# Accident from overexertion
    icd9 == "1730" ~ "C000",    # Malignant neoplasm of lip
    icd9 == "1738" ~ "C449",    # Skin, other specified sites
    icd9 == "5118" ~ "J90",     # Pleural effusion
    icd9 == "9984" ~ "T8110XA", # Postoperative shock
    icd9 == "1731" ~ "C4459",   # Skin, trunk
    icd9 == "1739" ~ "C4499",   # Skin, unspecified
    icd9 == "3108" ~ "F079",    # Nonpsychotic mental disorder, organic
    icd9 == "9708" ~ "T4396XA", # CNS stimulant poisoning
    icd9 == "E9961" ~ "Y3621XA",# War injury
    icd9 == "2399" ~ "D489",    # Neoplasm NOS
    icd9 == "1737" ~ "C4469",   # Malignant neoplasm skin, lower limb/foot, unspecified
    icd9 == "2398" ~ "D489",    # Neoplasm of uncertain behavior, other specified sites
    icd9 == "2766" ~ "E871",    # Electrolyte imbalance, NEC
    icd9 == "2841" ~ "D6109",   # Pancytopenia due to antineoplastic chemotherapy / other cause
    icd9 == "2865" ~ "D6832",   # Hemorrhagic disorder due to circulating anticoagulants
    icd9 == "2880" ~ "D709",    # Neutropenia, unspecified
    icd9 == "3488" ~ "G939",    # Other brain conditions
    icd9 == "4440" ~ "I7401",   # Arterial embolism, upper extremity
    icd9 == "4538" ~ "I8290",   # Other venous embolism/thrombosis
    icd9 == "5163" ~ "J8410",   # Idiopathic pulmonary hemosiderosis / fibrosis
    icd9 == "6118" ~ "N639",    # Other specified disorders of breast
    icd9 == "7473" ~ "Q282",    # Other anomalies of great veins
    icd9 == "7806" ~ "R509",    # Fever
    icd9 == "7889" ~ "R399",    # Unspecified urinary symptom
    icd9 == "7931" ~ "R912",    # Nonspecific abnormal radiological finding, lung field
    icd9 == "7955" ~ "R748",    # Abnormal tuberculin reaction, without active TB
    icd9 == "7992" ~ "R5381",   # Physical debility / fatigue
    icd9 == "9973" ~ "T8102XA", # Postoperative hemorrhage
    icd9 == "9974" ~ "T8109XA", # Postoperative infection
    icd9 == "9980" ~ "T8189XA", # Other surgical complications NEC
    icd9 == "E927"  ~ "X509XXA",# Excessive exertion, unspecified
    icd9 == "E995"  ~ "Y360XXA",# War operations involving explosion of marine weapons
    icd9 == "V109"  ~ "Z859",   # Personal history of malignancy (site unspecified)
    icd9 == "V138"  ~ "Z8789",  # Personal history of other diseases
    icd9 == "V251"  ~ "Z30402", # Encounter for contraceptive management, oral
    icd9 == "V403"  ~ "Z8719",  # Personal history of other mental disorder
    icd9 == "V854"  ~ "Z6814",  # BMI 40–44.9
    TRUE ~ icd10
  )) 

dig_b <- diagnoses_clean_cancer |>
  filter(icd_version == 9) |>
  left_join(ICD_map_cancer, by = c("icd_code" = "icd9")) |>
  select(1:3, 7, 5, 6) |>
  rename(icd_code = icd10)

cancer_dig_norm10 <- diagnoses_clean_cancer |>
  filter(icd_version == 10) |>
  bind_rows(dig_b)

diag_groups <- cancer_dig_norm10 %>%
  inner_join(select(cancer_attr, subject_id, hadm_id, cancer_stage),
             by = c("subject_id", "hadm_id")) %>%
  filter(cancer_stage %in% c('pre-cancer','onset','post-cancer')) |>
  mutate(icd_code_norm = str_sub(icd_code, 1, 4))


normalize_icd <- function(code) {
  code <- toupper(gsub("[^A-Z0-9]", "", code))  # 去掉非字母数字
  ifelse(
    grepl("^C", code),  # 若以 C 开头（肿瘤）
    substr(code, 1, 4), # 保留前4位
    substr(code, 1, 3)  # 否则保留前3位
  )
}

diag_groups_clean <- diag_groups %>%
  mutate(icd_code_norm = normalize_icd(icd_code)) %>%
  # 同一次住院去重
  distinct(subject_id, hadm_id, icd_code_norm, cancer_stage, .keep_all = FALSE)

# 频次过滤（减少稀疏；阈值可调：20/50/100）
keep_codes <- diag_groups_clean %>%
  count(icd_code_norm) %>%
  #filter(n >= 50) %>%
  pull(icd_code_norm)

diag_groups_clean <- diag_groups_clean %>% filter(icd_code_norm %in% keep_codes)

message("剩余 ICD token 数：", n_distinct(diag_groups_clean$icd_code_norm))
message("样本量（patients）：", n_distinct(diag_groups_clean$subject_id))
message("住院量（admissions）：", n_distinct(diag_groups_clean$hadm_id))

