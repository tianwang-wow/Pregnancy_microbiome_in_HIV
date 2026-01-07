library(phyloseq)
library(microbiome)
library(lme4)
library(lmerTest)
library(matrixStats)
library(ggplot2)
library(dplyr)
dec = function(x, k = 3) {return(trimws(format(round(x,k),nsmall=k)))}
pvalue_dec = function(pvalue, threshold = 0.001, k = 3) { if (is.na(pvalue) | is.nan(pvalue)) { out = ' ' } else if (pvalue >= threshold) { out = paste0('p=', dec(pvalue, k)) } else if (pvalue < threshold) { out = paste0('p<', threshold) }; return(out)  }
value_perc = function(x, k1 = 3, k2 = 3) { out = paste0(dec(x,k1),' (',dec(x/sum(x)*100,k2),'%)'); names(out) = names(x); return(out) }
UniFrac_BC = function(OTUtab = OTUtab, w = 2) { return(as.matrix(stats::dist(OTUtab, method = 'manhattan'))/w) }
Long_func = function(Sample_keys = Sample_keys,
                     Variable = 'sample_visit',
                     Values = unique(Sample_keys$sample_visit))
{ 
  out_ls = list() 
  out_id = 0 
  all_samples = unique(Sample_keys$first_half_id) 
  for (s_id in 1:length(all_samples)) { 
    current_subject = all_samples[s_id] 
    subdata = subset(Sample_keys, first_half_id == current_subject) 
    temp_ids = which(subdata[, Variable] %in% Values) 
    if (length(temp_ids) > 0) { 
      out_id = out_id + 1 
      subdata = subdata[temp_ids, ] 
      out_ls[[out_id]] = subdata 
      names(out_ls)[out_id] = current_subject 
    }
  }
  return(out_ls)
}
remove_col_0_counts = function(tab) { 
  col_ids_0_counts = which(colSums(tab) == 0) 
  if (length(col_ids_0_counts) > 0) { tab = tab[, - col_ids_0_counts] } 
  return(tab)
}

Meta_data = Meta_data_BACKUP = as.data.frame(readxl::read_xlsx("./Raw_data/Meta Data (Tian).xlsx"))
rownames(Meta_data) = Meta_data$first_half_id
Sample_keys = Sample_keys_BACKUP = as.data.frame(readxl::read_xlsx("./Raw_data/Sample Keys (Updated) 03202023.xlsx"))
rownames(Sample_keys) = Sample_keys$sample_id
Sample_keys$sample_visit_details = paste0(Sample_keys$sample_visit, '_', 
                                          Sample_keys$visit_type)

temp_ids = which(Sample_keys$visit_type == '6 PPt') 
Sample_keys[which(Sample_keys$first_half_id == Sample_keys$first_half_id[temp_ids]), ] 
Sample_keys$visit_type[temp_ids] = '24 PPt' 
Sample_keys$sample_visit_details[temp_ids] = 'post partum_24 PPt' 

Sample_keys = Sample_keys[order(Sample_keys$first_half_id, 
                                Sample_keys$specimen_GA), ]

raw_ASV_tab = data.frame(read.table(file = './Raw_data/2022_RUPAK_STOOL_SILVA_asvs_Jan_26_2023.txt', 
                                    header = T, sep = '\t')) 
raw_ASV_tab = raw_ASV_tab[which(!is.na(raw_ASV_tab$Global.Spec.ID)), ] 

raw_ASV_tab = raw_ASV_tab[which(raw_ASV_tab$Global.Spec.ID %in% Sample_keys$sample_id), ] 

dup_IDs = raw_ASV_tab$Global.Spec.ID[ duplicated(raw_ASV_tab$Global.Spec.ID) ] 
for (d_id in 1:length(dup_IDs)) { 
  current_rows = which(raw_ASV_tab$Global.Spec.ID == dup_IDs[d_id]) 
  current_OTU_total = rowSums(raw_ASV_tab[current_rows, 10:ncol(raw_ASV_tab)]) 
  raw_ASV_tab$Global.Spec.ID[ current_rows[ which.min(current_OTU_total) ] ] = NA 
}
raw_ASV_tab = raw_ASV_tab[which(!is.na(raw_ASV_tab$Global.Spec.ID)), ] 
rownames(raw_ASV_tab) = raw_ASV_tab$Global.Spec.ID 
raw_ASV_tab = raw_ASV_tab[, 10:ncol(raw_ASV_tab)] 

ASV_classification = data.frame(read.table(file = './Raw_data/SILVA132_raw_01232023.csv', 
                                           header = T, sep = ',')) 
rownames(ASV_classification) = ASV_classification$X 
ASV_classification = ASV_classification[, -1] 

ps = phyloseq(otu_table(raw_ASV_tab, taxa_are_rows = F), tax_table(as.matrix(ASV_classification))) 
Genus_ct_tab = t(otu_table(microbiome::aggregate_taxa(ps, level = 'Genus', verbose = F))) 

row2remove = which(rowSums(Genus_ct_tab) < 1000) 
Genus_ct_tab = Genus_ct_tab[- row2remove, ] 

col2remove = which(colSums(Genus_ct_tab) == 0) 
Genus_ct_tab = Genus_ct_tab[, - col2remove] 

Genus_abun_tab = Genus_ct_tab / rowSums(Genus_ct_tab) 

Sample_keys = Sample_keys[which(Sample_keys$sample_id %in% rownames(Genus_abun_tab)), ] 


alpha_div = microbiome::alpha(t(Genus_ct_tab),
                              index = c('shannon', 'dominance_simpson'))
alpha_div$dominance_simpson = 1 - alpha_div$dominance_simpson 
other_alpha_div = estimate_richness(phyloseq(otu_table(Genus_ct_tab, taxa_are_rows = F))) 
rownames(other_alpha_div) = rownames(Genus_ct_tab) 
new_alpha_div = other_alpha_div[rownames(alpha_div), c("Chao1", "Fisher")] 
alpha_div = cbind(alpha_div, new_alpha_div) 
colnames(alpha_div)[which(colnames(alpha_div) == 'diversity_shannon')] = "Shannon" 
colnames(alpha_div)[which(colnames(alpha_div) == 'dominance_simpson')] = "Simpson"


alpha_div = alpha_div[Sample_keys$sample_id, ] 
Sample_keys = cbind(Sample_keys, alpha_div) 



Meta_data$education = factor(Meta_data$education, 
                             levels = c('None to primary', 'Middle to high school', 'Post high school to post graduate'))
Meta_data$post_high_school = NA
Meta_data$post_high_school[ which(Meta_data$education == 'Post high school to post graduate')] = 1
Meta_data$post_high_school[ which(Meta_data$education %in% c('None to primary', 'Middle to high school'))] = 0


Meta_data$muac_6mo_pp = as.numeric(Meta_data$muac_6mo_pp)

Meta_data$muacttrim = as.numeric(Meta_data$muacttrim)
Meta_data$muacstrim = as.numeric(Meta_data$muacstrim) 
Meta_data$bmi_zscore_6mo_cat = NA
Meta_data$bmi_zscore_6mo_cat[ which(Meta_data$bmi_zscore_6mo >= -2)] = 0
Meta_data$bmi_zscore_6mo_cat[ which(Meta_data$bmi_zscore_6mo < -2)] = 1

sample_ls = list()
SubjectID = Meta_data$first_half_id 
case_ids = which(Meta_data$hivst == 'Positive') 
control_ids = which(Meta_data$hivst == 'Negative') 
preg_ids = which(Meta_data$preg == 'yes') 
nonpreg_ids = which(Meta_data$preg == 'no') 

sample_ls$"Study 2" = SubjectID[which(Meta_data$study == 'sub-study')]
sample_ls$"Study 3" = SubjectID[which(Meta_data$study == 'stand-alone')]
sample_ls$"Study 2 & 3" = unique(union(sample_ls$"Study 2", sample_ls$"Study 3"))

sample_ls$"HIV+" = SubjectID[case_ids] 
sample_ls$"HIV-" = SubjectID[control_ids] 

sample_ls$"Pregnant" = SubjectID[preg_ids] 
sample_ls$"Nonpregnant" = SubjectID[nonpreg_ids] 

sample_ls$"HIV+ pregnant" = intersect(sample_ls$`HIV+`, sample_ls$Pregnant)
sample_ls$"HIV+ nonpregnant" = intersect(sample_ls$`HIV+`, sample_ls$Nonpregnant)
sample_ls$"HIV- pregnant" = intersect(sample_ls$`HIV-`, sample_ls$Pregnant)
sample_ls$"HIV- nonpregnant" = intersect(sample_ls$`HIV-`, sample_ls$Nonpregnant)

sample_ls$"Intense HIV+ pregnant" = intersect(sample_ls$"HIV+ pregnant", sample_ls$`Study 2 & 3`)
sample_ls$"Intense HIV- pregnant" = intersect(sample_ls$"HIV- pregnant", sample_ls$`Study 2 & 3`)

full_long_ls = Long_func(Sample_keys = Sample_keys, 
                         Variable = 'sample_visit')
for (s_id in 1:length(full_long_ls)) { 
  full_long_ls[[s_id]]$timepoint = NA 
  full_long_ls[[s_id]]$timepoint[ full_long_ls[[s_id]]$sample_visit %in% 'second trimester enrollment'] = 1
  full_long_ls[[s_id]]$timepoint[ full_long_ls[[s_id]]$sample_visit %in% c('third trimester enrollment', 'third trimester visit')] = 2
  full_long_ls[[s_id]]$timepoint[ full_long_ls[[s_id]]$sample_visit %in% 'post partum'] = 3
  full_long_ls[[s_id]]$timepoint[ full_long_ls[[s_id]]$sample_visit %in% 'non pregnant'] = 5
}
sample_ls$"2nd Trimester vs 3rd Trimester" = names(which(sapply(full_long_ls, function(x) { sum(c(1, 2) %in% x[, 'timepoint']) == 2 }))) 
sample_ls$"2nd Trimester vs 6 Month Ppt" = names(which(sapply(full_long_ls, function(x) { sum(c(1, 3) %in% x[, 'timepoint']) == 2 }))) 
sample_ls$"3rd Trimester vs 6 Month Ppt" = names(which(sapply(full_long_ls, function(x) { sum(c(2, 3) %in% x[, 'timepoint']) == 2 }))) 
sample_ls$"2nd Trimester vs 3rd Trimester common" = 
  sample_ls$"2nd Trimester vs 6 Month Ppt common" = 
  sample_ls$"3rd Trimester vs 6 Month Ppt common" = 
  intersect(intersect(sample_ls$"2nd Trimester vs 3rd Trimester", sample_ls$"2nd Trimester vs 6 Month Ppt"), sample_ls$"3rd Trimester vs 6 Month Ppt")

sample_ls$"2nd Trimester" = names(which(sapply(full_long_ls, function(x) { sum(c(1) %in% x[, 'timepoint']) == 1 })))
sample_ls$"3rd Trimester" = names(which(sapply(full_long_ls, function(x) { sum(c(2) %in% x[, 'timepoint']) == 1 })))
sample_ls$"6 Month Ppt" = names(which(sapply(full_long_ls, function(x) { sum(c(3) %in% x[, 'timepoint']) == 1 })))
sample_ls$"Nonpregnant" = names(which(sapply(full_long_ls, function(x) { sum(c(5) %in% x[, 'timepoint']) == 1 })))

sample_ls$infant = Meta_data$baby_externalid[ !is.na(Meta_data$baby_externalid) ] 
Meta_data$externalid
sample_ls$infant_HEU = intersect(sample_ls$infant,
                                 Meta_data$baby_externalid[ Meta_data$hivst[ match(sample_ls$infant, Meta_data$baby_externalid) ] == "Positive" ])
sample_ls$infant_HUU = intersect(sample_ls$infant,
                                 Meta_data$baby_externalid[ Meta_data$hivst[ match(sample_ls$infant, Meta_data$baby_externalid) ] == "Negative" ])
sample_ls$infant = substr(sample_ls$infant, 1, 9) 
sample_ls$infant_HEU = substr(sample_ls$infant_HEU, 1, 9) 
sample_ls$infant_HUU = substr(sample_ls$infant_HUU, 1, 9) 


### Analysis using 71 subjects
colnames(Genus_ct_tab)[which(colnames(Genus_ct_tab) == 'Unknown')] = "Other"
colnames(Genus_abun_tab)[which(colnames(Genus_abun_tab) == 'Unknown')] = "Other"

full_long_ls = Long_func(Sample_keys = Sample_keys, 
                         Variable = 'sample_visit', 
                         Values = c('second trimester enrollment',
                                    'MB enrollment',
                                    'gmb',
                                    'third trimester enrollment',
                                    'third trimester visit'
                         ))
intense_long_ls = full_long_ls[c(sample_ls$`Intense HIV+ pregnant`, 
                                 sample_ls$`Intense HIV- pregnant`)]
intense_pos_ls = full_long_ls[c(sample_ls$`Intense HIV+ pregnant`)] 
intense_neg_ls = full_long_ls[c(sample_ls$`Intense HIV- pregnant`)] 
unlist_mat = do.call(rbind, intense_long_ls) 
rownames(unlist_mat) = unlist_mat$sample_id 

BC_dist = UniFrac_BC(Genus_abun_tab[rownames(unlist_mat), ]) 
Current_Genus_tab = Genus_abun_tab[rownames(unlist_mat), ] 

### Figure 1
set.seed(123)
Stab_data = Current_Genus_tab 
dim(Stab_data) 

Stab_data = remove_col_0_counts(Stab_data) 
dim(Stab_data) 

temp_Genera_vec = colnames(Stab_data) 
temp_Genera_vec = gsub("/", "_", temp_Genera_vec)

colnames(Stab_data) = temp_Genera_vec 
Stab_data = Stab_data[, sort(colnames(Stab_data))] 
Other_id = which(colnames(Stab_data) == 'Other') 
Stab_data = cbind(Stab_data[, - Other_id], Stab_data[, Other_id, drop = F]) 
Genera_vec = colnames(Stab_data) 

Alpha_div_vec = c("Shannon", "Simpson", "Chao1", "Fisher") 
Stab_data = cbind(Stab_data, 
                  unlist_mat[rownames(Stab_data), Alpha_div_vec]) 
Stab_data$Sample_ID = rownames(Stab_data) 
Stab_data$Subject_ID = Sample_keys$first_half_id[ match(Stab_data$Sample_ID, Sample_keys$sample_id) ] 

Quant_vec = c(Alpha_div_vec) 

Stab_data$specimen_GA = Sample_keys$specimen_GA[ match(Stab_data$Sample_ID, Sample_keys$sample_id) ] 
Stab_data$specimen_GA_rounded = round(Stab_data$specimen_GA) 

Stab_data = as.data.frame(Stab_data)
Stab_data$HIV_pos = NA
Stab_data$HIV_pos[ which(Stab_data$Subject_ID %in% sample_ls$"HIV+") ] = 1
Stab_data$HIV_pos[ which(Stab_data$Subject_ID %in% sample_ls$"HIV-") ] = 0

model_result_mat_ls = list() 
for (model_id in 1:4) { 
  
  coef_2_save_vec = c('Intercept',
                      't_coef', 't_pvalue', 
                      'HIV_coef', 'HIV_pvalue',
                      'HIVxt_coef', 'HIVxt_pvalue') 
  model_result_mat = matrix(NA, 
                            nrow = length(Quant_vec), 
                            ncol = length(coef_2_save_vec)) 
  rownames(model_result_mat) = Quant_vec
  colnames(model_result_mat) = coef_2_save_vec
  
  for (q_id in 1:length(Quant_vec)) { 
    current_quant = Quant_vec[q_id] 
    
    if (model_id == 1) { 
      FORMULA = as.formula(paste(current_quant, " ~ specimen_GA * HIV_pos + (1 | Subject_ID)"))
      LMM_data = Stab_data 
    } else if (model_id == 2) { 
      FORMULA = as.formula(paste(current_quant, " ~ specimen_GA + HIV_pos + (1 | Subject_ID)"))
      LMM_data = Stab_data 
    } else if (model_id == 3) { 
      FORMULA = as.formula(paste(current_quant, " ~ specimen_GA + (1 | Subject_ID)"))
      LMM_data = Stab_data[which(Stab_data$HIV_pos == 1), ] 
    } else if (model_id == 4) { 
      FORMULA = as.formula(paste(current_quant, " ~ specimen_GA + (1 | Subject_ID)"))
      LMM_data = Stab_data[which(Stab_data$HIV_pos == 0), ] 
    }
    LMM_model = lmerTest::lmer(FORMULA, 
                               LMM_data, 
                               REML = T) 
    LMM_model_summary = summary(LMM_model) 
    
    
    model_result_mat[q_id, c('Intercept', 't_coef', 't_pvalue')] =
      c(LMM_model_summary$coefficients['(Intercept)', 'Estimate'],
        LMM_model_summary$coefficients['specimen_GA', 'Estimate'],
        LMM_model_summary$coefficients['specimen_GA', 'Pr(>|t|)'])
    if (model_id %in% c(1,2)) { 
      model_result_mat[q_id, c('HIV_coef', 'HIV_pvalue')] =
        c(LMM_model_summary$coefficients['HIV_pos', 'Estimate'],
          LMM_model_summary$coefficients['HIV_pos', 'Pr(>|t|)'])
    }
    if (model_id == 1) { 
      model_result_mat[q_id, c('HIVxt_coef', 'HIVxt_pvalue')] =
        c(LMM_model_summary$coefficients['specimen_GA:HIV_pos', 'Estimate'],
          LMM_model_summary$coefficients['specimen_GA:HIV_pos', 'Pr(>|t|)'])
    }
  } 
  model_result_mat_ls[[ model_id ]] = model_result_mat
} 
names(model_result_mat_ls) = c('LMM_with_interaction',
                               'LMM_without_interaction',
                               'LMM_without_interaction_HIVpos_only',
                               'LMM_without_interaction_HIVneg_only')
### For sentence "The coefficient of gestational age within SN is 0.005 (p=0.30) and within WHIV is <0.001 (p=0.93)"
model_result_mat_ls$LMM_without_interaction_HIVpos_only
model_result_mat_ls$LMM_without_interaction_HIVneg_only


Stab_data_ls = split(Stab_data, 
                     Stab_data$Subject_ID)
par(mfrow = c(1, 1))
for (q_id in 1) { 
  current_quant = Quant_vec[q_id] 
  
  plot(x = Stab_data_ls[[1]]$specimen_GA, 
       y = Stab_data_ls[[1]][, current_quant],
       pch = 16,
       cex = 0.6,
       xlab = 'gestational age',
       ylab = current_quant,
       cex.lab = 1.35,
       ylim = c(min(Stab_data[, current_quant]) * 0.8,
                max(Stab_data[, current_quant]) * 1.2),
       xlim = c(min(Stab_data$specimen_GA) * 0.98,
                max(Stab_data$specimen_GA) * 1.02),
       type = 'o',
       col = ifelse(Stab_data_ls[[1]]$HIV_pos[1] == 0, 'dodgerblue', 'tomato'),
       main = paste0('coefficient of HIV positive status = ', 
                     dec(model_result_mat_ls$LMM_without_interaction[q_id, 'HIV_coef'], 3), ' (',
                     pvalue_dec(model_result_mat_ls$LMM_without_interaction[q_id, 'HIV_pvalue']), ')')
  )
  for (sub_id in 2:length(Stab_data_ls)) {
    lines(x = Stab_data_ls[[ sub_id ]]$specimen_GA, 
          y = Stab_data_ls[[ sub_id ]][, current_quant],
          pch = 16,
          cex = 0.6,
          cex.lab = 1.35,
          type = 'o',
          col = ifelse(Stab_data_ls[[ sub_id ]]$HIV_pos[1] == 0, 'dodgerblue', 'tomato'))
  } 
  
  
  t_range = Stab_data$specimen_GA 
  new_t = seq(min(t_range), max(t_range), length.out = 101)
  lines(x = new_t, 
        y = model_result_mat_ls$LMM_without_interaction[q_id, 'Intercept'] +
          new_t * model_result_mat_ls$LMM_without_interaction[q_id, 't_coef'],
        col = 'blue', lwd = 4)
  lines(x = new_t, 
        y = model_result_mat_ls$LMM_without_interaction[q_id, 'Intercept'] +
          new_t * model_result_mat_ls$LMM_without_interaction[q_id, 't_coef'] +
          model_result_mat_ls$LMM_without_interaction[q_id, 'HIV_coef'],
        col = 'red', lwd = 4)
  
  
  if (q_id == 1) {
    HIV_status_vec = sapply(Stab_data_ls, function(x) x[1, 'HIV_pos']) 
    legend('topleft',
           legend = paste0(c('WHIV', 'SN'), ' n=', c(sum(HIV_status_vec == 1), sum(HIV_status_vec == 0))),
           pch = 16,
           lty = 1,
           col = c('tomato', 'dodgerblue'),
           bty = 'n',
           cex = 1.5)
  }
} 



current_strim_ids = which(Stab_data$specimen_GA < 28) 
current_ttrim_ids = which(Stab_data$specimen_GA >= 28) 
selected_cols = c('Subject_ID', 'Shannon', 'Simpson', 'Chao1', 'Fisher', 'specimen_GA', 'HIV_pos')
Stab_data_2T = Stab_data[current_strim_ids, selected_cols] 
Stab_data_3T = Stab_data[current_ttrim_ids, selected_cols] 


Stab_data_2T_mean = aggregate(. ~ Subject_ID, data = Stab_data_2T, FUN = mean) 
Stab_data_3T_mean = aggregate(. ~ Subject_ID, data = Stab_data_3T, FUN = mean) 
rownames(Stab_data_2T_mean) = Stab_data_2T_mean$Subject_ID
rownames(Stab_data_3T_mean) = Stab_data_3T_mean$Subject_ID

Stab_data_3T_mean = Stab_data_3T_mean[rownames(Stab_data_2T_mean), ]
Stab_data_3T_mean$Subject_ID = rownames(Stab_data_3T_mean) = rownames(Stab_data_2T_mean) 
Stab_data_3T_mean$HIV_pos = Stab_data_2T_mean$HIV_pos 

### Supplementary Table 4
for (group_id in 1:3) { 
  if (group_id == 1) {
    current_samples = rownames(Stab_data_2T_mean) 
    cat('Use all samples:\n')
  } else if (group_id == 2) {
    current_samples = rownames(Stab_data_2T_mean)[ which(Stab_data_2T_mean$HIV_pos == 1) ] 
    cat('Use HIV+ samples:\n')
  } else if (group_id == 3) {
    current_samples = rownames(Stab_data_2T_mean)[ which(Stab_data_2T_mean$HIV_pos == 0) ] 
    cat('Use HIV- samples:\n')
  }
  
  paired_test_mat = matrix(NA, nrow = 1, ncol = 5)
  rownames(paired_test_mat) = c('Shannon')
  
  colnames(paired_test_mat) = c('Sample size', '2T mean', '3T mean', 'test statistic', 'pvalue')
  
  
  for (row_id in 1:nrow(paired_test_mat)) { 
    current_alpha = rownames(paired_test_mat)[ row_id ] 
    values_2T = Stab_data_2T_mean[current_samples, current_alpha] 
    values_3T = Stab_data_3T_mean[current_samples, current_alpha] 
    current_test = t.test(values_2T,
                          values_3T,
                          paired = TRUE)
    paired_test_mat[row_id, c('Sample size', '2T mean', '3T mean', 'test statistic', 'pvalue')] =
      c(sum(!is.na(values_3T)),
        mean(values_2T, na.rm = T),
        mean(values_3T, na.rm = T),
        current_test$statistic,
        current_test$p.value)
    
  } 
  
  print(paired_test_mat)
} 





### Load raw data
colnames(Genus_ct_tab)[which(colnames(Genus_ct_tab) == 'Unknown')] = "Other"
colnames(Genus_abun_tab)[which(colnames(Genus_abun_tab) == 'Unknown')] = "Other"


full_long_ls = Long_func(Sample_keys = Sample_keys, 
                         Variable = 'sample_visit', 
                         Values = c('second trimester enrollment',
                                    'MB enrollment',
                                    'gmb',
                                    'third trimester enrollment',
                                    'third trimester visit',
                                    'post partum',
                                    'infant'))


for (sub_id in 1:length(full_long_ls)) { 
  current_data = full_long_ls[[ sub_id ]] 
  temp_matched_ids = match(current_data$sample_id, Sample_keys$sample_id) 
  current_data$Subject_ID = Sample_keys$first_half_id[ temp_matched_ids ] 
  
  
  pp_row_ids = which(current_data$sample_visit == 'post partum') 
  if (length(pp_row_ids) > 0) {
    current_data[pp_row_ids, 'specimen_GA'] = 
      Meta_data[current_data$Subject_ID[1], 'ges_delivery'] + 
      Sample_keys[temp_matched_ids[ pp_row_ids ], 'specimen_postpartum_mother'] 
  }
  current_data = current_data[order(current_data$specimen_GA), ] 
  
  
  current_data$HIV_pos = NA
  current_data$HIV_pos[ which(current_data$Subject_ID %in% c(sample_ls$"HIV+",
                                                             sample_ls$infant_HEU)) ] = 1
  current_data$HIV_pos[ which(current_data$Subject_ID %in% c(sample_ls$"HIV-",
                                                             sample_ls$infant_HUU)) ] = 0
  
  
  temp_ids_meta_data = match(current_data$first_half_id, Meta_data$first_half_id) 
  
  current_data$age = Meta_data$age[ temp_ids_meta_data ] 
  
  current_data$post_high_school = Meta_data$post_high_school[ temp_ids_meta_data ] 
  
  current_strim_ids = which(current_data$specimen_GA < 28) 
  current_ttrim_ids = which(current_data$specimen_GA >= 28 & (! current_data$sample_visit %in% c('post partum', 'infant'))) 
  current_6mopp_ids = which(current_data$sample_visit == 'post partum') 
  current_infant_ids = which(current_data$sample_visit == 'infant') 
  current_data$sample_visit_2T3TPP = NA 
  current_data$sample_visit_2T3TPP[ current_strim_ids ] = '2T'
  current_data$sample_visit_2T3TPP[ current_ttrim_ids ] = '3T'
  current_data$sample_visit_2T3TPP[ current_6mopp_ids ] = 'PP'
  current_data$sample_visit_2T3TPP[ current_infant_ids ] = 'infant'
  
  current_data$muac = NA
  if (length(current_strim_ids) > 0) { current_data[current_strim_ids, 'muac'] = Meta_data$muacstrim[ temp_ids_meta_data[ current_strim_ids ] ] }
  if (length(current_ttrim_ids) > 0) { current_data[current_ttrim_ids, 'muac'] = Meta_data$muacttrim[ temp_ids_meta_data[ current_ttrim_ids ] ] }
  if (length(current_6mopp_ids) > 0) { current_data[current_6mopp_ids, 'muac'] = Meta_data$muac_6mo_pp[ temp_ids_meta_data[ current_6mopp_ids ] ] }
  
  
  if (length(current_infant_ids) > 0) { current_data[current_infant_ids, 'muac'] = Meta_data$bmi_zscore_6mo[ temp_ids_meta_data[ current_infant_ids ] ] }
  
  
  full_long_ls[[ sub_id ]] = current_data
} 
Stab_data = do.call(rbind, full_long_ls)
rownames(Stab_data) = Stab_data$sample_id

### Supplementary Figure 1
Stab_data_2T_HIV_pos = Stab_data[which(Stab_data$sample_visit_2T3TPP == '2T' & 
                                         Stab_data$HIV_pos == 1), ]
Stab_data_2T_HIV_neg = Stab_data[which(Stab_data$sample_visit_2T3TPP == '2T' & 
                                         Stab_data$HIV_pos == 0), ]
Stab_data_3T_HIV_pos = Stab_data[which(Stab_data$sample_visit_2T3TPP == '3T' & 
                                         Stab_data$HIV_pos == 1), ]
Stab_data_3T_HIV_neg = Stab_data[which(Stab_data$sample_visit_2T3TPP == '3T' & 
                                         Stab_data$HIV_pos == 0), ]
Stab_data_PP_HIV_pos = Stab_data[which(Stab_data$sample_visit_2T3TPP == 'PP' & 
                                         Stab_data$HIV_pos == 1), ]
Stab_data_PP_HIV_neg = Stab_data[which(Stab_data$sample_visit_2T3TPP == 'PP' & 
                                         Stab_data$HIV_pos == 0), ]
Quant_vec = c("Shannon", "Simpson", "Chao1", "Fisher")
for (q_id in 1) { 
  current_quant = Quant_vec[q_id] 
  
  boxplot(list(Stab_data_2T_HIV_pos[, current_quant], 
               Stab_data_2T_HIV_neg[, current_quant], 
               Stab_data_3T_HIV_pos[, current_quant], 
               Stab_data_3T_HIV_neg[, current_quant], 
               Stab_data_PP_HIV_pos[, current_quant], 
               Stab_data_PP_HIV_neg[, current_quant]),
          at = c(1,2, 4,5, 7,8),
          names = c("2T WHIV", "2T SN", 
                    "3T WHIV", "3T SN", 
                    "PP WHIV", "PP SN"),
          col = c("pink2", "skyblue"),
          main = current_quant,
          xlab = "",
          ylab = "",
          pch = 16,
          cex = 0.4)
  
  points(c(1,2, 4,5, 7,8),
         c(mean(Stab_data_2T_HIV_pos[, current_quant]),
           mean(Stab_data_2T_HIV_neg[, current_quant]),
           mean(Stab_data_3T_HIV_pos[, current_quant]),
           mean(Stab_data_3T_HIV_neg[, current_quant]),
           mean(Stab_data_PP_HIV_pos[, current_quant]),
           mean(Stab_data_PP_HIV_neg[, current_quant])
         ),
         col = c('red', 'blue'),
         pch = 15)
  
}


library(ANCOMBC)
ANCOMBC_meta_ls = list() 
ANCOMBC_data_ls = list() 
ANCOMBC_fit_ls = list() 
ANCOMBC_result_mat_ls = list() 
for (set_id in 1:3) { 
  set.seed(set_id)
  if (set_id == 1) { 
    CURRENT_Stab_data = Stab_data[which(Stab_data$sample_visit_2T3TPP %in% c('2T', '3T')), ]
    rand_formula = "(1 | Subject_ID)"
  } else if (set_id == 2) { 
    CURRENT_Stab_data = Stab_data[which(Stab_data$sample_visit_2T3TPP %in% c('PP')), ]
    rand_formula = NULL
  } else if (set_id == 3) { 
    CURRENT_Stab_data = Stab_data[which(Stab_data$sample_visit_2T3TPP %in% c('infant')), ]
    rand_formula = NULL
  }
  
  
  ANCOMBC_meta = CURRENT_Stab_data[, c('Subject_ID', 'specimen_GA', 
                                       'HIV_pos', 'age', 'post_high_school', 'muac')] 
  ANCOMBC_meta$HIV_pos = factor(ANCOMBC_meta$HIV_pos,
                                levels = c(0, 1))
  ANCOMBC_data = as.data.frame( t(Genus_ct_tab[rownames(ANCOMBC_meta), ]) ) 
  ANCOMBC_meta_ls[[ set_id ]] = ANCOMBC_meta 
  ANCOMBC_data_ls[[ set_id ]] = ANCOMBC_data 
  
  
  ANCOMBC_fit_ls[[ set_id ]] = list()
  ANCOMBC_result_mat_ls[[ set_id ]] = list()
  for (model_id in 1:2) { 
    if (model_id == 1) {
      fix_formula = "HIV_pos"
    } else if (model_id == 2) {
      fix_formula = "HIV_pos + age + post_high_school + muac"
    }
    output_model0_BH = ancombc2(data = as.data.frame(ANCOMBC_data),
                                meta_data = ANCOMBC_meta,
                                assay_name = "counts", 
                                tax_level = "Genus",
                                fix_formula = fix_formula,
                                rand_formula = rand_formula,
                                group = "HIV_pos", 
                                p_adj_method = "BH", 
                                
                                pseudo = 0.1,
                                pseudo_sens = FALSE,
                                
                                prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                                struc_zero = TRUE, neg_lb = TRUE,
                                alpha = 0.05, n_cl = 1, verbose = TRUE,
                                iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE))
    closeAllConnections()
    ANCOMBC_fit_ls[[ set_id ]][[ model_id ]] = output_model0_BH
    
    
    ANCOMBC_result_mat = output_model0_BH$res 
    ANCOMBC_result_mat = ANCOMBC_result_mat[order(- sign(ANCOMBC_result_mat$lfc_HIV_pos1), 
                                                  (ANCOMBC_result_mat$p_HIV_pos1)), ]
    ANCOMBC_result_mat_ls[[ set_id ]][[ model_id ]] = ANCOMBC_result_mat
  }
  names(ANCOMBC_fit_ls[[ set_id ]]) = names(ANCOMBC_result_mat_ls[[ set_id ]]) = 
    c('ANCOM_without_covariates',
      'ANCOM_with_covariates')
}
names(ANCOMBC_fit_ls) = names(ANCOMBC_result_mat_ls) = names(ANCOMBC_meta_ls) = names(ANCOMBC_data_ls) =
  c('2T3T_samples',
    'PP_samples',
    'Infant_samples')

c_1_96 = qnorm(0.975)
current_ANCOM_results_COMBINED_ls = list()
taxa_structural_zero_ls = list()
for (set_id in 1:3) { 
  current_ANCOM_results = ANCOMBC_result_mat_ls[[ set_id ]]$ANCOM_with_covariates 
  
  current_ANCOM_results = current_ANCOM_results[, c('taxon', 'lfc_HIV_pos1', 'p_HIV_pos1', 'q_HIV_pos1', 'se_HIV_pos1')] 
  
  current_ANCOM_results$lower_lfc_HIV_pos1 = current_ANCOM_results$lfc_HIV_pos1 - c_1_96 * current_ANCOM_results$se_HIV_pos1
  current_ANCOM_results$upper_lfc_HIV_pos1 = current_ANCOM_results$lfc_HIV_pos1 + c_1_96 * current_ANCOM_results$se_HIV_pos1
  
  current_ANCOM_results$significance = NA
  current_ANCOM_results$significance[ which(current_ANCOM_results$q_HIV_pos1 <= 0.05) ] = 'Significant'
  
  
  current_ANCOM_results = current_ANCOM_results[order(- sign(current_ANCOM_results$lfc_HIV_pos1),
                                                      current_ANCOM_results$q_HIV_pos1), ]
  
  current_ANCOM_results_COMBINED_ls[[ set_id ]] = current_ANCOM_results 
  
  
  temp_mat = ANCOMBC_fit_ls[[ set_id ]]$ANCOM_with_covariates$zero_ind 
  tab = table(temp_mat$`structural_zero (HIV_pos = 0)`,
              temp_mat$`structural_zero (HIV_pos = 1)`)
  names(attributes(tab)$dimnames) = c('HIV-', 'HIV+')
  
  taxa_structural_zero = temp_mat[which(temp_mat$`structural_zero (HIV_pos = 0)` + 
                                          temp_mat$`structural_zero (HIV_pos = 1)` == 1), 1:3] 
  taxa_structural_zero = taxa_structural_zero[order(taxa_structural_zero[,2],
                                                    taxa_structural_zero[,3],
                                                    taxa_structural_zero[,1]), ]
  taxa_structural_zero_ls[[ set_id ]] = taxa_structural_zero
  
} 
names(current_ANCOM_results_COMBINED_ls) = 
  names(taxa_structural_zero_ls) =
  names(ANCOMBC_result_mat_ls)

current_ANCOM_results_COMBINED_ls
taxa_structural_zero_ls$PP_samples # Table 2
taxa_structural_zero_ls$Infant_samples # Table 3



### Supplementary Table 7 & Figure 2
temp = current_ANCOM_results_COMBINED_ls$`2T3T_samples`
rownames(temp) = temp$taxon
temp[c('Megamonas',
       'Lachnoclostridium',
       'Lachnospiraceae_NK4A136_group',
       'Fusobacterium'), ]

current_ANCOM_results = ANCOMBC_result_mat_ls$`2T3T_samples`$ANCOM_without_covariates 
current_ANCOM_results = current_ANCOM_results[, c('taxon', 'lfc_HIV_pos1', 'p_HIV_pos1', 'q_HIV_pos1', 'se_HIV_pos1')] 

current_ANCOM_results$lower_lfc_HIV_pos1 = current_ANCOM_results$lfc_HIV_pos1 - c_1_96 * current_ANCOM_results$se_HIV_pos1
current_ANCOM_results$upper_lfc_HIV_pos1 = current_ANCOM_results$lfc_HIV_pos1 + c_1_96 * current_ANCOM_results$se_HIV_pos1

rownames(current_ANCOM_results) = current_ANCOM_results$taxon
temp_out = current_ANCOM_results[c('Megamonas',
                                   'Lachnoclostridium',
                                   'Lachnospiraceae_NK4A136_group',
                                   'Fusobacterium'), 
                                 c('taxon', 'lfc_HIV_pos1', 'lower_lfc_HIV_pos1', 'upper_lfc_HIV_pos1', 'p_HIV_pos1')]
temp_out



PERMANOVA_pvalue_mat = matrix(NA, nrow = 2, ncol = 4)
rownames(PERMANOVA_pvalue_mat) = c('Model 0', 'Model 1')
colnames(PERMANOVA_pvalue_mat) = c('3T', '2T', 'PP', 'Infant')
for (set_id in 1:4) { 
  if (set_id == 1) { 
    CURRENT_Stab_data = Stab_data[which(Stab_data$sample_visit %in% c('third trimester enrollment', 
                                                                      'third trimester visit')), ]
    CURRENT_Stab_data[which(CURRENT_Stab_data$first_half_id == 'PRA-11463'), ]
    CURRENT_Stab_data = CURRENT_Stab_data[- which(CURRENT_Stab_data$first_half_id == 'PRA-11463' &
                                                    CURRENT_Stab_data$specimen_GA == 31.7), ]
  } else if (set_id == 2) { 
    CURRENT_Stab_data = Stab_data[which(Stab_data$sample_visit %in% c('second trimester enrollment')), ]
  } else if (set_id == 3) { 
    CURRENT_Stab_data = Stab_data[which(Stab_data$sample_visit %in% c('post partum')), ]
  } else if (set_id == 4) { 
    CURRENT_Stab_data = Stab_data[which(Stab_data$sample_visit %in% c('infant')), ]
  }
  if (F) table(sapply(split(CURRENT_Stab_data, CURRENT_Stab_data$Subject_ID), nrow)) 
  print(table(CURRENT_Stab_data$HIV_pos))
  n_distinct(CURRENT_Stab_data$first_half_id)
  
  colnames(PERMANOVA_pvalue_mat)[ set_id ] = paste0(colnames(PERMANOVA_pvalue_mat)[ set_id ],
                                                    ' n=', nrow(CURRENT_Stab_data))
  
  
  PERMANOVA_meta = CURRENT_Stab_data[, c('Subject_ID', 'specimen_GA', 
                                         'HIV_pos', 'age', 'post_high_school', 'muac')] 
  PERMANOVA_meta$HIV_pos = factor(PERMANOVA_meta$HIV_pos,
                                  levels = c(0, 1))
  PERMANOVA_data = as.data.frame(Genus_ct_tab[rownames(PERMANOVA_meta), ] ) 
  PERMANOVA_data = PERMANOVA_data / rowSums(PERMANOVA_data) 
  
  current_BC_dist = UniFrac_BC(PERMANOVA_data)
  
  library(vegan) 
  n_perm = 9999 
  Group = CURRENT_Stab_data$HIV_pos
  Age = CURRENT_Stab_data$age
  Post_high_school = CURRENT_Stab_data$post_high_school
  Muac = CURRENT_Stab_data$muac
  
  set.seed(set_id)
  PM_results = vegan::adonis2(current_BC_dist ~ Group, permutations = n_perm, by = "terms")
  PM_pvalue = PM_results$`Pr(>F)`[1]
  
  set.seed(set_id)
  PM_results_model1 = vegan::adonis2(current_BC_dist ~ Group + Age + Post_high_school + Muac, 
                                     permutations = n_perm, na.action = na.omit, by = "terms")
  PM_pvalue_model1 = PM_results_model1$`Pr(>F)`[1]
  
  PERMANOVA_pvalue_mat[, set_id] = c(PM_pvalue, PM_pvalue_model1)
} 

noquote(dec(PERMANOVA_pvalue_mat, 3)) # Supplementary Table 5


for (set_id in 3:4) { 
  if (set_id == 1) { 
    CURRENT_Stab_data = Stab_data[which(Stab_data$sample_visit %in% c('third trimester enrollment', 
                                                                      'third trimester visit')), ]
    CURRENT_Stab_data[which(CURRENT_Stab_data$first_half_id == 'PRA-11463'), ]
    CURRENT_Stab_data = CURRENT_Stab_data[- which(CURRENT_Stab_data$first_half_id == 'PRA-11463' &
                                                    CURRENT_Stab_data$specimen_GA == 31.7), ]
  } else if (set_id == 2) { 
    CURRENT_Stab_data = Stab_data[which(Stab_data$sample_visit %in% c('second trimester enrollment')), ]
  } else if (set_id == 3) { 
    CURRENT_Stab_data = Stab_data[which(Stab_data$sample_visit %in% c('post partum')), ]
  } else if (set_id == 4) { 
    CURRENT_Stab_data = Stab_data[which(Stab_data$sample_visit %in% c('infant')), ]
  }
  
  
  alpha_div_pvalue_mat = matrix(NA, nrow = 1, ncol = 4)
  rownames(alpha_div_pvalue_mat) = 'Shannon'
  colnames(alpha_div_pvalue_mat) = c('Model0', 'Model1', 'n_Model0', 'n_Model1')
  for (alpha_id in 1) { 
    for (model_id in 1:2) {
      if (model_id == 1) {
        
        test_results = t.test(CURRENT_Stab_data[which(CURRENT_Stab_data$HIV_pos == 1), 'Shannon'],
                              CURRENT_Stab_data[which(CURRENT_Stab_data$HIV_pos == 0), 'Shannon'])
        alpha_div_pvalue_mat[alpha_id, 'Model0'] = test_results$p.value
        alpha_div_pvalue_mat['Shannon', paste0('n_Model', model_id-1)] = nrow(CURRENT_Stab_data)
      } else if (model_id == 2) { 
        glm_fit = glm(Shannon ~ HIV_pos + age + post_high_school + muac, 
                      data = CURRENT_Stab_data,
                      family = 'gaussian')
        coef_mat = summary(glm_fit)$coefficients[, c(1,4)]
        colnames(coef_mat)[which(colnames(coef_mat) == 'Pr(>|t|)')] = 'pvalue'
        alpha_div_pvalue_mat[alpha_id, 'Model1'] = coef_mat['HIV_pos', 'pvalue']
        alpha_div_pvalue_mat['Shannon', paste0('n_Model', model_id-1)] = nobs(glm_fit)
      } 
      
      
    }
  }
  alpha_div_pvalue_mat
  
  
  
  par(mfrow = c(1, 1), las = 0, mar = c(5, 4, 3, 2) + 0.1, cex.main = 0.8)
  main = paste0("Shannon", "\n", 
                'Model 0: ', pvalue_dec(alpha_div_pvalue_mat[alpha_id, 'Model0'], threshold = 0.001, k = 3), 
                " (n=", alpha_div_pvalue_mat[alpha_id, 'n_Model0'], ")\n",
                'Model 1: ', pvalue_dec(alpha_div_pvalue_mat[alpha_id, 'Model1'], threshold = 0.001, k = 3), 
                " (n=", alpha_div_pvalue_mat[alpha_id, 'n_Model1'], ")")
  CURRENT_HIV_0_ids = which(CURRENT_Stab_data$HIV_pos == 0)
  CURRENT_HIV_1_ids = which(CURRENT_Stab_data$HIV_pos == 1)
  
  
  
  
  
  
  
  
  bb = boxplot(Shannon ~ HIV_pos, 
               data = CURRENT_Stab_data,
               col = "white",
               main = main,
               ylab = '',
               xlab = '',
               boxwex = 0.5,
               pch = 16, 
               cex = 0.4,
               plot = F
               
               
  )
  temp_ids = match(bb$names, c('0', '1'))
  if (set_id == 3) {
    bb$names = c(paste0('SN n=', length(CURRENT_HIV_0_ids)),
                 paste0('WHIV n=', length(CURRENT_HIV_1_ids)))[ temp_ids ]
  } else if (set_id == 4) {
    bb$names = c(paste0('CHUU n=', length(CURRENT_HIV_0_ids)),
                 paste0('CHEU n=', length(CURRENT_HIV_1_ids)))[ temp_ids ]
  }
  bxp(bb, main = main,
      boxwex = 0.5,
      pch = 16, 
      cex = 0.4
  )
  points(c(1, 2), 
         c(mean(CURRENT_Stab_data[CURRENT_HIV_0_ids, 'Shannon']),
           mean(CURRENT_Stab_data[CURRENT_HIV_1_ids, 'Shannon']))[ temp_ids ], 
         pch = 15, col = 'red') 
} 


### Supplementary Table 6
for (set_id in 1:3) { 
  set.seed(set_id)
  
  
  model_result_mat_ls = list() 
  for (model_id in 2) { 
    
    coef_2_save_vec = c('Intercept',
                        't_coef', 't_pvalue', 
                        'HIV_coef', 'HIV_pvalue',
                        'HIVxt_coef', 'HIVxt_pvalue') 
    model_result_mat = matrix(NA, 
                              nrow = length(Quant_vec), 
                              ncol = length(coef_2_save_vec)) 
    rownames(model_result_mat) = Quant_vec
    colnames(model_result_mat) = coef_2_save_vec
    
    
    for (q_id in 1:1) { 
      current_quant = Quant_vec[q_id] 
      
      FORMULA = as.formula(paste(current_quant, " ~ HIV_pos + age + post_high_school + muac + (1 | Subject_ID)"))
      LMM_data = Stab_data 
      if (set_id == 1) {
        LMM_data = LMM_data[which(LMM_data$sample_visit_2T3TPP %in% c('2T', '3T')), ]
      } else if (set_id == 2) {
        LMM_data = LMM_data[which(LMM_data$sample_visit_2T3TPP %in% c('2T')), ]
      } else if (set_id == 3) {
        LMM_data = LMM_data[which(LMM_data$sample_visit_2T3TPP %in% c('3T')), ]
      }
      LMM_model = lmerTest::lmer(FORMULA, 
                                 LMM_data, 
                                 REML = T) 
      LMM_model_summary = summary(LMM_model) 
      
      
      if (model_id %in% c(1,2)) { 
        model_result_mat[q_id, c('HIV_coef', 'HIV_pvalue')] =
          c(LMM_model_summary$coefficients['HIV_pos', 'Estimate'],
            LMM_model_summary$coefficients['HIV_pos', 'Pr(>|t|)'])
      }
      if (model_id == 1) { 
        model_result_mat[q_id, c('HIVxt_coef', 'HIVxt_pvalue')] =
          c(LMM_model_summary$coefficients['specimen_GA:HIV_pos', 'Estimate'],
            LMM_model_summary$coefficients['specimen_GA:HIV_pos', 'Pr(>|t|)'])
      }
    } 
    model_result_mat_ls[[ model_id ]] = model_result_mat
  } 
  names(model_result_mat_ls)[2] = 'LMM_without_interaction'
  print(model_result_mat_ls$LMM_without_interaction[1, , drop = F])
  
} 


### Load raw data
colnames(Genus_ct_tab)[which(colnames(Genus_ct_tab) == 'Unknown')] = "Other"
colnames(Genus_abun_tab)[which(colnames(Genus_abun_tab) == 'Unknown')] = "Other"


full_long_ls = Long_func(Sample_keys = Sample_keys, 
                         Variable = 'sample_visit', 
                         Values = c('second trimester enrollment',
                                    'MB enrollment',
                                    'gmb',
                                    'third trimester enrollment',
                                    'third trimester visit'
                                    
                         ))
intense_long_ls = full_long_ls[c(sample_ls$`Intense HIV+ pregnant`, 
                                 sample_ls$`Intense HIV- pregnant`)]
intense_pos_ls = full_long_ls[c(sample_ls$`Intense HIV+ pregnant`)] 
intense_neg_ls = full_long_ls[c(sample_ls$`Intense HIV- pregnant`)] 
unlist_mat = do.call(rbind, intense_long_ls) 
rownames(unlist_mat) = unlist_mat$sample_id 



library(readxl)

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

raw_data = read_excel_allsheets('./Raw_data/COLU-05-19VW Data Tables.xlsx')
raw_data_batch = raw_data$`Batch-normalized Data`
sample_names = substr(raw_data$`Sample Meta Data`$CLIENT_SAMPLE_ID, 1, 9) 
rownames(raw_data_batch) = sample_names
raw_data_batch = raw_data_batch[, -1]
raw_data_ls = list(batch = raw_data_batch)

meta_impute_new = raw_data_ls$batch 
for (chem_id in 1:ncol(meta_impute_new)) { 
  old_values = raw_data_ls$batch[, chem_id] 
  meta_impute_new[is.na(old_values), chem_id] = min(old_values, na.rm = T) / 2 
}
meta_final = log(meta_impute_new) 


sample_data = raw_data$`Sample Meta Data` 
rownames(sample_data) = sample_data$PARENT_SAMPLE_NAME 
chem_annotation = raw_data$`Chemical Annotation` 
rownames(chem_annotation) = paste0('CHEM_', chem_annotation$CHEM_ID) 
colnames(meta_final) = paste0('CHEM_', colnames(meta_final)) 


chem_annotation$missing_rate = colMeans(is.na(raw_data_ls$batch)) 
chem_annotation$"n_nonmissing" = NA
chem_annotation$"n_HIV+" = NA
chem_annotation$"n_HIV-" = NA
for (d_id in 1:nrow(chem_annotation)) { 
  temp_samples = rownames(raw_data_ls$batch)[ which(!is.na(raw_data_ls$batch[, d_id])) ] 
  temp_HIV_status = Meta_data[temp_samples, 'hivst'] 
  chem_annotation[d_id, 'n_nonmissing'] = length(temp_samples) 
  chem_annotation[d_id, 'n_HIV+'] = sum(temp_HIV_status == 'Positive')
  chem_annotation[d_id, 'n_HIV-'] = sum(temp_HIV_status == 'Negative')
}

chem_to_keep = rownames(chem_annotation)[ which(chem_annotation$missing_rate <= 0.7) ] 
meta_final_full = meta_final 
meta_final = meta_final[, chem_to_keep] 
# meta_final_METABOLON = meta_final_METABOLON[, chem_to_keep]

META_pat_vec = rownames(meta_final) 
n_sam_vec = unlist(sapply(full_long_ls[ META_pat_vec ], nrow)) 
full_long_3T_ls = sapply(full_long_ls[ META_pat_vec ], 
                         function(x) {x[which(x$sample_visit %in% c('third trimester enrollment',
                                                                    'third trimester visit')), , drop = F]})
full_long_3T = do.call(rbind, full_long_3T_ls) 

Meta_data_new_var = Meta_data[, which(!colnames(Meta_data) %in% colnames(full_long_3T))] 
full_long_3T = cbind(full_long_3T, 
                     Meta_data_new_var[rownames(full_long_3T), ])


METAB_data = meta_final[rownames(full_long_3T), ] 
MICRO_data = Genus_abun_tab[full_long_3T$sample_id, ] 
rownames(MICRO_data) = full_long_3T$first_half_id 
MICRO_data = remove_col_0_counts(MICRO_data) 
MM_data = cbind(METAB_data, 
                MICRO_data)
n_metab = ncol(METAB_data)


cor_mat = cor(MM_data)
idx1 = 1:n_metab
idx2 = (n_metab + 1):ncol(MM_data)
clust1 = hclust(dist(cor_mat[idx2, idx1]))
sub_cor_mat = cor_mat[idx2, idx1] 
row_dist = as.dist(1 - cor(t(sub_cor_mat)))   
row_clust = hclust(row_dist)
col_dist = as.dist(1 - cor(sub_cor_mat))      
col_clust = hclust(col_dist)
sub_cor_mat_ordered = sub_cor_mat[row_clust$order, col_clust$order]


pos_id = as.data.frame(which(sub_cor_mat_ordered >= 0.6, arr.ind = TRUE))
pos_id$correlation = sub_cor_mat_ordered[ as.matrix(pos_id) ]
pos_id = pos_id[order(- pos_id$correlation), ]
pos_id[, 'row'] = noquote(rownames(sub_cor_mat_ordered)[pos_id[, 'row']])
pos_id[, 'col'] = noquote(chem_annotation[colnames(sub_cor_mat_ordered)[pos_id[, 'col']], 'CHEMICAL_NAME'])
rownames(pos_id) = NULL
colnames(pos_id) = c('genus', 'chemical', 'correlation')


neg_id = as.data.frame(which(sub_cor_mat_ordered <= - 0.6, arr.ind = TRUE))
neg_id$correlation = sub_cor_mat_ordered[ as.matrix(neg_id) ]
neg_id = neg_id[order(neg_id$correlation), ]
neg_id[, 'row'] = noquote(rownames(sub_cor_mat_ordered)[neg_id[, 'row']])
neg_id[, 'col'] = noquote(chem_annotation[colnames(sub_cor_mat_ordered)[neg_id[, 'col']], 'CHEMICAL_NAME'])
rownames(neg_id) = NULL
colnames(neg_id) = c('genus', 'chemical', 'correlation')



Cases = full_long_3T$first_half_id[ which(full_long_3T$hivst == 'Positive') ]
Controls = full_long_3T$first_half_id[ which(full_long_3T$hivst == 'Negative') ]
Case_Group_Name = 'HIV+'
Control_Group_Name = 'HIV-'
Case_Group_Name_with_n = paste0(Case_Group_Name, ' (n=', length(Cases), ')')
Control_Group_Name_with_n = paste0(Control_Group_Name, ' (n=', length(Controls), ')')
quant_vec = c(paste0('mean ', Case_Group_Name),
              paste0('sd ', Case_Group_Name),
              paste0('mean ', Control_Group_Name),
              paste0('sd ', Control_Group_Name),
              paste0('mean ', Case_Group_Name, ' - mean ', Control_Group_Name),
              'test statistic',
              'raw pvalue', 
              'BF adjusted pvalue', 
              'BH adjusted pvalue')

pvalue_mat_ls = list() 
ancom_model_fit_ls = list() 
for (comp_id in 1:5) { 
  cat(comp_id)
  
  
  
  
  
  if (comp_id == 1) {
    CURRENT_DATA = METAB_data
    n_feature = ncol(CURRENT_DATA)
    pvalue_mat = matrix(NA, nrow = n_feature, ncol = length(quant_vec))
    rownames(pvalue_mat) = colnames(CURRENT_DATA)
    colnames(pvalue_mat) = quant_vec
    
    
    for (feature_id in 1:n_feature) {
      case_values = as.vector(CURRENT_DATA[Cases, feature_id])
      control_values = as.vector(CURRENT_DATA[Controls, feature_id])
      if (length(unique(c(case_values, control_values))) == 1) { 
        test_stat = 0; raw_pvalue = 1
      } else {
        if (min(length(unique(case_values)),
                length(unique(case_values))) == 1) 
        {
          test_results = t.test(case_values, control_values)
        } else {
          test_results = t.test(case_values, control_values, var.equal = TRUE)
        }
        test_stat = test_results$statistic
        raw_pvalue = test_results$p.value
      }
      pvalue_mat[feature_id, quant_vec] =
        c(mean(case_values),
          sd(case_values),
          mean(control_values),
          sd(control_values),
          mean(case_values) - mean(control_values),
          test_stat,
          raw_pvalue,
          NA, NA)
    }
    pvalue_mat[, 'BF adjusted pvalue'] = p.adjust(pvalue_mat[, 'raw pvalue'], method = "bonferroni")
    pvalue_mat[, 'BH adjusted pvalue'] = p.adjust(pvalue_mat[, 'raw pvalue'], method = "BH")
    pvalue_mat = as.data.frame(pvalue_mat)
    pvalue_mat$"+/-" = NA
    pvalue_mat$"+/-"[ which(pvalue_mat$`test statistic` > 0) ] = paste0(Case_Group_Name, ' > ', Control_Group_Name) 
    pvalue_mat$"+/-"[ which(pvalue_mat$`test statistic` == 0) ] = paste0(Case_Group_Name, ' = ', Control_Group_Name)
    pvalue_mat$"+/-"[ which(pvalue_mat$`test statistic` < 0) ] = paste0(Case_Group_Name, ' < ', Control_Group_Name)
    pvalue_mat$significance = ''
    pvalue_mat$significance[ pvalue_mat$`BH adjusted pvalue` <= 0.05 ] = 'Significant'
    
    
    
    pvalue_mat = pvalue_mat[order(pvalue_mat$`+/-`,
                                  pvalue_mat$`raw pvalue`), ]
    
    
  } else if (comp_id >= 2) { 
    
    ancom_model_fit = NULL
    
    library(ANCOMBC)
    CURRENT_DATA = as.data.frame( t(Genus_ct_tab[full_long_3T$sample_id, ]) ) 
    colnames(CURRENT_DATA) = full_long_3T$first_half_id 
    full_long_3T$hivst
    
    if (comp_id == 2) { 
      ancom_model_fit = ancombc2(data = as.data.frame(CURRENT_DATA),
                                 meta_data = full_long_3T,
                                 assay_name = "counts", 
                                 tax_level = "Genus",
                                 fix_formula = "hivst", rand_formula = NULL, 
                                 p_adj_method = "BH", pseudo_sens = TRUE, 
                                 prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                                 group = "hivst", struc_zero = TRUE, neg_lb = TRUE,
                                 alpha = 0.05, n_cl = 2, verbose = TRUE,
                                 
                                 iter_control = list(tol = 1e-2, max_iter = 20, 
                                                     verbose = TRUE))
    } else { 
      ancom_model_fit = ancombc2(data = as.data.frame(CURRENT_DATA),
                                 meta_data = full_long_3T,
                                 assay_name = "counts", 
                                 tax_level = "Genus",
                                 fix_formula = "hivst", rand_formula = NULL, 
                                 p_adj_method = "BH", 
                                 
                                 pseudo = c(0.1, 0.5, 1)[ comp_id - 2 ],
                                 pseudo_sens = FALSE, 
                                 
                                 prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                                 group = "hivst", struc_zero = TRUE, neg_lb = TRUE,
                                 alpha = 0.05, n_cl = 2, verbose = TRUE,
                                 
                                 iter_control = list(tol = 1e-2, max_iter = 20, 
                                                     verbose = TRUE))
      ancom_model_fit$res$passed_ss_hivstPositive = NA
    }
    closeAllConnections()
    
    pvalue_mat = ancom_model_fit$res[, c('taxon', 'lfc_hivstPositive', 'se_hivstPositive', 'p_hivstPositive', 'q_hivstPositive', 'passed_ss_hivstPositive')]
    pvalue_mat$significance = ''
    pvalue_mat$significance[ pvalue_mat$q_hivstPositive <= 0.05 ] = 'Significant'
    pvalue_mat = pvalue_mat[order(sign(pvalue_mat$lfc_hivstPositive),
                                  pvalue_mat$p_hivstPositive), ]
    rownames(pvalue_mat) = pvalue_mat$taxon
    
    
    pvalue_mat$"# of samples used" = colSums(t(CURRENT_DATA)[, pvalue_mat$taxon] > 0)    
    
  } 
  
  
  
  
  
  pvalue_mat_ls[[ comp_id ]] = pvalue_mat
  ancom_model_fit_ls[[ comp_id ]] = ancom_model_fit
} 
names(pvalue_mat_ls) = names(ancom_model_fit_ls) = 
  c('Metabolomics', 'Microbiome',
    paste0('Microbiome_pseudo_', c(0.1, 0.5, 1)))


pvalue_mat_ls$Metabolomics = cbind(chem_annotation[rownames(pvalue_mat_ls$Metabolomics), 
                                                   c('CHEM_ID', 'CHEMICAL_NAME', 'SUPER_PATHWAY', 'SUB_PATHWAY')],
                                   pvalue_mat_ls$Metabolomics)
pvalue_mat_ls$Microbiome = cbind(rownames(pvalue_mat_ls$Microbiome), 
                                 pvalue_mat_ls$Microbiome)
colnames(pvalue_mat_ls$Microbiome)[1] = 'Genus'


temp_mat = ancom_model_fit_ls$"Microbiome_pseudo_0.1"$zero_ind 
tab = table(temp_mat$`structural_zero (hivst = Negative)`,
            temp_mat$`structural_zero (hivst = Positive)`)
names(attributes(tab)$dimnames) = c('HIV-', 'HIV+')
tab

taxa_structural_zero = temp_mat[which(temp_mat$`structural_zero (hivst = Negative)` + temp_mat$`structural_zero (hivst = Positive)` == 1), 1:3] 
taxa_structural_zero = taxa_structural_zero[order(taxa_structural_zero[,2],
                                                  taxa_structural_zero[,3],
                                                  taxa_structural_zero[,1]), ]


sig_METAB_vec = rownames(pvalue_mat_ls$Metabolomics)[ which(pvalue_mat_ls$Metabolomics$`BH adjusted pvalue` <= 0.05) ] 
sig_MICRO_vec = rownames(pvalue_mat_ls$"Microbiome_pseudo_0.1")[ which(pvalue_mat_ls$"Microbiome_pseudo_0.1"$p_hivstPositive <= 0.1) ] 


temp_cor_mat = cor_mat[c(sig_MICRO_vec, taxa_structural_zero$taxon), 
                       c(sig_METAB_vec)]

temp_pos_id = as.data.frame(which(temp_cor_mat >= 0.3, arr.ind = TRUE))
temp_pos_id$correlation = temp_cor_mat[ as.matrix(temp_pos_id) ]
temp_pos_id = temp_pos_id[order(- temp_pos_id$correlation), ]
temp_pos_id[, 'row'] = noquote(rownames(temp_cor_mat)[temp_pos_id[, 'row']])
temp_pos_id[, 'col'] = noquote(chem_annotation[colnames(temp_cor_mat)[temp_pos_id[, 'col']], 'CHEMICAL_NAME'])
rownames(temp_pos_id) = NULL
colnames(temp_pos_id) = c('genus', 'chemical', 'correlation')
temp_pos_id$genus_pvalue = pvalue_mat_ls$"Microbiome_pseudo_0.1"[match(temp_pos_id$genus, 
                                                                       pvalue_mat_ls$"Microbiome_pseudo_0.1"$taxon), 'p_hivstPositive']
temp_pos_id$genus_BH_adj_pvalue = pvalue_mat_ls$"Microbiome_pseudo_0.1"[match(temp_pos_id$genus, 
                                                                              pvalue_mat_ls$"Microbiome_pseudo_0.1"$taxon), 'q_hivstPositive']
temp_pos_id$chemical_pvalue = pvalue_mat_ls$Metabolomics[match(temp_pos_id$chemical, 
                                                               pvalue_mat_ls$Metabolomics$CHEMICAL_NAME), 'raw pvalue']
temp_pos_id$chemical_BH_adj_pvalue = pvalue_mat_ls$Metabolomics[match(temp_pos_id$chemical, 
                                                                      pvalue_mat_ls$Metabolomics$CHEMICAL_NAME), 'BH adjusted pvalue']


temp_neg_id = as.data.frame(which(temp_cor_mat <= - 0.3, arr.ind = TRUE))
temp_neg_id$correlation = temp_cor_mat[ as.matrix(temp_neg_id) ]
temp_neg_id = temp_neg_id[order(temp_neg_id$correlation), ]
temp_neg_id[, 'row'] = noquote(rownames(temp_cor_mat)[temp_neg_id[, 'row']])
temp_neg_id[, 'col'] = noquote(chem_annotation[colnames(temp_cor_mat)[temp_neg_id[, 'col']], 'CHEMICAL_NAME'])
rownames(temp_neg_id) = NULL
colnames(temp_neg_id) = c('genus', 'chemical', 'correlation')
temp_neg_id$genus_pvalue = pvalue_mat_ls$"Microbiome_pseudo_0.1"[match(temp_neg_id$genus, 
                                                                       pvalue_mat_ls$"Microbiome_pseudo_0.1"$taxon), 'p_hivstPositive']
temp_neg_id$genus_BH_adj_pvalue = pvalue_mat_ls$"Microbiome_pseudo_0.1"[match(temp_neg_id$genus, 
                                                                              pvalue_mat_ls$"Microbiome_pseudo_0.1"$taxon), 'q_hivstPositive']
temp_neg_id$chemical_pvalue = pvalue_mat_ls$Metabolomics[match(temp_neg_id$chemical, 
                                                               pvalue_mat_ls$Metabolomics$CHEMICAL_NAME), 'raw pvalue']
temp_neg_id$chemical_BH_adj_pvalue = pvalue_mat_ls$Metabolomics[match(temp_neg_id$chemical, 
                                                                      pvalue_mat_ls$Metabolomics$CHEMICAL_NAME), 'BH adjusted pvalue']


pvalue_mat_ls$detailed_correlation = rbind(temp_pos_id,
                                           temp_neg_id)
detaled_cor_mat_part_1 = pvalue_mat_ls$detailed_correlation[which(!is.na(pvalue_mat_ls$detailed_correlation$genus_pvalue)), ] 
detaled_cor_mat_part_2 = pvalue_mat_ls$detailed_correlation[which(is.na(pvalue_mat_ls$detailed_correlation$genus_pvalue)), ] 
### Supplementary Table 9
detaled_cor_mat_part_1
detaled_cor_mat_part_2

### Figure 3b
library(pheatmap); library(RColorBrewer)

par(mfrow = c(1,1))
my_colors = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)
temp_mat_plot = temp_cor_mat[sig_MICRO_vec, ] 
colnames(temp_mat_plot) = noquote(chem_annotation[match(colnames(temp_mat_plot), rownames(chem_annotation)), 
                                                  'CHEMICAL_NAME']) 

row_lab_vec = rep("", nrow(temp_mat_plot))
rownames(temp_mat_plot)
row_lab_vec[ which(rownames(temp_mat_plot) == 'Fusobacterium') ] = "Fusobacterium"
row_lab_vec[ which(rownames(temp_mat_plot) == 'Ruminococcaceae_UCG-003') ] = "Ruminococcaceae_UCG-003"

col_lab_vec = rep("", ncol(temp_mat_plot))
colnames(temp_mat_plot)
col_lab_vec[ which(colnames(temp_mat_plot) == '2-hydroxyglutarate') ] = "2-HG"
col_lab_vec[ which(colnames(temp_mat_plot) == 'phenol glucuronide') ] = "PhG"
col_lab_vec[ which(colnames(temp_mat_plot) == 'hydroxypalmitoyl sphingomyelin (d18:1/16:0(OH))**') ] = "HPSM(d18:1/16:0(OH))"
col_lab_vec[ which(colnames(temp_mat_plot) == 'glycosyl ceramide (d18:1/20:0, d16:1/22:0)*') ] = "glycosyl Cer (d18:1/20:0, d16:1/22:0)"
col_lab_vec[ which(colnames(temp_mat_plot) == 'palmitoyl sphingomyelin (d18:1/16:0)') ] = "PSM (d18:1/16:0)"
col_lab_vec[ which(colnames(temp_mat_plot) == 'sphingomyelin (d18:2/23:1)*') ] = "SM (d18:2/23:1)"
col_lab_vec[ which(colnames(temp_mat_plot) == '5alpha-pregnan-3beta,20alpha-diol monosulfate (2)') ] = "5α-pregnan-3β,20α-ylene sulfate"
col_lab_vec[ which(colnames(temp_mat_plot) == 'pregnenolone sulfate') ] = "PREGS"
pheatmap(temp_mat_plot,
         cluster_cols = FALSE, cluster_rows = FALSE,
         show_rownames = TRUE, show_colnames = TRUE,
         fontsize_number = 8,
         labels_row = row_lab_vec,
         labels_col = col_lab_vec,
         main = 'Correlations between 14 genera and 87 metabolites',
         fontsize = 8,
         color = my_colors,
         breaks = seq(-0.51, 0.51, length.out = 101)
)


### Supplementary Figure 2 
library(pheatmap); library(RColorBrewer)
par(mfrow = c(1,1))
my_colors = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)
temp_mat_plot = temp_cor_mat[taxa_structural_zero$taxon, ] 
colnames(temp_mat_plot) = noquote(chem_annotation[match(colnames(temp_mat_plot), rownames(chem_annotation)), 
                                                  'CHEMICAL_NAME']) 


row_lab_vec = rep("", nrow(temp_mat_plot))
rownames(temp_mat_plot)
row_lab_vec[ which(rownames(temp_mat_plot) == 'Tyzzerella_4') ] = "Tyzzerella_4"
row_lab_vec[ which(rownames(temp_mat_plot) == 'Peptoniphilus') ] = "Peptoniphilus"
row_lab_vec[ which(rownames(temp_mat_plot) == 'Rikenellaceae_RC9_gut_group') ] = "Rikenellaceae_RC9_gut_group"

col_lab_vec = rep("", ncol(temp_mat_plot))
sort(colnames(temp_mat_plot))
col_lab_vec[ which(colnames(temp_mat_plot) == '2-hydroxyglutarate') ] = "2-HG"
col_lab_vec[ which(colnames(temp_mat_plot) == 'palmitoyl sphingomyelin (d18:1/16:0)') ] = "PSM (d18:1/16:0)"
col_lab_vec[ which(colnames(temp_mat_plot) == 'glycosyl ceramide (d18:1/20:0, d16:1/22:0)*') ] = "glycosyl Cer (d18:1/20:0, d16:1/22:0)"
col_lab_vec[ which(colnames(temp_mat_plot) == 'hydroxypalmitoyl sphingomyelin (d18:1/16:0(OH))**') ] = "HPSM(d18:1/16:0(OH))"
pheatmap(temp_mat_plot,
         cluster_cols = FALSE, cluster_rows = FALSE,
         show_rownames = TRUE, show_colnames = TRUE,
         fontsize_number = 8,
         labels_row = row_lab_vec,
         labels_col = col_lab_vec,
         main = 'Correlations between 17 structural zero genera and 87 metabolites',
         fontsize = 8,
         color = my_colors,
         breaks = seq(-0.51, 0.51, length.out = 101)
)



