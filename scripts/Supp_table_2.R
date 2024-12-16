library(tidyverse)


mag_taxonomy <- load_mag_taxonomy()

checkm2 <- read_tsv('data/checkm2/filter_checkm2_quality_report.tsv') %>% 
  mutate(bin_type =case_when(
    Completeness >= 70 & Contamination <= 10 ~ 'High quality',
    Completeness >= 50 & Contamination <= 10 ~ 'Medium quality',
    TRUE ~ 'Low quality'
  )) %>% 
  select(Name, Completeness, Contamination, Contig_N50, Coding_Density, 
         Average_Gene_Length, Genome_Size)

growth_pred <- read_csv('data/grodon/grodon_results.csv') %>% 
  select(bin, CUBHE, CUB, GC, Doubling_Time = d)

final_table <- checkm2 %>% 
  left_join(growth_pred, by = c('Name' = 'bin')) %>% 
  left_join(mag_taxonomy, by = c('Name' = 'bin')) %>% 
  rename('Bin_name' = 'Name')

write_csv(final_table, 'output_tables/Supplementary_table_2.csv')
