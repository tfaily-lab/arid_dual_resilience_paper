library(tidyverse)

source('scripts/custom_functions.R')

mag_taxonomy <- load_mag_taxonomy()

checkm2 <- read_tsv('data/checkm2/filter_checkm2_quality_report.tsv') %>% 
  select(Name, Completeness, Contamination, Contig_N50, Coding_Density, 
         Average_Gene_Length, Genome_Size)

growth_pred <- read_csv('data/grodon/grodon_results.csv') %>% 
  select(bin, CUBHE, CUB, GC, Doubling_Time = d)


trnas <- read_tsv('data/dram_annotations/trnas_merged.tsv') %>% 
  filter(is.na(note)) %>% 
  group_by(scaffold) %>% 
  count() %>% 
  rename(bin = scaffold, `# tRNAs` = n)


rrnas <- read_tsv('data/dram_annotations/rrnas_merged.tsv') %>% 
  select(fasta, scaffold, type) %>% 
  pivot_wider(names_from = 'type', values_from = 'scaffold', values_fn = ~paste0(.x, collapse = ','))

final_table <- checkm2 %>% 
  left_join(growth_pred, by = c('Name' = 'bin')) %>% 
  left_join(mag_taxonomy, by = c('Name' = 'bin')) %>% 
  left_join(trnas, by = c('Name' = 'bin')) %>% 
  left_join(rrnas, by = c('Name' = 'fasta')) %>% 
  rename('Bin_name' = 'Name') %>% 
  mutate(bin_type =case_when(
    (Completeness >= 90 & Contamination <= 5 & `# tRNAs` >= 18 &
      !is.na(`5S rRNA`) & !is.na(`16S rRNA`) & !is.na(`23S rRNA`))~ 'High quality',
    Completeness >= 50 & Contamination <= 10 ~ 'Medium quality',
    TRUE ~ 'Low quality'
  )) 

write_csv(final_table, 'updated_tables/Supplementary_table_6.csv')
