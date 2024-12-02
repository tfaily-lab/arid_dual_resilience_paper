library(tidyverse)
source('scripts/custom_functions.R')

env_data <- load_environmental_data() %>% 
  rownames_to_column(var = 'SampleID')

metadata <- load_metadata() %>% 
  select(SampleID = asv_sample, month, stage)

supp_table3 <- env_data %>% 
  left_join(metadata, by = 'SampleID')

write_csv(supp_table3, 'output_tables/Supplementary_table_3.csv')  
