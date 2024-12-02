library(tidyverse)
source('scripts/custom_functions.R')

fticr_transf <- load_fticr_transformations() %>% 
  select(-sampling, -Transformation) %>% 
  select(SampleID = sample, everything()) 

write_csv(fticr_transf, 'output_tables/Supplementary_table_4.csv')  
