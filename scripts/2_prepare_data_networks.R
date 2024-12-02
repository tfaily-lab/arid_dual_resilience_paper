library(tidyverse)
source('scripts/custom_functions.R')

# Load data

## CoverM

coverm_data <- load_mag_abundance(prev_min_counts = 0,
                                  prev_min_samples = 2)

mag_taxonomy <- load_mag_taxonomy() %>% 
  rename(true_id = bin)

bin_code <- tibble(Genome = colnames(coverm_data)) %>% 
  mutate(ID = paste0('bin', n():1))

write_csv(bin_code, 'data/networks/mag_network_names_updated.csv')

coverm_dry <- coverm_data %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'Genome') %>% 
  inner_join(bin_code) %>% 
  select(ID, contains("SAMP1"), contains('SAMP8'), -contains('22')) %>% 
  rowwise() %>% 
  mutate(sum = sum(c_across(-ID) > 0),
         across(everything(), ~ifelse(.x == 0, NA, .x))) %>% 
  ungroup() %>% 
  filter(sum > 3) %>% 
  select(-sum) 

write_tsv(coverm_dry, 'data/networks/mags_dry_updated.tsv',
          na = '')

coverm_wet <- coverm_data %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'Genome') %>% 
  inner_join(bin_code) %>% 
  select(ID, contains("SAMP5"), contains('SAMP6')) %>% 
  rowwise() %>% 
  mutate(sum = sum(c_across(-ID) > 0),
         across(everything(), ~ifelse(.x == 0, NA, .x))) %>% 
  ungroup() %>% 
  filter(sum > 3) %>% 
  select(-sum)

write_tsv(coverm_wet, 'data/networks/mags_wet_updated.tsv',
          na = '')
