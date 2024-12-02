library(tidyverse)
library(vegan)

# Creating IDs for singleM data

singlem <- read_tsv('data/singlem_output/reads_otu_table.tsv')

singlem_rplP <- singlem %>% 
  filter(gene == 'S3.4.ribosomal_protein_L16_L10E_rplP',
         !str_detect(taxonomy, 'Eukaryota'))

singlem_id <- singlem_rplP %>% 
  select(sequence) %>% 
  distinct() %>% 
  mutate(otu_id = paste0('singlem_rplP_OTU_', str_pad(n():1, 6, pad = '0')))

singlem_rplP_id <- singlem_rplP %>% 
  inner_join(singlem_id) %>% 
  select(otu_id, everything())

write_csv(singlem_rplP_id, 'data/preprocessed_data/singlem_rplP.csv')

# Creating ASV rarefied table

asv <- read_csv('data/dada2_output/AVS_matrix_final_2024.csv')

asv_wide <- asv %>% 
  select(ASV, sampleid, abundance) %>% 
  pivot_wider(names_from = ASV, values_from = abundance) %>% 
  column_to_rownames(var = 'sampleid')

asv_annot <- asv %>% 
  select(-sampleid, -abundance) %>% 
  distinct()

counts_sample <- rowSums(asv_wide)

asv_rare <- rrarefy(asv_wide, 25000) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'sampleid') %>% 
  pivot_longer(!sampleid, names_to = 'ASV', values_to = 'abundance') %>% 
  left_join(asv_annot)

#write_csv(asv_rare, 'data/preprocessed_data/ASV_matrix_rarefied.csv')