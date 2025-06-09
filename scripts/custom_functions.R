library(vegan)
library(tidyverse)

# Function to load MAG taxonomy
load_mag_taxonomy <- function(){
  gtdb_bac <- read_tsv('data/gtdb_res/gtdbtk.bac120.decorated.tree-taxonomy',
                       col_names = FALSE) %>% 
    rename(bin = X1) %>% 
    filter(str_detect(bin, 'bin')) %>% 
    separate_wider_delim(X2,
                         delim = '; ',
                         names = c('Domain', 'Phylum', 'Class', 'Order',
                                   'Family', 'Genus', 'Species'))
  
  
  gtdb_ar <- read_tsv('data/gtdb_res/gtdbtk.ar53.decorated.tree-taxonomy',
                      col_names = FALSE)  %>% 
    rename(bin = X1) %>% 
    filter(str_detect(bin, 'bin')) %>% 
    separate_wider_delim(X2,
                         delim = '; ',
                         names = c('Domain', 'Phylum', 'Class', 'Order',
                                   'Family', 'Genus', 'Species'))
  
  gtdb_all <- rbind(gtdb_bac, gtdb_ar) %>% 
    mutate(across(Domain:Species, ~str_remove(.x, '.__')))
  
  return(gtdb_all)
}

# Function to load MAG abundaces (CoverM)
load_mag_abundance <- function(prev_min_counts = 0,
                               prev_min_samples = 2){
  
  metadata <- load_metadata()
  
  mags_coverm <- read_tsv('data/coverm/coverage_coverm_aln.tsv')  %>% 
    pivot_longer(!Genome, names_to = 'sample', values_to = 'tmm') %>% 
    mutate(sample = str_remove(sample, '_S.*')) %>% 
    group_by(sample, Genome) %>% 
    summarise(tmm = mean(tmm)) %>% 
    ungroup() %>% 
    inner_join(metadata, by = 'sample') %>% 
    select(asv_sample, Genome, tmm) %>% 
    pivot_wider(names_from = Genome, values_from = tmm) %>% 
    column_to_rownames(var = 'asv_sample')
  
  mags_coverm <- mags_coverm[, colSums(mags_coverm > prev_min_counts) >= prev_min_samples]
  
  return(mags_coverm)
  
}

# Function to load ASV data
# For prevalence filtering: 
# prev_min_counts   min number of counts
# prev_min_samples  min number of samples an ASV must appear to be considered present
# normalize_sum     whether to perform sum normalization
load_asv_data <- function(prev_min_counts = 5,
                          prev_min_samples = 2,
                          normalize_sum = TRUE){
  metadata <- read_csv('data/metadata.csv')
  
  asv_data <- read_csv('data/asv_abundances_16S.csv') %>% 
    select(-sequence) %>% 
    rename_with(~str_remove(.x, '-16S.*')) %>%
    pivot_longer(!ASV_ID, names_to = 'asv_sample', values_to = 'counts') %>% 
    inner_join(metadata, by = 'asv_sample') %>% 
    select(asv_sample, ASV_ID, counts) %>% 
    pivot_wider(names_from = ASV_ID, values_from = counts) %>% 
    column_to_rownames(var = 'asv_sample')
  
  asv_data <- asv_data[, colSums(asv_data > prev_min_counts) >= prev_min_samples]
  
  if(normalize_sum){
    asv_data <- decostand(asv_data, method = 'total')
  }
  
  return(asv_data)
}

# Function to load ASV rarefied data
load_asv_data_rarefied <- function(prev_min_counts = 5,
                                   prev_min_samples = 2,
                                   normalize_sum = TRUE){
  metadata <- read_csv('data/metadata.csv')
  
  asv_data <- read_csv('data/preprocessed_data/ASV_matrix_rarefied.csv') %>% 
    select(sampleid, ASV, abundance) %>% 
    pivot_wider(names_from = ASV, values_from = abundance) %>% 
    column_to_rownames(var = 'sampleid')
  
  asv_data <- asv_data[, colSums(asv_data > prev_min_counts) >= prev_min_samples]
  
  if(normalize_sum){
    asv_data <- decostand(asv_data, method = 'total')
  }
  
  return(asv_data)
}

# Load ASV taxonomy
load_asv_taxonomy <- function(){
  
  df <- read_csv('data/preprocessed_data/ASV_matrix_rarefied.csv') %>% 
    select(-ASV_seq, -abundance, -sampleid) %>% 
    distinct()
  
  return(df)
}

load_singlem_data <- function(which,
                              prev_min_counts = 0,
                              prev_min_samples = 2,
                              normalize_sum = TRUE){
  
  metadata <- read_csv('data/metadata.csv')
  
  singlem <- read_csv('data/preprocessed_data/singlem_rplP.csv') %>% 
    filter(!str_detect(sample, 'S0'))  %>% 
    mutate(sample = str_remove(sample, '_L.*')) 
  
  if(which == 'counts'){
    
    singlem <- singlem %>% 
      mutate(sample = str_remove(sample, '_S.*')) %>% 
      select(otu_id, sample, value = num_hits, sequence)
    
  } else if (which == 'coverage') {
    singlem <- singlem %>% 
      mutate(sample = str_remove(sample, '_S.*')) %>% 
      select(otu_id, sample, value = coverage, sequence)
  }
  
  singlem <- singlem %>%  
    group_by(sample, otu_id) %>%
    summarise(value = sum(value)) %>%
    ungroup() %>% 
    inner_join(metadata, by = 'sample') %>% 
    select(asv_sample, otu_id, value) %>% 
    pivot_wider(names_from = otu_id, values_from = value, values_fill = 0) %>% 
    column_to_rownames(var = 'asv_sample')
  
  singlem_filt <- singlem[, colSums(singlem >= prev_min_counts) >= prev_min_samples]
  
  if(normalize_sum){
    singlem_filt <- decostand(singlem_filt, method = 'total')
  }
  
  return(singlem_filt)
}

load_singlem_taxonomy <- function(){
  
  metadata <- read_csv('data/metadata.csv')
  
  singlem_tax <- read_tsv('data/singlem_output/reads_profile.tsv') %>% 
    filter(!str_detect(sample, 'S0'))  %>% 
    mutate(sample = str_remove(sample, '_L.*'),
           sample = str_remove(sample, '_S.*')) %>% 
    group_by(sample, taxonomy) %>%
    summarise(coverage = sum(coverage)) %>%
    separate(taxonomy, 
             into = c('Root', 'Domain', 'Phylum', 'Class', 'Order', 
                      'Family', 'Genus', 'Species'),
             sep = '; ') %>% 
    filter(!is.na(Phylum),
           Domain != 'd__Eukaryota') %>% 
    inner_join(metadata, by = 'sample') %>% 
    select(sample = asv_sample, coverage, Domain, Phylum, Class, Order, Family, Genus, Species) %>% 
    mutate(across(Domain:Species, ~str_remove(.x, '.__')))
  
  return(singlem_tax)
}

load_fticr_data <- function(mode,
                            prev_min_counts = 0,
                            prev_min_samples = 2,
                            normalize_sum = TRUE){
  if(mode == 'positive'){
    file <- 'data/metabodirect_output/fticr_pos_report_processed_noNorm.csv'
  } else if(mode == 'negative'){
    file <- 'data/metabodirect_output/fticr_neg_report_processed_noNorm.csv'
  }
  
  metadata_fticr <- read_csv('data/metadata_fticr.csv') %>% 
    select(SampleID, treatment_code)
  
  fticr_data <- read_csv(file) %>% 
    select(Mass, contains('Pos'), contains('Neg')) %>% 
    rename_with(~str_remove(.x, '_(P|N).*')) %>% 
    pivot_longer(!Mass, names_to = 'SampleID', values_to = 'value') %>% 
    inner_join(metadata_fticr) %>% 
    select(-SampleID) %>% 
    pivot_wider(names_from = Mass, values_from = value) %>% 
    column_to_rownames(var = 'treatment_code')
  
  fticr_data <- fticr_data[, colSums(fticr_data > prev_min_counts) >= prev_min_samples]
  
  if(normalize_sum){
    fticr_data <- decostand(fticr_data, method = 'total')
  }
  
  return(fticr_data)
  
}

load_fticr_annot <- function(mode){
  
  if(mode == 'positive'){
    file <- 'data/metabodirect_output/Report_processed_pos.csv'
  } else if(mode == 'negative'){
    file <- 'data/metabodirect_output/Report_processed_neg.csv'
  }
  
  fticr_data <- read_csv(file) %>% 
    select(Mass, MolecularFormula, El_comp, Class, OC, HC,
           NOSC, GFE, DBE, DBE_O, AI, AI_mod, DBE_AI) %>% 
    drop_na(MolecularFormula)
  
  return(fticr_data)
}


load_fticr_transformations <- function(mode,
                                       prev_min_counts = 0,
                                       prev_min_samples = 2,
                                       normalize_sum = TRUE){
  
  
  metadata_fticr <- read_csv('data/metadata_fticr.csv') %>% 
    select(SampleID, treatment_code)
  
  metadata <- load_metadata()
  
  transf <- read_csv('data/metabodirect_output/transformations_fticr.csv') %>% 
    mutate(SampleID = str_remove(SampleID, '_(P|N).*')) %>% 
    inner_join(metadata_fticr) %>% 
    select(-SampleID) %>% 
    rename(sample = treatment_code) %>% 
    inner_join(metadata, by = c('sample' = 'asv_sample')) %>% 
    select(-sample.y)
  
  return(transf)
  
}


load_environmental_data <- function(){
  
  ndvi <- read_csv('data/environmental_data/ndvi_data.csv') %>% 
    mutate(date = parse_date(`system:time_start`, '%b %d, %Y'))
  
  noaa_data <- read_csv('data/environmental_data/precipitation_data_daily_tucson_11_updated_2021-2023.csv') 
  
  prec <- noaa_data %>% 
    select(LST_DATE, P_DAILY_CALC) %>% 
    mutate(date = ymd(LST_DATE),
           P_DAILY_CALC = ifelse(P_DAILY_CALC < 0, NA, P_DAILY_CALC)) %>% 
    mutate(sampling = case_when(
      between(date, ymd(20210503), ymd(20210506)) ~ 'samp1',
      between(date, ymd(20210612), ymd(20210615)) ~ 'samp2',
      between(date, ymd(20210702), ymd(20210705)) ~ 'samp3',
      between(date, ymd(20210712), ymd(20210715)) ~ 'samp4',
      between(date, ymd(20210724), ymd(20210727)) ~ 'samp5',
      between(date, ymd(20210812), ymd(20210815)) ~ 'samp6',
      between(date, ymd(20210912), ymd(20210915)) ~ 'samp7',
      between(date, ymd(20211015), ymd(20211018)) ~ 'samp8',
      between(date, ymd(20220525), ymd(20220528)) ~ 'samp9'
    )) %>% 
    filter(!is.na(sampling)) %>% 
    group_by(sampling) %>% 
    summarise(date = max(date),
              precipitation = sum(P_DAILY_CALC, na.rm = TRUE))
  
  temp <-  noaa_data %>% 
    select(LST_DATE, SOIL_TEMP_10_DAILY, T_DAILY_AVG) %>% 
    mutate(date = ymd(LST_DATE),
           air_temperature = ifelse(T_DAILY_AVG < 0, NA, T_DAILY_AVG),
           soil_temperature = ifelse(SOIL_TEMP_10_DAILY < 0, NA, 
                                     SOIL_TEMP_10_DAILY)) %>% 
    mutate(sampling = case_when(
      date == ymd(20210506) ~ 'samp1',
      date == ymd(20210615) ~ 'samp2',
      date == ymd(20210705) ~ 'samp3',
      date == ymd(20210715) ~ 'samp4',
      date == ymd(20210727) ~ 'samp5',
      date == ymd(20210815) ~ 'samp6',
      date == ymd(20210915) ~ 'samp7',
      date == ymd(20211018) ~ 'samp8',
      date == ymd(20220528) ~ 'samp9'
    )) %>% 
    filter(!is.na(sampling)) %>% 
    select(date, air_temperature, soil_temperature)
  
  ndvi_4mod <- ndvi %>% 
    mutate(time = as.numeric(difftime(date, ymd(20210509), units = 'days')))  %>% 
    filter(!is.na(NDVI),
           date <= ymd(20220528))
  
  ndvi_mod <- loess(NDVI ~ time, data = ndvi_4mod)
  
  ndvi_ready <- tibble(date = temp$date) %>% 
    mutate(time = as.numeric(difftime(date, ymd(20210509), units = 'days')),
           NDVI = predict(ndvi_mod, time),
           origin = 'prediction') %>% 
    select(-time) %>% 
    mutate(NDVI = case_when(date == ymd(20210506) ~ 1386,
                            date == ymd(20220528) ~ 1669,
                            TRUE ~ NDVI))
  
  env_data <- read_csv('data/environmental_data/envs_matrix_11_18_24.csv') %>% 
    mutate(sample = str_replace(Sample_ID, '_', '-')) %>% 
    left_join(prec) %>% 
    left_join(temp) %>% 
    left_join(ndvi_ready %>% select(-origin)) %>%
    column_to_rownames(var = 'sample') %>% 
    select(-date, -Sample_ID, -sampling, -stage, -month, -site)
  
  return(env_data)
}

load_metadata <- function(){
  metadata <- read_csv('data/metadata.csv') %>% 
    mutate(month = factor(month, levels = c('May', 'June', 'July_1',
                                            'July_2', 'July_3', 'August',
                                            'September', 'October', 'May_22')),
           stage = factor(stage, levels = c('Pre monsoon', 'During monsoon', 'Post monsoon')))
}

load_colors <- function(variable){
  
  mt <- suppressMessages(load_metadata())
  
  colors <- list('stage' = set_names(c('#FEDD7F', '#9bc1bc', '#ed6a5a'),
                                     nm = levels(mt$stage)),
                 'month' = set_names(c('orange', '#FEDD7F',
                                       get_palette('PuBuGn', 5)[2:5],
                                       '#d55f51', 'sienna4', '#f4f1bb'),
                                     nm = levels(mt$month)))
  
  return(colors[[variable]])
  
}

create_network_zipi <- function(edges, nodes){
  
  graph <- graph_from_data_frame(edges, 
                                 directed = FALSE, 
                                 vertices = nodes)
  
  node_par <- net_par(graph, mode = 'all')$v_index
  
  greedy_modules <- module_detect(graph, method = 'cluster_fast_greedy')
  
  nodes_zipi <- zp_analyse(greedy_modules) %>% 
    MetaNet::get_v() %>% 
    mutate(module = paste0('Module_', module)) %>% 
    left_join(node_par)
  
  hubs <- nodes_zipi %>% 
    filter(roles == 2) %>% 
    pull(name)
  
  updated_edges <- edges %>% 
    mutate(connection_with = ifelse(to %in% hubs | from %in% hubs, 
                                    'Link with hubs', 'Link with other MAG'))
  
  updated_graph <- graph_from_data_frame(updated_edges, 
                                         directed = FALSE, 
                                         vertices = nodes_zipi)
  
  return(updated_graph)
  
}

network_stats <- function(graph){
  
  cluster <- cluster_fast_greedy(graph)
  
  df <- net_properties(graph) %>%
    as.data.frame %>% 
    rownames_to_column(var = 'index') %>% 
    add_row(index = 'transitivity',
            value = transitivity(graph)) %>% 
    add_row(index = 'num.modules',
            value = length(cluster)) %>%
    add_row(index = 'modularity',
            value = modularity(cluster)) 
  
  return(df)
  
}

custom_theme <- theme_bw() +
  theme(plot.title = element_text(size = 7, face = 'bold', hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        strip.text = element_text(size = 5),
        strip.background = element_blank(),
        legend.key.size = unit(0.3, 'lines'))


