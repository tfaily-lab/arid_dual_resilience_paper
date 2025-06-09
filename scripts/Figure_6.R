
## Author: Viviana Freire-Zapata

library(patchwork)
library(ggh4x)
library(rstatix)
library(tidyverse)

# Loading metaT raw data

metadata <- read_csv('data/expression/metadata.csv') 

metadata_ready <- metadata %>% 
  mutate(Month = factor(Month, levels = c("May", "July_1", "July_2", "July_3",
                                          "August", "October", "May_22")))

## Phylum colors

mag_taxonomy <- read_csv('data/expression/mag_taxonomy_denovo.csv') %>% 
  mutate(Phylum_lump = fct_lump_n(Phylum, n = 10))


## Expression calculated reaction level

module_exp <- read_csv("data/expression/expression_per_sample_cnps_10-31_fixed.csv")

def <- read_csv("data/expression/definitions_C_N_S_P_11_14_24.csv")

stress_def <- read_csv("data/expression/traits_definitions_update_10_9.csv")

stress <- read_csv("data/expression/expression_per_sample_traits_10-31.csv")

## Joining expression with annotation of modules


stress_ready <- stress %>%
  left_join(stress_def, by = c('module_id' = 'def_id'))

kegg_ready <- module_exp %>% 
  left_join(def, by = c('module_id' = 'Module_id'))


## Loading normalized expression

## geTMM

geTMM <- read_csv("data/expression/geTMM_normalized_matrix_2024.csv")

geTMM_long <- geTMM %>% 
  pivot_longer(!Geneid, names_to = "sampleid", values_to = "geTMM")

## Annotations matrix

annotations <- read_tsv("data/expression/annotations_merged.tsv") %>% 
  rename(Geneid = `...1`)

## Environmental factors

env <- read_csv("data/expression/env_data_ready_11-19-24.csv") %>% 
  mutate(sample = str_replace(sample, "-", "_"))


metadata_env <- read_csv("data/assembly_microbes/metadata_arid_field.csv")


env_ready <- env %>% 
  inner_join(metadata_env, by = c("sample" = "sampleid")) %>% 
  filter(!month == "june" & !month == "september") %>% 
  mutate(sample = str_replace(sample, "_", "-"),
         C_N = total_C/total_N)
  
###### KEGG per reaction definitions #######

## Step 1 = filter those modules with expression zero in all samples

rx_ft_kegg <- module_exp %>% 
  left_join(metadata_ready, by = 'SampleID') %>% 
  group_by(module_id, bin, Month) %>% 
  summarise(expression = mean(expression)) %>% 
  filter(expression > 0) %>% 
  left_join(mag_taxonomy, by = 'bin') %>% 
  left_join(def, by = c('module_id' = 'Module_id'))

# # Stress 

stress_ready <- stress_ready %>%
  left_join(metadata_ready, by = 'SampleID') %>%
  group_by(module_id, bin, Name, Month, Category_1, Category_2) %>%
  summarise(expression_mean = mean(expression)) %>%
  filter(expression_mean > 0) %>%
  left_join(mag_taxonomy, by = 'bin')


## Nitrogen cycle

keep_path <- c("Nitrification", "Denitrification", "Dissimilatory nitrate reduction",
               "Assimilatory nitrate reduction")

rx_nitrogen <- rx_ft_kegg %>% 
  filter(Category_2 == 'Nitrogen metabolism') %>% 
  filter(Phylum == 'Thermoproteota') %>% 
  filter(Category_1 %in% keep_path)

## Heatmap Nitrogen pathways

color_strips <- strip_nested(
  background_y = elem_list_rect(fill = c("white")),
  background_x = elem_list_rect(fill = c("white")))

heat_thermo_nitro <- rx_nitrogen %>% 
  mutate(log_getmm = log10(expression),
         Category_1 = factor (Category_1, levels = c("Nitrification", "Denitrification", 
                                                     "Dissimilatory nitrate reduction",
                                                     "Assimilatory nitrate reduction"))) %>% 
  ggplot()+
  geom_tile(aes(x = Month,
                y = bin,
                fill = log_getmm),
            color = 'white')+
  scale_fill_distiller(palette = 'PuBuGn', direction = 1)+
  labs(y = 'MAG',
       fill = 'log10(geTMM)')+
  facet_grid2(Genus  ~ Category_1 + Name,
              scales = 'free',
              space = 'free',
              switch = 'y',
              strip = color_strips,
              labeller = label_wrap_gen(15)) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        strip.text.y.left = element_text(angle = 90),
        strip.placement = 'outer',
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        #axis.text.y = element_blank(),
        legend.key.size = unit(0.5, 'lines'),
        legend.key.height = unit(1, 'lines'),
        legend.box.spacing = unit(0.5, 'line'))


heat_thermo_nitro

## Stress tolerance 

selected_paths <- c('Oxidative stress', 'pH stress', 'Transcription regulation')

genes_vip <- c('urease', 'TFB')


heat_thermo_stress <- stress_ready %>% 
  filter(Phylum == 'Thermoproteota') %>% 
  filter(Category_1 %in% selected_paths) %>% 
  filter(Name %in% genes_vip) %>% 
  mutate(log_getmm = log10(expression_mean),
         Name = str_replace(Name, 'urease', 'urease complex')) %>% 
  ggplot()+
  geom_tile(aes(x = Month,
                y = bin,
                fill = log_getmm),
            color = 'white')+
  scale_fill_distiller(palette = 'PuBuGn', direction = 1)+
  labs(y = 'MAG',
       fill = 'log10(geTMM)')+
  facet_grid2(Genus  ~ Name,
              scales = 'free',
              space = 'free',
              switch = 'y',
              strip = color_strips) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        strip.text.y.left = element_text(angle = 90),
        strip.placement = 'outer',
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        legend.key.size = unit(0.5, 'lines'),
        legend.key.height = unit(1, 'lines'),
        legend.box.spacing = unit(0.5, 'line'))


heat_thermo_stress


## Correlations 

## Genes of interest

nitrification <- c('K10944', 'K10945', 'K10946')

## Nitrification

nitri_ft <- annotations %>% 
  filter(ko_id %in% nitrification) %>% 
  select(Geneid, fasta, ko_id, kegg_hit, contains("peptidase")) %>% 
  inner_join(mag_taxonomy, by = c('fasta' = 'bin')) %>% 
  inner_join(geTMM_long, by = 'Geneid') %>% 
  inner_join(metadata_ready, by = c('sampleid' = 'SampleID')) %>% 
  filter(Phylum == "Thermoproteota") %>%
  group_by(ko_id, fasta, sampleid, Month, Phylum, Class, Order, Family, Genus) %>% 
  summarise(geTMM_sum = sum(geTMM))


## Preparing table for correlation 


nitrification_ready <- nitri_ft %>% 
  filter(!Month == 'May_22') %>% 
  mutate(sampleid_match = str_replace(sampleid, "RNA-", "SAMP")) %>% 
  inner_join(env_ready, by = c('sampleid_match' = 'sample')) %>% 
  split(.$fasta)


nitrification_complex <- imap(nitrification_ready, function(mag, name){
  
  mean <- mag %>% 
    group_by(sampleid, fasta, Month, sampleid_match, total_C , total_N,
             ph_mean, moisture, C_N, NDVI) %>% 
    summarise(geTMM_mean = mean(geTMM_sum))
  
  return(mean)
  
})

nitri_complex_ready <- reduce(nitrification_complex, rbind) %>% 
  filter(!fasta == 'bco_bin_78_1')## delete mag that only have one subunit

## Correlation with ENVS factors

nitri_complx_corre <- nitri_complex_ready %>% 
  group_by(fasta) %>% 
  cor_test(vars = "geTMM_mean",
           vars2 = c("total_C", "total_N", "ph_mean", 
                     "moisture", "C_N", "NDVI"),
           method = "pearson", 
           alternative = "two.sided") %>% 
  adjust_pvalue(method = 'fdr')  


## Plotting only significants


complx_sig <- nitri_complx_corre %>% 
  filter(p.adj < 0.05)

label_position <- tribble(
  ~fasta, ~env, ~x, ~y,
  'bco_bin_91', 'Moisture [%]', 5, 4100,
  'bco_bin_91', 'NDVI', 2250, 4100,
  'DNA-1_bin_53_1', 'Moisture [%]', 7.5, 1750,
  'DNA-1_bin_53_1', 'NDVI', 2500, 1750,
  'bcn_bin_142_1', 'Moisture [%]', 10, 0,
  'bcn_bin_142_1', 'NDVI', 1500, 0,
)

label_data <- complx_sig %>% 
  select(fasta, env = var2, cor, p.adj) %>% 
  mutate(env = ifelse(env == 'moisture', 'Moisture [%]', env)) %>% 
  mutate(cor = paste0('Pearson r = ', round(cor, 2)),
         p.adj = paste0('P adjusted = ', round(p.adj, 4)),
         lab = glue::glue("{cor}\n{p.adj}")) %>%
  left_join(label_position)


nitri_plot <- nitri_complex_ready %>% 
  filter(fasta %in% complx_sig$fasta)



plot_1cx <- nitri_plot %>% 
  ungroup() %>% 
  filter(fasta %in% complx_sig$fasta) %>% 
  mutate(Month = factor(Month, levels = c("May", "July_1", "July_2",
                                          "July_3", "August", "October"))) %>% 
  select(fasta, geTMM_mean, moisture, NDVI) %>% 
  pivot_longer(c(moisture, NDVI), names_to = 'env', values_to = 'value') %>% 
  mutate(env = ifelse(env == 'moisture', 'Moisture [%]', env)) %>% 
  ggplot() +
  geom_point(aes(x = value,
                 y = geTMM_mean,
                 color = fasta))+
  geom_smooth(aes(x = value,
                  y = geTMM_mean,
                  color = fasta),
              method = lm) +
  facet_wrap(~ env,
             scales = 'free_x',
             switch = 'x',
             ncol = 2) +
  scale_color_brewer(palette = 'Dark2') +
  geom_label(data = label_data,
            aes(x = x,
                y = y,
                label = lab,
                color = fasta),
            show.legend = FALSE,
            hjust = 0,
            vjust = 0,
            size = 2.2) +
  guides(color = guide_legend(nrow = 1)) +
  # annotate('text', x = 18, y = 2000, label = "Pearson r = -0.66")+
  # annotate('text', x = 18, y = 1800, label = "P.adj = 0.00395")+
  labs(y = "geTMM - amoABC",
       color = "MAG")+
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.title.x = element_blank(),
        strip.placement = 'outer',
        legend.position = 'bottom',
        legend.key.size = unit(0.5, 'lines'),
        panel.grid = element_blank())


plot_1cx



layout_plot <- "
A
B
C
"

final <- heat_thermo_nitro + plot_1cx + heat_thermo_stress +
  plot_layout(nrow = 2, design = layout_plot) +
  plot_annotation(tag_levels = 'A') &
  theme(text = element_text(size = 7))

final


ggsave('output_figures/Figure_6.svg', final,
       dpi = 300, height = 210, width = 180, units = "mm")

ggsave('output_figures/Figure_6.png', final,
       dpi = 300, height = 210, width = 180, units = "mm")


## Thermoproteota MAGs expression table - KOs level


keeg_hit_ready <- annotations %>% 
  dplyr::select(ko_id, kegg_hit) %>% 
  filter(!is.na(ko_id)) %>% 
  distinct()

kos_thermo <- annotations %>% 
  filter(!is.na(ko_id)) %>% 
  select(Geneid, fasta, ko_id, kegg_hit) %>% 
  inner_join(mag_taxonomy, by = c('fasta' = 'bin')) %>% 
  inner_join(geTMM_long, by = 'Geneid') %>% 
  inner_join(metadata_ready, by = c('sampleid' = 'SampleID')) %>% 
  filter(Phylum == "Thermoproteota") %>%
  group_by(ko_id, fasta, sampleid, Month, Phylum, Class, Order, Family, Genus) %>% 
  summarise(geTMM_sum = sum(geTMM)) %>% 
  group_by(ko_id, fasta, Month, Phylum, Class, Order, Family, Genus) %>% 
  summarise(geTMM_mean = mean(geTMM_sum)) %>% 
  filter(geTMM_mean > 0) %>% 
  left_join(keeg_hit_ready, by = 'ko_id')
  
#write_csv(kos_thermo, 'output_tables/Supplementary_table_5.csv')

