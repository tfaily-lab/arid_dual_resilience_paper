
## Author: Viviana Freire-Zapata

# Load in required libraries
require(reshape2)
require(ggplot2)
require(ggthemes)
require(vegan)
library(dbplyr)
library(ggh4x)
library(ggsci)
library(ggrepel)
library(rstatix)
library(ggpubr)
library(patchwork)
library(tidyverse)


#### Color palette

month_color <- set_names(c('orange', '#FEDD7F',
                           get_palette('PuBuGn', 5)[2:5],
                           '#d55f51', 'sienna4', '#f4f1bb'),
                         nm = c('may', 'june', 'july_1', 'july_2', 'july_3', 
                                'august', 'september', 'october', 'may_22'))

eco_colors <- set_names(c('#e76f51', '#f4a261', '#264653', '#2a9d8f', '#f1faee'),
                        nm = c('Variable selection', 'Homogeneous selection',
                               'Dispersal limitation', 'Homogenizing dispersal', 
                               'Undominated'))
#### Load in data ####

#### Microbiome - ASV 

# Load in bNTI
char_microbe = read_csv("data/assembly_microbes/weighted_BNTI_ASV_rel_abundance.csv") %>% 
  column_to_rownames(var = '...1')


# Load in RCBC
char.rcbc_microbe = read_csv("data/assembly_microbes/RCBC_ASV_arid_rel_abundance.csv") %>% 
  column_to_rownames(var = 'sampleid')

# ######################## #
#### Data Preprocessing ####
# ######################## #


# Matching the order of RCBC to bNTI results
char.rcbc_microbe = char.rcbc_microbe[row.names(char_microbe), colnames(char_microbe), drop = F]


# Setting the RCBC diagonals to NA
diag(char.rcbc_microbe) = NA


# Reflecting null matrices
char_microbe[upper.tri(char_microbe)] = t(char_microbe)[upper.tri(char_microbe)]
char.rcbc_microbe[upper.tri(char.rcbc_microbe)] = t(char.rcbc_microbe)[upper.tri(char.rcbc_microbe)]



# Removing significant bNTI results from the RCBC results
char.rcbc_microbe[abs(char_microbe) > 2] = NA


# ##################################### #
#### Preprocessing Data for plotting ####
# ##################################### #

# Melting data
char_microbe = melt(as.matrix(char_microbe))

char.rcbc_microbe = melt(as.matrix(char.rcbc_microbe))


# Removing null values
char_microbe = char_microbe[!is.na(char_microbe$value),]

char.rcbc_microbe = char.rcbc_microbe[!is.na(char.rcbc_microbe$value),]

### pLOTTING
metadata <- read_csv("data/assembly_microbes/metadata_arid_field.csv") %>% 
  filter(!month == "may_22") 


metadata$month <- factor(metadata$month, levels = c("may","june","july_1", "july_2",
                                                    "july_3", "august", "september", "october"))

###### Ecological processes per whole system ###############

# Determining ecological processes fractionation 
## This is for the whole system

## Renaming columns

bnti_data <- char_microbe %>% 
  rename(BNTI = value)

rcbc_data <- char.rcbc_microbe %>% 
  rename(RCBC = value)

## Joining BNTI and RCBC datasets

data_full_microbe <- full_join(bnti_data, rcbc_data)

## removing duplicated comparisons

temp_data <- data_full_microbe

for(i in 1:nrow(temp_data)){
  temp_data$comb[i] <- paste0(sort(c(as.character(temp_data$Var1[i]), 
                                     as.character(temp_data$Var2[i]))), 
                              collapse = ';')
}

data_full_filt_micro <- temp_data %>% 
  group_by(comb) %>% 
  summarise(BNTI = mean(BNTI),
            RCBC = mean(RCBC, na.rm = TRUE)) %>% 
  separate(comb, into = c('Var1', 'Var2'), sep = ';')


## Determining ecological processes fractionation

eco_microbe <- data_full_filt_micro %>% 
  mutate(ecological = case_when(BNTI <= -2 ~ "Homogeneous selection",
                                BNTI >= 2 ~ "Variable selection",
                                abs(BNTI) <2 & RCBC < -0.95 ~ "Homogenizing dispersal",
                                abs(BNTI) <2 & RCBC > 0.95 ~ "Dispersal limitation",
                                abs(BNTI) <2 & abs(RCBC) <= 0.95 ~ "Undominated")) %>% 
  dplyr::count(ecological) %>% 
  mutate(percentage = n/sum(n)*100)

################################## eco processes per sample


eco_microbe_full <- data_full_microbe %>% 
  mutate(ecological = case_when(BNTI <= -2 ~ "Homogeneous selection",
                                BNTI >= 2 ~ "Variable selection",
                                abs(BNTI) <2 & RCBC < -0.95 ~ "Homogenizing dispersal",
                                abs(BNTI) <2 & RCBC > 0.95 ~ "Dispersal limitation",
                                abs(BNTI) <2 & abs(RCBC) <= 0.95 ~ "Undominated"))


#### Plot by sampling - month  ####


may <- metadata %>% 
  filter(month == "may") %>% 
  pull(sampleid)

june <- metadata %>% 
  filter(month == "june") %>% 
  pull(sampleid)

july_1 <- metadata %>% 
  filter(month == "july_1") %>% 
  pull(sampleid)

july_2 <- metadata %>% 
  filter(month == "july_2") %>% 
  pull(sampleid)

july_3 <- metadata %>% 
  filter(month == "july_3") %>% 
  pull(sampleid)

aug <- metadata %>% 
  filter(month == "august") %>% 
  pull(sampleid)

sept <- metadata %>% 
  filter(month == "september") %>% 
  pull(sampleid)

oct <- metadata %>% 
  filter(month == "october") %>% 
  pull(sampleid)



# Boxplots by like samples (example: Palsa vs Palsa)

bnti_data_may <- data_full_filt_micro %>% 
  filter(Var1 %in% may & Var2 %in% may) %>% 
  mutate(month = "may")

bnti_data_june <- data_full_filt_micro  %>% 
  filter(Var1 %in% june & Var2 %in% june) %>% 
  mutate(month = "june")

bnti_data_july_1 <- data_full_filt_micro %>% 
  filter(Var1 %in% july_1 & Var2 %in% july_1) %>% 
  mutate(month = "july_1")

bnti_data_july_2 <- data_full_filt_micro %>% 
  filter(Var1 %in% july_2 & Var2 %in% july_2) %>% 
  mutate(month = "july_2")

bnti_data_july_3 <- data_full_filt_micro %>% 
  filter(Var1 %in% july_3 & Var2 %in% july_3) %>% 
  mutate(month = "july_3")

bnti_data_aug <- data_full_filt_micro %>% 
  filter(Var1 %in% aug & Var2 %in% aug) %>% 
  mutate(month = "august")


bnti_data_sep <- data_full_filt_micro %>% 
  filter(Var1 %in% sept & Var2 %in% sept) %>% 
  mutate(month = "september")

bnti_data_oct <- data_full_filt_micro %>% 
  filter(Var1 %in% oct & Var2 %in% oct) %>% 
  mutate(month = "october")



# Joining data

bnti_final_microbe <- rbind(bnti_data_may, bnti_data_june, bnti_data_july_1, bnti_data_july_2, bnti_data_july_3,
                    bnti_data_aug, bnti_data_sep, bnti_data_oct)


bnti_final_microbe$month <- factor(bnti_final_microbe$month, levels = c("may", "june", "july_1", "july_2",
                                                        "july_3", "august", "september","october"))


### All samples


all_plot <- data_full_microbe %>% 
  mutate(ecological = case_when(BNTI <= -2 ~ "Homogeneous selection",
                                BNTI >= 2 ~ "Variable selection",
                                abs(BNTI) <2 & RCBC < -0.95 ~ "Homogenizing dispersal",
                                abs(BNTI) <2 & RCBC > 0.95 ~ "Dispersal limitation",
                                abs(BNTI) <2 & abs(RCBC) <= 0.95 ~ "Undominated")) %>% 
  left_join(metadata, by = c('Var1' = 'sampleid')) %>% 
  mutate(month = factor(month, levels = c('may', 'june', 'july_1',
                                          'july_2', 'july_3', 'august',
                                          'september','october')))

new_plot_all <- all_plot %>% 
  mutate(ecological = factor(ecological, levels = c('Variable selection', 
                                                    'Homogeneous selection',
                                                    'Dispersal limitation', 
                                                    'Homogenizing dispersal',
                                                    'Undominated'))) %>% 
  ggplot() +
  # geom_boxplot(aes(x = month,
  #                 y = BNTI,
  #                 fill = month),
  #             alpha = 0.3,
  #             show.legend = FALSE) +
  # ggbeeswarm::geom_beeswarm(aes(x = month,
  #                               y = BNTI,
  #                               color = ecological)) +
  ggridges::geom_density_ridges(aes(x= BNTI,
                                    y = month,
                                    fill = month),
                                alpha = 0.5) +
  geom_vline(xintercept = c(-2,2), color = "red", lty = 2)+
  # geom_point(aes(x = month,
  #                y = BNTI,
  #                color = ecological),
  #            size = 2) +
  scale_y_discrete(limits = rev) +
  scale_fill_manual(values = month_color)+
  #scale_color_manual(values = assembly_color) +
  labs(title = 'Microbiome',
       y = expression(bold(paste(beta,"-nearest taxon index"))),
       color = "Ecological processes") +
  theme_bw() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'none') 

new_plot_all  



########## Ecological processes by month

eco_month <- bnti_final_microbe %>% 
  mutate(ecological = case_when(BNTI <= -2 ~ "Homogeneous selection",
                                BNTI >= 2 ~ "Variable selection",
                                abs(BNTI) <2 & RCBC < -0.95 ~ "Homogenizing dispersal",
                                abs(BNTI) <2 & RCBC > 0.95 ~ "Dispersal limitation",
                                abs(BNTI) <2 & abs(RCBC) <= 0.95 ~ "Undominated")) %>% 
  group_by(month) %>% 
  dplyr::count(ecological) %>% 
  mutate(percentage = n/sum(n)*100) %>% 
  mutate(Analysis = "Month") %>% 
  rename(set = month)

## Bar plot -ecological processees - Month

ecobar_month <- ggplot(data = eco_month, aes(x = set, y = percentage, fill = ecological))+
  geom_bar(stat = "identity", position = "dodge")+
  labs(title = expression(bold("Ecological processes - Microbiome")),
       y = expression(bold("Assembly patterns (%)")),
       fill = "Assembly process")+
  xlab(NULL)+
  ylim(0,100)+
  scale_fill_manual(values = c('purple','pink','indianred2', 'steelblue'))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        axis.title.y = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white', face = 'bold', size = 8))

ecobar_month


###################################### FT-ICR #################################

## Positive mode

# Load in bNTI
char_icr = read_csv("data/assembly_metabolites/Field_bulk_TWCD_bNTI_999_pos_7_11_22.csv") %>% 
  column_to_rownames(var = '...1')


# Load in RCBC
char.rcbc_icr = read_csv("data/assembly_metabolites/Field_bulk_pos_TWCD_RCBC_9999_pos_7-11-22.csv") %>% 
  column_to_rownames(var = '...1')

# ######################## #
#### Data Preprocessing ####
# ######################## #


# Matching the order of RCBC to bNTI results
char.rcbc_icr = char.rcbc_icr[row.names(char_icr), colnames(char_icr), drop = F]


# Setting the RCBC diagonals to NA
diag(char.rcbc_icr) = NA


# Reflecting null matrices
char_icr[upper.tri(char_icr)] = t(char_icr)[upper.tri(char_icr)]


# Removing significant bNTI results from the RCBC results
char.rcbc_icr[abs(char_icr) > 2] = NA


# ##################################### #
#### Preprocessing Data for plotting ####
# ##################################### #

# Melting data
char_icr = melt(as.matrix(char_icr))

char.rcbc_icr = melt(as.matrix(char.rcbc_icr))


# Removing null values
char_icr = char_icr[!is.na(char_icr$value),]

char.rcbc_icr = char.rcbc_icr[!is.na(char.rcbc_icr$value),]

### PLOTTING


###### Ecological processes per whole system ###############

# Determining ecological processes fractionation 
## This is for the whole system

## Renaming columns

bnti_data_icr <- char_icr %>% 
  rename(BNTI = value)

rcbc_data_icr <- char.rcbc_icr %>% 
  rename(RCBC = value)

## Joining BNTI and RCBC datasets

data_full_icr <- full_join(bnti_data_icr, rcbc_data_icr)

## removing duplicated comparisons

temp_data <- data_full_icr

for(i in 1:nrow(temp_data)){
  temp_data$comb[i] <- paste0(sort(c(as.character(temp_data$Var1[i]), 
                                     as.character(temp_data$Var2[i]))), 
                              collapse = ';')
}

data_full_filt_icr <- temp_data %>% 
  group_by(comb) %>% 
  summarise(BNTI = mean(BNTI),
            RCBC = mean(RCBC, na.rm = TRUE)) %>% 
  separate(comb, into = c('Var1', 'Var2'), sep = ';') 


## Determining ecological processes fractionation

eco_icr <- data_full_filt_icr %>% 
  mutate(ecological = case_when(BNTI <= -2 ~ "Homogeneous selection",
                                BNTI >= 2 ~ "Variable selection",
                                abs(BNTI) <2 & RCBC < -0.95 ~ "Homogenizing dispersal",
                                abs(BNTI) <2 & RCBC > 0.95 ~ "Dispersal limitation",
                                abs(BNTI) <2 & abs(RCBC) <= 0.95 ~ "Undominated")) %>% 
  dplyr::count(ecological) %>% 
  mutate(percentage = n/sum(n)*100)

#######################
#### Plot by sampling - month  ####


may <- metadata %>% 
  filter(month == "may") %>% 
  pull(sampleid_icr)

june <- metadata %>% 
  filter(month == "june") %>% 
  pull(sampleid_icr)

july_1 <- metadata %>% 
  filter(month == "july_1") %>% 
  pull(sampleid_icr)

july_2 <- metadata %>% 
  filter(month == "july_2") %>% 
  pull(sampleid_icr)

july_3 <- metadata %>% 
  filter(month == "july_3") %>% 
  pull(sampleid_icr)

aug <- metadata %>% 
  filter(month == "august") %>% 
  pull(sampleid_icr)

sept <- metadata %>% 
  filter(month == "september") %>% 
  pull(sampleid_icr)

oct <- metadata %>% 
  filter(month == "october") %>% 
  pull(sampleid_icr)


# Boxplots by like samples (example: Palsa vs Palsa)

bnti_data_may <- data_full_filt_icr %>% 
  filter(Var1 %in% may & Var2 %in% may) %>% 
  mutate(month = "may")

bnti_data_june <- data_full_filt_icr %>% 
  filter(Var1 %in% june & Var2 %in% june) %>% 
  mutate(month = "june")

bnti_data_july_1 <- data_full_filt_icr %>% 
  filter(Var1 %in% july_1 & Var2 %in% july_1) %>% 
  mutate(month = "july_1")

bnti_data_july_2 <- data_full_filt_icr %>% 
  filter(Var1 %in% july_2 & Var2 %in% july_2) %>% 
  mutate(month = "july_2")

bnti_data_july_3 <- data_full_filt_icr %>% 
  filter(Var1 %in% july_3 & Var2 %in% july_3) %>% 
  mutate(month = "july_3")

bnti_data_aug <- data_full_filt_icr %>% 
  filter(Var1 %in% aug & Var2 %in% aug) %>% 
  mutate(month = "august")


bnti_data_sep <- data_full_filt_icr %>% 
  filter(Var1 %in% sept & Var2 %in% sept) %>% 
  mutate(month = "september")

bnti_data_oct <- data_full_filt_icr %>% 
  filter(Var1 %in% oct & Var2 %in% oct) %>% 
  mutate(month = "october")



# Joining data

bnti_final_icr <- rbind(bnti_data_may, bnti_data_june, bnti_data_july_1, bnti_data_july_2, bnti_data_july_3,
                    bnti_data_aug, bnti_data_sep, bnti_data_oct)


bnti_final_icr$month <- factor(bnti_final_icr$month, levels = c("may", "june", "july_1", "july_2",
                                                        "july_3", "august", "september","october"))


#########################################################

### All samples


all_plot_icr <- data_full_icr %>% 
  mutate(ecological = case_when(BNTI <= -2 ~ "Homogeneous selection",
                                BNTI >= 2 ~ "Variable selection",
                                abs(BNTI) <2 & RCBC < -0.95 ~ "Homogenizing dispersal",
                                abs(BNTI) <2 & RCBC > 0.95 ~ "Dispersal limitation",
                                abs(BNTI) <2 & abs(RCBC) <= 0.95 ~ "Undominated")) %>% 
  left_join(metadata, by = c('Var1' = 'sampleid_icr')) %>% 
  mutate(month = factor(month, levels = c('may', 'june', 'july_1',
                                          'july_2', 'july_3', 'august',
                                          'september','october')))


new_plot_all_icr <- all_plot_icr %>% 
  mutate(ecological = factor(ecological, levels = c('Variable selection', 
                                                    'Homogeneous selection',
                                                    'Dispersal limitation', 
                                                    'Homogenizing dispersal',
                                                    'Undominated'))) %>% 
  ggplot() +
  ggridges::geom_density_ridges(aes(x= BNTI,
                                    y = month,
                                    fill = month),
                                alpha = 0.5) +
  geom_vline(xintercept = c(-2,2), color = "red", lty = 2)+
  scale_y_discrete(limits = rev) +
  scale_fill_manual(values = month_color)+
  labs(title = 'Metabolome',
       y = expression(bold(paste(beta,"-nearest taxon index"))),
       color = "Ecological processes") +
  theme_bw() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'none') 

new_plot_all_icr  

########### Loading metabolome - fractions

# Lipids

char_lipid = read_csv("data/assembly_metabolites/Field_lipids_pos_TWCD_bNTI_999_8_11_22.csv") %>% 
  column_to_rownames(var = '...1')


# Load in RCBC
char.rcbc_lipid = read_csv("data/assembly_metabolites//Field_lipids_pos_TWCD_RCBC_9999_8-12-22.csv") %>% 
  column_to_rownames(var = '...1')

# Sugars

char_sugars = read_csv("data/assembly_metabolites//Field_sugars_pos_TWCD_bNTI_999_8_11_22.csv") %>% 
  column_to_rownames(var = '...1')


# Load in RCBC
char.rcbc_sugars = read_csv("data/assembly_metabolites//Field_sugars_pos_TWCD_RCBC_9999_8-12-22.csv") %>% 
  column_to_rownames(var = '...1')

# TLC

char_tlc = read_csv("data/assembly_metabolites/Field_TLC_pos_TWCD_bNTI_999_8_11_22.csv") %>% 
  column_to_rownames(var = '...1')


# Load in RCBC
char.rcbc_tlc = read_csv("data/assembly_metabolites/Field_TLC_pos_TWCD_RCBC_9999_8-12-22.csv") %>% 
  column_to_rownames(var = '...1')


# ######################## #
#### Data Preprocessing ####
# ######################## #


##Lipids

# Matching the order of RCBC to bNTI results
char.rcbc_lipid = char.rcbc_lipid[row.names(char_lipid), colnames(char_lipid), drop = F]


# Setting the RCBC diagonals to NA
diag(char.rcbc_lipid) = NA


# Reflecting null matrices
char_lipid[upper.tri(char_lipid)] = t(char_lipid)[upper.tri(char_lipid)]


# Removing significant bNTI results from the RCBC results
char.rcbc_lipid[abs(char_lipid) > 2] = NA


# ##################################### #
#### Preprocessing Data for plotting ####
# ##################################### #

# Melting data
char_lipid = melt(as.matrix(char_lipid))

char.rcbc_lipid = melt(as.matrix(char.rcbc_lipid))


# Removing null values
char_lipid = char_lipid[!is.na(char_lipid$value),]

char.rcbc_lipid = char.rcbc_lipid[!is.na(char.rcbc_lipid$value),]


###### Ecological processes per whole system ###############

# Determining ecological processes fractionation 
## This is for the whole system

## Renaming columns

bnti_data_lipid <- char_lipid %>% 
  rename(BNTI = value)

rcbc_data_lipid <- char.rcbc_lipid %>% 
  rename(RCBC = value)

## Joining BNTI and RCBC datasets

data_full_lipid <- full_join(bnti_data_lipid, rcbc_data_lipid)


## removing duplicated comparisons

temp_data <- data_full_lipid

for(i in 1:nrow(temp_data)){
  temp_data$comb[i] <- paste0(sort(c(as.character(temp_data$Var1[i]), 
                                     as.character(temp_data$Var2[i]))), 
                              collapse = ';')
}

data_filt_lipid <- temp_data %>% 
  group_by(comb) %>% 
  summarise(BNTI = mean(BNTI),
            RCBC = mean(RCBC, na.rm = TRUE)) %>% 
  separate(comb, into = c('Var1', 'Var2'), sep = ';')


## Determining ecological processes fractionation

eco_lipid <- data_filt_lipid %>% 
  mutate(ecological = case_when(BNTI <= -2 ~ "Homogeneous selection",
                                BNTI >= 2 ~ "Variable selection",
                                abs(BNTI) <2 & RCBC < -0.95 ~ "Homogenizing dispersal",
                                abs(BNTI) <2 & RCBC > 0.95 ~ "Dispersal limitation",
                                abs(BNTI) <2 & abs(RCBC) <= 0.95 ~ "Undominated")) %>% 
  dplyr::count(ecological) %>% 
  mutate(percentage = n/sum(n)*100)


#### Plot by sampling - month  ####

# Boxplots by like samples (example: May vs May)

bnti_data_may <- data_filt_lipid %>% 
  filter(Var1 %in% may & Var2 %in% may) %>% 
  mutate(month = "may")

bnti_data_june <- data_filt_lipid %>% 
  filter(Var1 %in% june & Var2 %in% june) %>% 
  mutate(month = "june")

bnti_data_july_1 <- data_filt_lipid %>% 
  filter(Var1 %in% july_1 & Var2 %in% july_1) %>% 
  mutate(month = "july_1")

bnti_data_july_2 <- data_filt_lipid %>% 
  filter(Var1 %in% july_2 & Var2 %in% july_2) %>% 
  mutate(month = "july_2")

bnti_data_july_3 <- data_filt_lipid %>% 
  filter(Var1 %in% july_3 & Var2 %in% july_3) %>% 
  mutate(month = "july_3")

bnti_data_aug <- data_filt_lipid %>% 
  filter(Var1 %in% aug & Var2 %in% aug) %>% 
  mutate(month = "august")


bnti_data_sep <- data_filt_lipid %>% 
  filter(Var1 %in% sept & Var2 %in% sept) %>% 
  mutate(month = "september")

bnti_data_oct <- data_filt_lipid %>% 
  filter(Var1 %in% oct & Var2 %in% oct) %>% 
  mutate(month = "october")



# Joining data

bnti_final_lipid <- rbind(bnti_data_may, bnti_data_june, bnti_data_july_1, bnti_data_july_2, bnti_data_july_3,
                        bnti_data_aug, bnti_data_sep, bnti_data_oct)


bnti_final_lipid$month <- factor(bnti_final_lipid$month, levels = c("may", "june", "july_1", "july_2",
                                                                "july_3", "august", "september","october"))

################## SUGARS

##sugarss

# Matching the order of RCBC to bNTI results
char.rcbc_sugars = char.rcbc_sugars[row.names(char_sugars), colnames(char_sugars), drop = F]


# Setting the RCBC diagonals to NA
diag(char.rcbc_sugars) = NA


# Reflecting null matrices
char_sugars[upper.tri(char_sugars)] = t(char_sugars)[upper.tri(char_sugars)]


# Removing significant bNTI results from the RCBC results
char.rcbc_sugars[abs(char_sugars) > 2] = NA


# ##################################### #
#### Preprocessing Data for plotting ####
# ##################################### #

# Melting data
char_sugars = melt(as.matrix(char_sugars))

char.rcbc_sugars = melt(as.matrix(char.rcbc_sugars))


# Removing null values
char_sugars = char_sugars[!is.na(char_sugars$value),]

char.rcbc_sugars = char.rcbc_sugars[!is.na(char.rcbc_sugars$value),]


###### Ecological processes per whole system ###############

# Determining ecological processes fractionation 
## This is for the whole system

## Renaming columns

bnti_data_sugars <- char_sugars %>% 
  rename(BNTI = value)

rcbc_data_sugars <- char.rcbc_sugars %>% 
  rename(RCBC = value)

## Joining BNTI and RCBC datasets

data_full_sugars <- full_join(bnti_data_sugars, rcbc_data_sugars)


## removing duplicated comparisons

temp_data <- data_full_sugars

for(i in 1:nrow(temp_data)){
  temp_data$comb[i] <- paste0(sort(c(as.character(temp_data$Var1[i]), 
                                     as.character(temp_data$Var2[i]))), 
                              collapse = ';')
}

data_filt_sugars <- temp_data %>% 
  group_by(comb) %>% 
  summarise(BNTI = mean(BNTI),
            RCBC = mean(RCBC, na.rm = TRUE)) %>% 
  separate(comb, into = c('Var1', 'Var2'), sep = ';')


## Determining ecological processes fractionation

eco_sugars <- data_filt_sugars %>% 
  mutate(ecological = case_when(BNTI <= -2 ~ "Homogeneous selection",
                                BNTI >= 2 ~ "Variable selection",
                                abs(BNTI) <2 & RCBC < -0.95 ~ "Homogenizing dispersal",
                                abs(BNTI) <2 & RCBC > 0.95 ~ "Dispersal limitation",
                                abs(BNTI) <2 & abs(RCBC) <= 0.95 ~ "Undominated")) %>% 
  dplyr::count(ecological) %>% 
  mutate(percentage = n/sum(n)*100)


#### Plot by sampling - month  ####

# Boxplots by like samples (example: May vs May)

bnti_data_may <- data_filt_sugars %>% 
  filter(Var1 %in% may & Var2 %in% may) %>% 
  mutate(month = "may")

bnti_data_june <- data_filt_sugars %>% 
  filter(Var1 %in% june & Var2 %in% june) %>% 
  mutate(month = "june")

bnti_data_july_1 <- data_filt_sugars %>% 
  filter(Var1 %in% july_1 & Var2 %in% july_1) %>% 
  mutate(month = "july_1")

bnti_data_july_2 <- data_filt_sugars %>% 
  filter(Var1 %in% july_2 & Var2 %in% july_2) %>% 
  mutate(month = "july_2")

bnti_data_july_3 <- data_filt_sugars %>% 
  filter(Var1 %in% july_3 & Var2 %in% july_3) %>% 
  mutate(month = "july_3")

bnti_data_aug <- data_filt_sugars %>% 
  filter(Var1 %in% aug & Var2 %in% aug) %>% 
  mutate(month = "august")


bnti_data_sep <- data_filt_sugars %>% 
  filter(Var1 %in% sept & Var2 %in% sept) %>% 
  mutate(month = "september")

bnti_data_oct <- data_filt_sugars %>% 
  filter(Var1 %in% oct & Var2 %in% oct) %>% 
  mutate(month = "october")



# Joining data

bnti_final_sugars <- rbind(bnti_data_may, bnti_data_june, bnti_data_july_1, bnti_data_july_2, bnti_data_july_3,
                          bnti_data_aug, bnti_data_sep, bnti_data_oct)


bnti_final_sugars$month <- factor(bnti_final_sugars$month, levels = c("may", "june", "july_1", "july_2",
                                                                    "july_3", "august", "september","october"))


########################## TLC


##tlcs

# Matching the order of RCBC to bNTI results
char.rcbc_tlc = char.rcbc_tlc[row.names(char_tlc), colnames(char_tlc), drop = F]


# Setting the RCBC diagonals to NA
diag(char.rcbc_tlc) = NA


# Reflecting null matrices
char_tlc[upper.tri(char_tlc)] = t(char_tlc)[upper.tri(char_tlc)]


# Removing significant bNTI results from the RCBC results
char.rcbc_tlc[abs(char_tlc) > 2] = NA


# ##################################### #
#### Preprocessing Data for plotting ####
# ##################################### #

# Melting data
char_tlc = melt(as.matrix(char_tlc))

char.rcbc_tlc = melt(as.matrix(char.rcbc_tlc))


# Removing null values
char_tlc = char_tlc[!is.na(char_tlc$value),]

char.rcbc_tlc = char.rcbc_tlc[!is.na(char.rcbc_tlc$value),]


###### Ecological processes per whole system ###############

# Determining ecological processes fractionation 
## This is for the whole system

## Renaming columns

bnti_data_tlc <- char_tlc %>% 
  rename(BNTI = value)

rcbc_data_tlc <- char.rcbc_tlc %>% 
  rename(RCBC = value)

## Joining BNTI and RCBC datasets

data_full_tlc <- full_join(bnti_data_tlc, rcbc_data_tlc)


## removing duplicated comparisons

temp_data <- data_full_tlc

for(i in 1:nrow(temp_data)){
  temp_data$comb[i] <- paste0(sort(c(as.character(temp_data$Var1[i]), 
                                     as.character(temp_data$Var2[i]))), 
                              collapse = ';')
}

data_filt_tlc <- temp_data %>% 
  group_by(comb) %>% 
  summarise(BNTI = mean(BNTI),
            RCBC = mean(RCBC, na.rm = TRUE)) %>% 
  separate(comb, into = c('Var1', 'Var2'), sep = ';')


## Determining ecological processes fractionation

eco_tlc <- data_filt_tlc %>% 
  mutate(ecological = case_when(BNTI <= -2 ~ "Homogeneous selection",
                                BNTI >= 2 ~ "Variable selection",
                                abs(BNTI) <2 & RCBC < -0.95 ~ "Homogenizing dispersal",
                                abs(BNTI) <2 & RCBC > 0.95 ~ "Dispersal limitation",
                                abs(BNTI) <2 & abs(RCBC) <= 0.95 ~ "Undominated")) %>% 
  dplyr::count(ecological) %>% 
  mutate(percentage = n/sum(n)*100)


#### Plot by sampling - month  ####

# Boxplots by like samples (example: May vs May)

bnti_data_may <- data_filt_tlc %>% 
  filter(Var1 %in% may & Var2 %in% may) %>% 
  mutate(month = "may")

bnti_data_june <- data_filt_tlc %>% 
  filter(Var1 %in% june & Var2 %in% june) %>% 
  mutate(month = "june")

bnti_data_july_1 <- data_filt_tlc %>% 
  filter(Var1 %in% july_1 & Var2 %in% july_1) %>% 
  mutate(month = "july_1")

bnti_data_july_2 <- data_filt_tlc %>% 
  filter(Var1 %in% july_2 & Var2 %in% july_2) %>% 
  mutate(month = "july_2")

bnti_data_july_3 <- data_filt_tlc %>% 
  filter(Var1 %in% july_3 & Var2 %in% july_3) %>% 
  mutate(month = "july_3")

bnti_data_aug <- data_filt_tlc %>% 
  filter(Var1 %in% aug & Var2 %in% aug) %>% 
  mutate(month = "august")


bnti_data_sep <- data_filt_tlc %>% 
  filter(Var1 %in% sept & Var2 %in% sept) %>% 
  mutate(month = "september")

bnti_data_oct <- data_filt_tlc %>% 
  filter(Var1 %in% oct & Var2 %in% oct) %>% 
  mutate(month = "october")



# Joining data

bnti_final_tlc <- rbind(bnti_data_may, bnti_data_june, bnti_data_july_1, bnti_data_july_2, bnti_data_july_3,
                           bnti_data_aug, bnti_data_sep, bnti_data_oct)


bnti_final_tlc$month <- factor(bnti_final_tlc$month, levels = c("may", "june", "july_1", "july_2",
                                                                      "july_3", "august", "september","october"))



######### Ecological processes by month

## Microbiome

eco_month_microbe <- bnti_final_microbe %>% 
  mutate(ecological = case_when(BNTI <= -2 ~ "Homogeneous selection",
                                BNTI >= 2 ~ "Variable selection",
                                abs(BNTI) <2 & RCBC < -0.95 ~ "Homogenizing dispersal",
                                abs(BNTI) <2 & RCBC > 0.95 ~ "Dispersal limitation",
                                abs(BNTI) <2 & abs(RCBC) <= 0.95 ~ "Undominated")) %>% 
  group_by(month) %>% 
  dplyr::count(ecological) %>% 
  mutate(percentage = n/sum(n)*100) %>% 
  mutate(Analysis = "Microbiome") %>% 
  rename(set = month)


## Metabolome - bulk


eco_month_icr <- bnti_final_icr %>% 
  mutate(ecological = case_when(BNTI <= -2 ~ "Homogeneous selection",
                                BNTI >= 2 ~ "Variable selection",
                                abs(BNTI) <2 & RCBC < -0.95 ~ "Homogenizing dispersal",
                                abs(BNTI) <2 & RCBC > 0.95 ~ "Dispersal limitation",
                                abs(BNTI) <2 & abs(RCBC) <= 0.95 ~ "Undominated")) %>% 
  group_by(month) %>% 
  dplyr::count(ecological) %>% 
  mutate(percentage = n/sum(n)*100) %>% 
  mutate(Analysis = "Metabolome") %>% 
  rename(set = month)

## Metabolome - lipid

eco_month_lipid <- bnti_final_lipid %>% 
  mutate(ecological = case_when(BNTI <= -2 ~ "Homogeneous selection",
                                BNTI >= 2 ~ "Variable selection",
                                abs(BNTI) <2 & RCBC < -0.95 ~ "Homogenizing dispersal",
                                abs(BNTI) <2 & RCBC > 0.95 ~ "Dispersal limitation",
                                abs(BNTI) <2 & abs(RCBC) <= 0.95 ~ "Undominated")) %>% 
  group_by(month) %>% 
  dplyr::count(ecological) %>% 
  mutate(percentage = n/sum(n)*100) %>% 
  mutate(Analysis = "Lipids") %>% 
  rename(set = month)

## Metabolome - sugars

eco_month_sugars <- bnti_final_sugars %>% 
  mutate(ecological = case_when(BNTI <= -2 ~ "Homogeneous selection",
                                BNTI >= 2 ~ "Variable selection",
                                abs(BNTI) <2 & RCBC < -0.95 ~ "Homogenizing dispersal",
                                abs(BNTI) <2 & RCBC > 0.95 ~ "Dispersal limitation",
                                abs(BNTI) <2 & abs(RCBC) <= 0.95 ~ "Undominated")) %>% 
  group_by(month) %>% 
  dplyr::count(ecological) %>% 
  mutate(percentage = n/sum(n)*100) %>% 
  mutate(Analysis = "Sugars") %>% 
  rename(set = month)


## Metabolome - tlc

eco_month_tlc <- bnti_final_tlc %>% 
  mutate(ecological = case_when(BNTI <= -2 ~ "Homogeneous selection",
                                BNTI >= 2 ~ "Variable selection",
                                abs(BNTI) <2 & RCBC < -0.95 ~ "Homogenizing dispersal",
                                abs(BNTI) <2 & RCBC > 0.95 ~ "Dispersal limitation",
                                abs(BNTI) <2 & abs(RCBC) <= 0.95 ~ "Undominated")) %>% 
  group_by(month) %>% 
  dplyr::count(ecological) %>% 
  mutate(percentage = n/sum(n)*100) %>% 
  mutate(Analysis = "TLC") %>% 
  rename(set = month)


########################## FINAL ################################################3
## Bar plot -ecological processes ALL

eco_final <- rbind(eco_month_microbe, eco_month_icr, eco_month_lipid, eco_month_sugars, eco_month_tlc) 


eco_colors <- set_names(c('#e76f51', '#f4a261', '#264653', '#2a9d8f', '#f1faee'),
                        nm = c('Variable selection', 'Homogeneous selection',
                               'Dispersal limitation', 'Homogenizing dispersal', 
                               'Undominated'))

ecobar_month_final <- eco_final %>% 
  mutate(Analysis = factor(Analysis, levels = c('Microbiome', 'Metabolome',
                                                'Lipids', 'Sugars', 'TLC')),
         ecological = factor(ecological, levels = c('Variable selection', 'Homogeneous selection',
                                                    'Dispersal limitation', 'Homogenizing dispersal', 
                                                    'Undominated'))) %>% 
ggplot(aes(y = set, x = percentage, fill = ecological))+
  geom_col(color = 'black', alpha = 0.8) +
  labs(x = expression(bold("Assembly patterns (%)")),
       fill = "Ecological processes")+
  scale_y_discrete(limits = rev) +
  scale_fill_manual(values = eco_colors)+
  facet_grid(cols = vars(Analysis)) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(color = 'black', size = 7),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1),
        legend.key.spacing = unit(.5, 'lines'),
        legend.key.size = unit(.5, 'lines'),
        legend.box.spacing = unit(.1, 'lines'),
        legend.position = 'bottom',
        legend.title.position = 'top')

ecobar_month_final


######################## Expression ##################################

## Metadata

metadata <- read_csv('data/expression/metadata.csv') 

metadata_ready <- metadata %>% 
  mutate(Month = factor(Month, levels = c("May", "July_1", "July_2", "July_3",
                                          "August", "October", "May_22"))) %>% 
  rename(sampleid = SampleID)

## Taxonomy

taxonomy <- read_csv("data/expression/mag_taxonomy_denovo.csv") %>% 
  rename(mag_id = bin)


## geTMM

geTMM <- read_csv("data/expression/geTMM_normalized_matrix_2024.csv")

geTMM_long <- geTMM %>% 
  pivot_longer(!Geneid, names_to = "sampleid", values_to = "geTMM")

## Expression calculated reaction level


stress <- read_csv("data/expression/expression_per_sample_traits_10-31.csv")

def_stress <- read_csv("data/expression/traits_definitions_update_10_9.csv")

stress_ready <- stress %>% 
  left_join(def_stress, by = c('module_id' = 'def_id'))


## Motility 


motility <- stress_ready %>% 
  filter(Category_1 == 'Motility')

## Plotting Motility 


motility_ready <- motility %>% 
  left_join(metadata, by = 'SampleID') %>% 
  group_by(module_id, bin, Name, Month, Category_1, Category_2) %>% 
  summarise(expression_mean = mean(expression)) %>% 
  filter(expression_mean > 0) %>% 
  left_join(taxonomy, by = c('bin' = 'mag_id')) %>% 
  ungroup() %>% 
  arrange(Phylum, bin) %>% 
  mutate(order = n():1)


## Heatmap flagella assembly

color_strips <- strip_nested(
  background_y = elem_list_rect(fill = c("white")))


motility_heat <- motility_ready %>% 
  mutate(Month = factor(Month, levels = c("May", "July_1", "July_2", "July_3",
                                          "August", "October", "May_22")),
         log_getmm = log10(expression_mean)) %>%
  ggplot(aes(x = Month,
             y = fct_reorder(bin, order),
             fill = log_getmm)) +
  geom_tile(color = 'white')+
  scale_fill_distiller(palette = 'PuBuGn', direction = 1)+
  facet_grid2(rows = vars(Category_2, Name),
              scale = 'free',
              space = 'free',
              strip = color_strips) +
  labs(fill = 'Log10(geTMM)') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.text.y = element_text(angle = 0, size = 6),
        title = element_text(face = 'bold', hjust = 0.5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(.8, 'lines'),
        legend.position = 'bottom',
        legend.text = element_text(size = 6),
        legend.box.spacing = unit(.01, 'lines'),
        legend.title.position = 'top',
        legend.title = element_text(size = 7))

motility_heat


## Phylum colors

phylum_colors <- set_names(get_palette('Paired', 18),
                           nm = levels(taxonomy$Phylum))

void_strips <- strip_nested(
  background_y = element_blank(),
  text_y = element_blank(),
  size = 'variable')

phylum_tile <-  motility_ready %>% 
  ggplot(aes(x = '1',
             y = fct_reorder(bin, order),
             fill = Phylum)) +
  geom_tile(color = 'white')+
  facet_grid2(rows = vars(Category_2, Name),
              scale = 'free',
              space = 'free',
              strip = void_strips) +
  labs(y = 'MAG',
       x = 'A') +
  scale_fill_manual(values = phylum_colors) +
  scale_x_discrete(expand = c(0, 0)) +
  theme_void() +
  theme(axis.title.y = element_text(size = 7, angle = 90, face = 'bold'),
        legend.key.spacing = unit(.5, 'lines'),
        legend.key.size = unit(.5, 'lines'),
        legend.box.spacing = unit(.1, 'lines'))

phylum_tile

heatmp_ready <- phylum_tile + motility_heat +
  plot_layout(guides = 'collect', nrow = 1,
              widths = c(.05, 1))+
  theme(legend.position = 'bottom')

heatmp_ready

## Final figure

layout <- "
AB
CC
DD
"


final_fig3 <- new_plot_all + new_plot_all_icr + ecobar_month_final + 
  free(heatmp_ready) +
  plot_layout(design = layout, heights = c(1, 1, 4)) +
  plot_annotation(tag_levels = list(c('A', 'B', 'C', 'D', ''))) &
  theme(text = element_text(size = 7))

final_fig3


ggsave('output_figures/Figure_4.svg', final_fig3,
       dpi = 300, height = 210, width = 180, units = "mm")

ggsave('output_figures/Figure_4.png', final_fig3,
       dpi = 300, height = 210, width = 180, units = "mm")


