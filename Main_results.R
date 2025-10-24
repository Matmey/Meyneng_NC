
##########################
# 0. Packages & options ##
##########################
# Install missing packages if needed (uncomment)
# required_pkgs <- c("phyloseq","stringr","dplyr","tidyr","ggplot2","patchwork",
#                    "ggpubr","forcats","vegan","ade4","adespatial","ComplexHeatmap",
#                    "circlize","metagMisc","lme4","sjPlot","tibble","rio","ggridges",
#                    "ggparty","partykit","magrittr","readr")
# missing <- setdiff(required_pkgs, rownames(installed.packages()))
# if(length(missing)) install.packages(missing)

library(phyloseq)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(forcats)
library(vegan)
library(ade4)
library(adespatial)
library(ComplexHeatmap)
library(circlize)
library(tibble)
library(readr)
library(partykit) 
library(metagMisc)
library(tseries)
library(rio)
library(lme4)
library(sjPlot)
library(ggridges)

options(stringsAsFactors = FALSE)
theme_set(theme_pubr(base_size = 10))

##########################
# 1. Data import #########
##########################
setwd(".")
if (dir.exists("Outputs/") == FALSE) {system("mkdir -p Outputs/")}

asv_tax_df <- read.table("rawdata/ASV_tax_table.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
metadata_df <- read.csv("rawdata/sample_metadata.csv", header = TRUE, sep = ";", stringsAsFactors = FALSE)
env_df <- read.csv("rawdata/data_env.csv", header = TRUE, dec = ",", sep = ";", stringsAsFactors = FALSE)
foram_df <- read.csv("rawdata/Foram_table.csv", header = TRUE, sep = ";", stringsAsFactors = FALSE)
rs_db_df <- read.csv("rawdata/database_RestSt.csv", header = TRUE, sep = ";", stringsAsFactors = FALSE)

##########################
# 2. Build phyloseq ######
##########################

# Assume asv_tax_df has first column ASV_ID and following columns are samples until taxonomy column at the end, as output of SAMBA workflow.
# Adjust if your file format differs.
# Find taxonomy column by name or pattern
tax_col_index <- grep("taxonomy", colnames(asv_tax_df), ignore.case = TRUE)
if(length(tax_col_index) == 0) stop("taxonomy column not found in ASV_tax_table.tsv")

# OTU table (counts)
otu_df_raw <- asv_tax_df[, -tax_col_index]
rownames(otu_df_raw) <- otu_df_raw[[1]]
otu_df_raw <- otu_df_raw[ , -1]
otu_mat <- as.matrix(otu_df_raw)
OTU_TABLE <- otu_table(otu_mat, taxa_are_rows = TRUE)

# TAX table
tax_df_raw <- asv_tax_df[, c(1, tax_col_index)]
colnames(tax_df_raw)[1] <- "ASV_ID"
# Split taxonomy into levels: up to 9
tax_split <- str_split_fixed(tax_df_raw[[2]], ";", 9)
tax_tab <- data.frame(ASV_ID = tax_df_raw$ASV_ID,
                      Domain = tax_split[,1],
                      Supergroup = tax_split[,2],
                      Division = tax_split[,3],
                      Subdivision = tax_split[,4],
                      Class = tax_split[,5],
                      Order = tax_split[,6],
                      Family = tax_split[,7],
                      Genus = tax_split[,8],
                      Species = tax_split[,9],
                      stringsAsFactors = FALSE)
rownames(tax_tab) <- tax_tab$ASV_ID
TAX_TABLE <- tax_table(as.matrix(tax_tab[,-1]))

# Sample data
rownames(metadata_df) <- metadata_df$Sample_ID
SAMPLE_DATA <- sample_data(metadata_df)

# Raw phyloseq
ps_raw <- phyloseq(OTU_TABLE, TAX_TABLE, SAMPLE_DATA)

##########################
# 3. Blank filtering #####
##########################

# identify blank samples by a column in sample_data (Depth_sampling == "Blc") used in original script
if("Depth_sampling" %in% colnames(sample_data(ps_raw))) {
  blanks_ps <- subset_samples(ps_raw, Depth_sampling == "Blc")
  blanks_ps <- prune_taxa(taxa_sums(blanks_ps) > 0, blanks_ps)
  asvs_in_blanks <- taxa_names(blanks_ps)
} else {
  asvs_in_blanks <- character(0)
  message("Warning: Depth_sampling column not present; skipping blank detection.")
}

##########################
# 4. Prepare analysis phyloseq objects
##########################

# 6.1 physeq_alpha: remove singletons/doubletons, remove blank ASVs, keep Eukaryota and filter off multicellular taxa
ps_alpha <- prune_taxa(taxa_sums(ps_raw) > 2, ps_raw)   # remove singletons & doubletons
if(length(asvs_in_blanks) > 0) ps_alpha <- prune_taxa(!(taxa_names(ps_alpha) %in% asvs_in_blanks), ps_alpha)

# Apply taxonomy filters (adjust exact strings/spaces if necessary)
ps_alpha <- subset_taxa(ps_alpha,
                        Domain == "Eukaryota" &
                          Subdivision != " Metazoa" &
                          Division != " Streptophyta" &
                          Class != " Ulvophyceae" )
ps_alpha <- subset_taxa(ps_alpha,
                        Class != " Florideophyceae" &
                          Genus != " Sargassum")
ps_alpha <- prune_samples(sample_sums(ps_alpha) > 0, ps_alpha)

# 6.2 physeq_rarefied: rarefy for comparative composition analyses
set.seed(12345)
# Define a reasonable sample.size: choose the minimum sample sum above zero, or use a fixed value (e.g. 40000) if available
ps_rarefied <- rarefy_even_depth(ps_raw, sample.size = 40000, trimOTUs = TRUE, replace = FALSE, rngseed = 12345)
if(length(asvs_in_blanks) > 0) ps_rarefied <- prune_taxa(!(taxa_names(ps_rarefied) %in% asvs_in_blanks), ps_rarefied)
ps_rarefied <- subset_taxa(ps_rarefied,
                           Domain == "Eukaryota" &
                             Subdivision != " Metazoa" &
                             Division != " Streptophyta" &
                             Class != " Ulvophyceae")
ps_rarefied <- subset_taxa(ps_rarefied,
                           Class != " Florideophyceae" &
                             Genus != " Sargassum")
ps_rarefied <- prune_samples(sample_sums(ps_rarefied) > 0, ps_rarefied)

# Save processed phyloseq objects if needed
#saveRDS(ps_alpha, file = "ps_alpha.rds")
#saveRDS(ps_rarefied, file = "ps_rarefied.rds")

##########################
# 5. Resting-stage subsets#
##########################
rs_asvs <- rs_db_df %>% filter(RestingStage == 1) %>% pull(1)
# rs_db may contain header row or an initial NA; remove NA and ensure ASV IDs match taxa names
rs_asvs <- rs_asvs[!is.na(rs_asvs) & rs_asvs != ""]

ps_rs_alpha <- prune_taxa(taxa_names(ps_alpha) %in% rs_asvs, ps_alpha)
ps_rs_rarefied <- prune_taxa(taxa_names(ps_rarefied) %in% rs_asvs, ps_rarefied)

#saveRDS(ps_rs_rarefied, "ps_rs_rarefied.rds")

##########################
# 6. Identification of sedimentary units
##########################
# Select relevant environmental variables 
env_subset <- env_df %>%
  select(
    Codification, Depth,
    Ni, Al
  )

#  Filter out rows with no sediment or metal data
env_subset <- env_subset %>%
  mutate(env_sum = rowSums(across(Ni:Al), na.rm = TRUE)) %>%
  filter(env_sum > 0) %>%
  select(-env_sum)

env_subset <- env_subset %>%
  mutate(Depth_neg = -Depth)   # easier for plotting or modeling

#  Conditional Inference Tree (ctree)
# Model depth partitioning based on Ni and Al concentration
env_ctree <- ctree(
  Depth_neg ~ Ni + Al,
  data = env_subset,
  control = ctree_control(mincriterion = 0.995)  # strict split criterion
)

# Print tree summary
print(env_ctree)

#  Plot the tree 
tree_plot <- plot(env_ctree, gp = gpar(fontsize = 10))

# Inspect samples matching a specific node rule
# (e.g. Unit IV)
env_subset_IV <- env_subset %>%
  filter(Al > 31846) %>%
  filter(Ni <= 1920)

print(env_subset_IV)

##########################
# 7. Figure 1: Environmental changes along core
##########################
# Prepare long format for selected variables
env_vars <- c("Clay","Sand","Ni","Al")
env_long <- env_df %>%
  select(Codification, Depth, all_of(env_vars)) %>%
  pivot_longer(cols = all_of(env_vars), names_to = "Variable", values_to = "Value") %>%
  mutate(Category = case_when(
    Variable %in% c("Clay", "Sand") ~ "Grain size",
    Variable == "Ni" ~ "Nickel",
    Variable == "Al" ~ "Aluminium",
    TRUE ~ "Other"
  )) %>%
  filter(!is.na(Value) & Value > 0) %>%
  mutate(Depth_neg = -Depth)

env_long$Category <- factor(env_long$Category, levels = c("Nickel","Aluminium","Grain size"))
env_long$Variable <- factor(env_long$Variable, levels = c("Ni","Al","Clay","Sand"))
color_variable <- c("#FF794D","#59BDC0","#BB593E","#D9AB37")

fig1_env_profile <- ggplot(env_long, aes(x = Depth_neg, y = Value, color = Variable, fill = Variable)) +
  geom_vline(xintercept = c(-46, -94, -138, -202), linetype = "dashed", color = "grey") + #sedimentary units established by ctree function (see extended data 1)
  facet_grid(~Category, scales = "free", space = "free_y") +
  geom_line(aes(group = Variable), size = 0.8, alpha = 0.6) +
  geom_point(size = 1) +
  scale_fill_manual(values = color_variable) + scale_color_manual(values = color_variable) +
  labs(x = "Depth (cm)", y = "Value") + coord_flip() +
  theme_pubr(base_size = 9)

# Source composition (UB, VS)
source_long <- env_df %>% select(Codification, Depth, UB, VS) %>%
  pivot_longer(cols = c("UB","VS"), names_to = "Source", values_to = "Value") %>%
  mutate(Depth_neg = -Depth) %>%
  filter(!is.na(Value) & Value > 0)
color_source <- c("#E75A1E","#2AADB6")

fig1_source <- ggplot(source_long, aes(x = Depth_neg, y = Value, fill = Source)) +
  geom_vline(xintercept = c(-46, -94, -138, -202), linetype = "dashed", color = "grey") +
  geom_area(alpha = 0.6, color = "white") +
  scale_fill_manual(values = color_source) +
  labs(y = "Value") + coord_flip() + theme_pubr(base_size = 9) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

fig1_combined <- fig1_env_profile + fig1_source + plot_layout(widths = c(3,1))
ggsave("Outputs/Figure1_environment_profiles.pdf", plot = fig1_combined,
       width = 180, height = 140, units = "mm", dpi = 600)

##########################
# 8. Figure 2: Alpha diversity (foraminiferal genera + ASV)
##########################

# --- Foraminifera summary by sedimentary units 
foram_df$Unit <- factor(foram_df$Unit, levels = c("V","IV","III","II","I"))

foram_summary <- foram_df %>%
  group_by(Unit) %>%
  summarise(mean_genus = mean(Sum_genus, na.rm = TRUE),
            sd_genus = sd(Sum_genus, na.rm = TRUE),
            min_date = min(date_median, na.rm = TRUE),
            max_date = max(date_median, na.rm = TRUE)) %>%
  mutate(med_date = (min_date + max_date)/2) %>% ungroup()

# Wilcoxon pairwise comparisons (successive parts)
pairwise_wilcox_foram <- list(
  V_vs_IV = wilcox.test(Sum_genus ~ Unit, data = subset(foram_df, Unit %in% c("V","IV"))),
  IV_vs_III = wilcox.test(Sum_genus ~ Unit, data = subset(foram_df, Unit %in% c("IV","III"))),
  III_vs_II = wilcox.test(Sum_genus ~ Unit, data = subset(foram_df, Unit %in% c("III","II"))),
  II_vs_I = wilcox.test(Sum_genus ~ Unit, data = subset(foram_df, Unit %in% c("II","I")))
)

# Foraminifera plot
fig2_foram_plot <- ggplot(foram_df, aes(y = Sum_genus, x = date_median)) +
  geom_vline(xintercept = c(1976, 1960, 1757, 1301), linetype = "dashed", color = "grey") +
  geom_ribbon(aes(ymin = 0, ymax = Sum_genus), alpha = 0.25) +
  geom_point(data = foram_summary, aes(x = med_date, y = mean_genus), inherit.aes = FALSE, size = 2) +
  geom_line(data = foram_summary, aes(x = med_date, y = mean_genus), inherit.aes = FALSE, size = 0.8) +
  geom_linerange(data = foram_summary, aes(x = med_date, ymin = mean_genus - sd_genus, ymax = mean_genus + sd_genus), inherit.aes = FALSE) +
  labs(x = "Median date", y = "Genus number") + coord_flip() + theme_pubr(base_size = 9)

# --- Microeukaryote alpha (Observed ASVs) using ps_alpha
alpha_richness_df <- estimate_richness(ps_alpha, measures = c("Observed"))
alpha_richness_df$SampleID <- rownames(alpha_richness_df)
# Merge with sample metadata (ensure columns exist)
meta_alpha <- as.data.frame(sample_data(ps_alpha))
meta_alpha$SampleID <- rownames(meta_alpha)
alpha_meta_df <- left_join(alpha_richness_df, meta_alpha, by = "SampleID")
# Ensure numeric columns
alpha_meta_df$Date_median <- as.numeric(alpha_meta_df$Date_median)
alpha_meta_df$Depth_reel <- as.numeric(alpha_meta_df$Depth_reel)
alpha_meta_df$alpha_observed <- as.numeric(alpha_meta_df$Observed)

alpha_summary_by_part <- alpha_meta_df %>%
  group_by(Unit) %>%
  summarise(mean_alpha = mean(alpha_observed, na.rm = TRUE),
            sd_alpha = sd(alpha_observed, na.rm = TRUE),
            min_date = min(Date_median, na.rm = TRUE),
            max_date = max(Date_median, na.rm = TRUE)) %>%
  mutate(med_date = c(1999,1968,1858,1520,1221)) %>% ungroup()

# pairwise Wilcoxon for alpha
pairwise_wilcox_alpha <- list(
  V_vs_IV = wilcox.test(alpha_observed ~ Unit, data = subset(alpha_meta_df, Unit %in% c("V","IV"))),
  IV_vs_III = wilcox.test(alpha_observed ~ Unit, data = subset(alpha_meta_df, Unit %in% c("IV","III"))),
  III_vs_II = wilcox.test(alpha_observed ~ Unit, data = subset(alpha_meta_df, Unit %in% c("III","II"))),
  II_vs_I = wilcox.test(alpha_observed ~ Unit, data = subset(alpha_meta_df, Unit %in% c("II","I")))
)

fig2_alpha_plot <- ggplot(alpha_meta_df, aes(y = alpha_observed, x = Date_median)) +
  geom_vline(xintercept = c(1976, 1960, 1757, 1301), linetype = "dashed", color = "grey") +
  geom_ribbon(aes(ymin = 0, ymax = alpha_observed), alpha = 0.25) +
  geom_point(data = alpha_summary_by_part, aes(x = med_date, y = mean_alpha), inherit.aes = FALSE, size = 2) +
  geom_line(data = alpha_summary_by_part, aes(x = med_date, y = mean_alpha), inherit.aes = FALSE, size = 0.8) +
  geom_linerange(data = alpha_summary_by_part, aes(x = med_date, ymin = mean_alpha - sd_alpha, ymax = mean_alpha + sd_alpha), inherit.aes = FALSE) +
  labs(x = "Median date", y = "ASV number") + coord_flip() + theme_pubr(base_size = 9)

# --- Podani partitioning (moving window) using ps_alpha abundance table
otu_alpha_mat <- t(as.data.frame(otu_table(ps_alpha))) # samples as rows
otu_alpha_mat <- otu_alpha_mat[order(rownames(otu_alpha_mat), decreasing = TRUE), ]

window_size <- 2
df_podani <- data.frame()
for(i in 1:(nrow(otu_alpha_mat) - window_size)) {
  tmp <- otu_alpha_mat[i:(i + window_size), ]
  tmp_beta <- adespatial::beta.div.comp(tmp, coef = "J", quant = FALSE)
  tmp_beta <- data.frame(start = rownames(tmp)[1], end = rownames(tmp)[nrow(tmp)], tmp_beta$part) %>%
    tibble::rownames_to_column("component")
  df_podani <- dplyr::bind_rows(df_podani, tmp_beta)
}

# Merge with depth
df_podani<- df_podani %>%
  mutate(period = paste(start, end, sep = "-")) %>%
  mutate(Sample_ID = paste(end)) 

sample_depth <- select(alpha_meta_df, "Sample_ID","Depth_sampling", "Unit")

df_podani <- merge(df_podani, sample_depth, by = "Sample_ID")
df_podani$Depth_sampling <- as.numeric(df_podani$Depth_sampling)

# Plot data
comp_plot <- c("BDtotal","Repl","RichDif")

fig2_podani_plot <- ggplot(df_podani %>% filter(component %in% comp_plot), aes(x = -Depth_sampling, y = tmp_beta.part, color = component, fill = component)) +
  geom_vline(xintercept = c(-41,-93,-137,-202,-226), linetype = "dashed", color = "grey") +
  geom_line(size = 0.9, alpha = 0.6) +
  scale_color_manual(values = c("#616161", "#72B678", "#407B7B")) +
  geom_point(size = 1) +
  coord_flip() + theme_pubr(base_size = 9) + 
  labs(
    x = "Depth (cm)",
    y = "Podani's indices"
  )

# Combine Figure 2 panels
fig2_combined <- fig2_foram_plot + fig2_alpha_plot + fig2_podani_plot + plot_layout(ncol = 3)
ggsave(filename = "Outputs/Figure2_alpha_diversity.pdf", plot = fig2_combined,
       width = 200, height = 120, units = "mm", dpi = 600)

##########################
# 9. Figure 3: Rarefied dataset - richness & composition
##########################
# Ensure ps_rarefied is present and filtered similarly (exclude multicellular)
ps_rarefied_prot <- subset_taxa(ps_rarefied,
                                Domain == "Eukaryota" &
                                  Subdivision != " Metazoa" &
                                  Division != " Streptophyta" &
                                  Class != " Ulvophyceae")
ps_rarefied_prot <- subset_taxa(ps_rarefied_prot,
                                Class != " Florideophyceae" &
                                  Genus != " Sargassum")
ps_rarefied_prot <- prune_samples(sample_sums(ps_rarefied_prot) > 0, ps_rarefied_prot)

# Resting-stage and unknown-genus subsets (using rs_db_df)
asv_rs_all <- rs_db_df %>% filter(RestingStage == 1) %>% pull(1) %>% na.omit() #based on functional database of Ramond et al., 2018 and additional literature when necessary
asv_unk <- rs_db_df %>% filter(is.na(RestingStage) & Unk_genus == 1) %>% pull(1) %>% na.omit() # i.e., no genus assignment

ps_rs_rar <- prune_taxa(taxa_names(ps_rarefied_prot) %in% asv_rs_all, ps_rarefied_prot)
ps_rs_unk_rar <- prune_taxa(taxa_names(ps_rarefied_prot) %in% asv_unk, ps_rarefied_prot)

# Richness observed
rich_rs <- estimate_richness(ps_rs_rar, measures = "Observed")
rich_rs$SampleID <- rownames(rich_rs)
colnames(rich_rs)[1] <- "Observed_rs"

reads_rs <- colSums(otu_table(ps_rs_rar))
reads_rs_df <- data.frame(SampleID = names(reads_rs), TotalReads_rs = reads_rs, stringsAsFactors = FALSE)

rich_tot <- estimate_richness(ps_rarefied_prot, measures = "Observed")
rich_tot$SampleID <- rownames(rich_tot)
colnames(rich_tot)[1] <- "Observed_tot"
reads_tot <- colSums(otu_table(ps_rarefied_prot))
reads_tot_df <- data.frame(SampleID = names(reads_tot), TotalReads_tot = reads_tot, stringsAsFactors = FALSE)

rich_unk <- estimate_richness(ps_rs_unk_rar, measures = "Observed")
rich_unk$SampleID <- rownames(rich_unk); colnames(rich_unk)[1] <- "Observed_rs_unk"
reads_unk <- colSums(otu_table(ps_rs_unk_rar))
reads_unk_df <- data.frame(SampleID = names(reads_unk), TotalReads_rs_unk = reads_unk, stringsAsFactors = FALSE)

# Merge richness/read summaries
rich_merge <- purrr::reduce(list(rich_tot[,c("Observed_tot","SampleID")], rich_rs[,c("Observed_rs","SampleID")],
                          rich_unk[,c("Observed_rs_unk","SampleID")], reads_tot_df, reads_rs_df, reads_unk_df),
                     function(x,y) full_join(x,y, by = "SampleID"))

# Add sample depth
colnames(sample_depth) <- c("SampleID","Depth_sampling","Unit")
rich_merge <- merge(rich_merge, sample_depth, by = "SampleID")

# Compute percentages
rich_merge <- rich_merge %>%
  mutate(Observed_add = coalesce(Observed_rs,0) + coalesce(Observed_rs_unk,0),
         TotalReads_add = coalesce(TotalReads_rs,0) + coalesce(TotalReads_rs_unk,0),
         perc_asv_rs = Observed_rs_unk / Observed_tot * 100,
         perc_read_rs = TotalReads_rs_unk / TotalReads_tot * 100)

# Visual panels (observed ASV numbers)
rich_merge <- rich_merge %>%
  arrange(Depth_sampling) %>%  # tri par profondeur croissante
  mutate(SampleID = factor(SampleID, levels = SampleID))  # fixe l’ordre pour ggplot

rich_long <- data.frame(SampleID = factor(rep(rich_merge$SampleID, 3), levels = rich_merge$SampleID),
                        Value = c(rich_merge$Observed_tot, rich_merge$Observed_rs, rich_merge$Observed_add),
                        Depth = rep(rich_merge$Depth_sampling, 3),
                        Type = factor(rep(c("Observed_tot","Observed_rs","Observed_rs_unk_add"), each = nrow(rich_merge))))
rich_long$Lower <- 0
rich_long$Depth <- as.numeric(rich_long$Depth)
rich_long$Type <- factor(rich_long$Type, levels = c("Observed_tot","Observed_rs_unk_add","Observed_rs"))

fig3_obs_plot <- ggplot(rich_long, aes(y = Value, x = -Depth, group = Type, color = Type, fill = Type)) +
  geom_ribbon(aes(ymin = Lower, ymax = Value), alpha = 1) +
  geom_line(alpha = 1, size = 1) +
  geom_vline(xintercept = c(-41,-93,-137), linetype = "dashed", color = "grey") +
  ylab("Number of ASV") + xlab("Depth (cm)") + coord_flip() + theme_pubr(base_size = 9) +
  scale_color_manual(values = c("black","#636363","#69C9C9")) +
  scale_fill_manual(values = c("white","#9C9C9C","#95DCDC")) + 
  theme(legend.position = "none")

# Reads panel
reads_long <- data.frame(SampleID = factor(rep(rich_merge$SampleID, 3), levels = rich_merge$SampleID),
                         Value = c(rich_merge$TotalReads_tot, rich_merge$TotalReads_rs, rich_merge$TotalReads_add),
                         Depth = rep(rich_merge$Depth_sampling, 3),
                         Type = factor(rep(c("TotalReads_tot","TotalReads_rs","TotalRead_rs_unk_add"), each = nrow(rich_merge))))
reads_long$Lower <- 0
reads_long$Depth <- as.numeric(reads_long$Depth)
reads_long$Type <- factor(reads_long$Type, levels = c("TotalReads_tot","TotalRead_rs_unk_add","TotalReads_rs"))

fig3_read_plot <- ggplot(reads_long, aes(y = Value, x = -Depth, group = Type, color = Type, fill = Type)) +
  geom_ribbon(aes(ymin = Lower, ymax = Value), alpha = 1) +
  geom_line(alpha = 1, size = 1) +
  geom_vline(xintercept = c(-41,-93,-137), linetype = "dashed", color = "grey") +
  ylab("Number of reads") + coord_flip() + theme_pubr(base_size = 9) +
  scale_color_manual(values = c("black","#636363","#69C9C9")) +
  scale_fill_manual(values = c("white","#9C9C9C","#95DCDC"))+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Composition of resting-stage taxa (relative abundance)
ps_rs_rel <- transform_sample_counts(ps_rs_rar, function(x) x * 100 / sum(x))
rs_melt <- psmelt(ps_rs_rel)

rs_melt <- rs_melt %>%
  mutate(Subdivision = if_else(Subdivision == "", paste0(.$Division, "_X"), .$Subdivision))%>%
  mutate(Subdivision = if_else(Subdivision == "_X", " Eukaryota_X", .$Subdivision))

rs_agg <- rs_melt %>%
  group_by(Sample, Depth_sampling, Subdivision, Part) %>%
  summarise(Abundance_per_Subdivision = sum(Abundance), .groups = "drop")

# Collapse rare subdivisions into "Others"
rare_subs <- rs_agg %>% 
  group_by(Subdivision) %>% 
  summarise(maxA = max(Abundance_per_Subdivision, na.rm = TRUE)) %>% 
  filter(maxA < 1) %>% pull(Subdivision)

rs_agg <- rs_agg %>% 
  mutate(Subdivision = fct_collapse(as.factor(Subdivision), Others = rare_subs))

# Add filler rows to get nice depth scale (only visualisation purpose)
filler_samples <- data.frame(Sample = c("A40", "A44", "A48", "A52", "A57", "A62", "A66", "A74", "A80", "A86", 
                                        "A108", "A112", "A117", "A122", "A126", "A130", "A134", "A138", "A142",
                                        "A150", "A156", "A162", "A166", "A170", "A174"),
                             Depth_sampling = c(40,44,48,52,57,62,66,74,80,86,108,112,117,122,126,130,134,138,142,150,156,162,166,170,174),
                             Subdivision = "", Part = "None", Abundance_per_Subdivision = 100,
                             stringsAsFactors = FALSE)
filler_samples$Depth_sampling <- as.numeric(filler_samples$Depth_sampling)
rs_agg$Depth_sampling <- as.numeric(rs_agg$Depth_sampling)
rs_plot_df <- bind_rows(rs_agg, filler_samples)

# Factor ordering by depth values (adjust levels if needed)
unique_depths <- sort(unique(as.numeric(rs_plot_df$Depth_sampling)))
rs_plot_df$Depth_sampling <- factor(rs_plot_df$Depth_sampling, levels = as.character(unique_depths))

# Color vector - adjust as desired
col_vector_rs <- c("#F3EB51","#EAE289","#FFA54B","#931C0F","#824600","#F5C17D",
                   "#6373DB","#01579E","#1F91B6","#7DD9F8","#BAF2E9","#FAACEB",
                   "#5A9451","#6BDB57","#979797","white")

rs_plot_df$Subdivision <- factor(rs_plot_df$Subdivision, 
                                   levels=c(" Evosea_X"," Tubulinea_X",
                                            " Chlorophyta_X",
                                            " Haptophyta_X",
                                            " Fungi"," Ichthyosporea",
                                            " Apicomplexa"," Chrompodellids"," Ciliophora"," Dinoflagellata"," Perkinsea",
                                            " Cercozoa",
                                            " Bigyra"," Gyrista",
                                            "Others",""
                                   ))

fig3_compo_rs_plot <- ggplot(rs_plot_df, aes(x = Abundance_per_Subdivision, y = fct_rev(Depth_sampling), fill = Subdivision)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = col_vector_rs, guide = guide_legend(ncol = 1)) +
  theme_minimal(base_size = 9) + labs(x = "Relative abundance (%)", y = "Depth") + theme(panel.grid = element_blank())

# For the complete community composition (rarefied)
ps_complet_rel <- transform_sample_counts(ps_rarefied, function(x) x * 100 / sum(x))
complet_melt <- psmelt(ps_complet_rel)

complet_melt <- complet_melt %>%
  mutate(Subdivision = if_else(Subdivision == "", paste0(.$Division, "_X"), .$Subdivision))%>%
  mutate(Subdivision = if_else(Subdivision == "_X", " Eukaryota_X", .$Subdivision))

complet_agg <- complet_melt %>%
  group_by(Sample, Depth_sampling, Subdivision, Part) %>%
  summarise(Abundance_per_Subdivision = sum(Abundance), .groups = "drop")

rare_subs2 <- complet_agg %>% 
  group_by(Subdivision) %>% 
  summarise(maxA = max(Abundance_per_Subdivision, na.rm = TRUE)) %>% 
  filter(maxA < 2) %>% 
  pull(Subdivision)

complet_agg <- complet_agg %>% 
  mutate(Subdivision = fct_collapse(as.factor(Subdivision), Others = rare_subs2))

filler_big <- data.frame(Sample = c("A40", "A44", "A48", "A52", "A57", "A62", "A66", "A74", "A80", "A86", 
                                        "A108", "A112", "A117", "A122", "A126", "A130", "A134", "A138", "A142",
                                        "A150", "A156", "A162", "A166", "A170", "A174"),
                             Depth_sampling = c(40,44,48,52,57,62,66,74,80,86,108,112,117,122,126,130,134,138,142,150,156,162,166,170,174),
                             Subdivision = "", Part = "None", Abundance_per_Subdivision = 100,
                             stringsAsFactors = FALSE)

complet_agg$Depth_sampling <- as.numeric(complet_agg$Depth_sampling)
filler_big$Depth_sampling <- as.numeric(filler_big$Depth_sampling)
complet_plot_df <- bind_rows(complet_agg, filler_big)

unique_depths2 <- sort(unique(as.numeric(complet_plot_df$Depth_sampling)))
complet_plot_df$Depth_sampling <- factor(as.character(complet_plot_df$Depth_sampling), levels = as.character(unique_depths2))

col_vector_complet <- c("#F3EB51","#FFA54B","#DD5757","#931C0F","#740B97","#824600","#F5C17D", "#573000",
                        "#6373DB","#01579E","#1F91B6","#7DD9F8","#023667","#FAACEB","#5A9451","#6BDB57",
                        "#096306","#000000","#979797","white")

complet_plot_df$Subdivision <- factor(complet_plot_df$Subdivision, 
                                    levels=c(" Evosea_X",
                                             " Chlorophyta_X",
                                             " Cryptophyta_X"," Haptophyta_X",
                                             " Breviatea_X",
                                             " Fungi"," Ichthyosporea"," Opisthokonta_X",
                                             " Apicomplexa"," Chrompodellids"," Ciliophora"," Dinoflagellata"," Alveolata_X",
                                             " Cercozoa",
                                             " Bigyra"," Gyrista"," Stramenopiles_X",
                                             " Eukaryota_X","Others",""
                                    ))

fig3_compo_rar_plot <- ggplot(complet_plot_df, aes(x = Abundance_per_Subdivision, y = fct_rev(Depth_sampling), fill = Subdivision)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = col_vector_complet, guide = guide_legend(ncol = 1)) +
  theme_minimal(base_size = 9) + labs(x = "Relative abundance (%)", y = "Depth") + theme(panel.grid = element_blank())

fig3_combined <- fig3_obs_plot + fig3_read_plot + fig3_compo_rs_plot + fig3_compo_rar_plot + plot_layout(ncol = 4,nrow=1, heights = c(1,1,1,1))
ggsave(filename = "Outputs/Figure3_rarefied_composition.pdf", plot = fig3_combined,
       width = 250, height = 200, units = "mm", dpi = 600)

##########################
# 10. Figure 4: Genus dynamics
##########################
# Prepare foram data
foram_df[is.na(foram_df)] <- 0
rownames(foram_df) <- foram_df$Sample
foram_div <- foram_df[, c(9:33)]   # select relevant columns (species presence/absence)

# Identify non-random distribution patterns using runs test 
# Convert abundance to presence/absence
abundance_data <- foram_div

# Identify ubiquitous taxa present in all samples
genus_all <- colnames(abundance_data)[colSums(abundance_data) == nrow(abundance_data)]

# Remove these columns temporarily to test only the variable taxa
abundance_data2 <- abundance_data[, !(colnames(abundance_data) %in% genus_all)]

# Identify taxa with significantly non-random distributions (p < 0.05)
significant_species <- c()
for (i in 1:ncol(abundance_data2)) {
  species_data <- abundance_data2[, i]
  test_result <- runs.test(as.factor(species_data > 0), alternative = "less")
  if (test_result$p.value < 0.05) {
    significant_species <- c(significant_species, colnames(abundance_data2)[i])
  }
}

# Subset abundance matrix to keep significant + ubiquitous species
subset_abundance_data <- abundance_data[, colnames(abundance_data) %in% c(significant_species, genus_all)]
rownames(subset_abundance_data) <- as.factor(foram_df$Depth)  # use depth as row labels

# Convert to long format for ggplot heatmap 
abundance_data_long <- subset_abundance_data %>%
  as.data.frame() %>%
  rownames_to_column("Row") %>%
  pivot_longer(cols = -Row, names_to = "Column", values_to = "Value")

# Order rows (depths) and columns (taxa)
abundance_data_long$Row <- factor(abundance_data_long$Row, levels = rev(rownames(subset_abundance_data)))
abundance_data_long$Column <- factor(abundance_data_long$Column, levels = c(
  "Ammonia","Quinqueloculina", "Nonion", "Haynesina",
  "Triloculina","Spiroloculina","Wiesnerella","Eggerella",
  "Amphistegina","Planorbulinella"
))

# Color scale
col_fun <- colorRamp2(c(-1, 0, 1), c("#662141", "#FBFBFB", "#524657"))

# Plot heatmap of selected foraminiferal genera
fig4_foram <- ggplot(abundance_data_long, aes(x = Column, y = Row, fill = Value)) +
  geom_tile(color = "white", linewidth = 0.1) +
  scale_fill_gradientn(colors = col_fun(seq(min(abundance_data_long$Value), max(abundance_data_long$Value), length.out = 100))) +
  labs(fill = "Foraminifera presence") +
  theme_minimal(base_size = 6) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title = element_blank(),
    legend.position = "none"
  )

### Barpolts of filtered microeukaryote genera

# Collapse ASVs to genus level
physeq_genus <- tax_glom(ps_rs_rar, taxrank = "Genus")

# Relative abundance per sample
physeq_genus_rel <- transform_sample_counts(physeq_genus, function(x) x * 100 / sum(x))

# Filter genera: present in >4 samples 
physeq_genus_n4 <- filter_taxa(
  physeq_genus_rel,
  function(x) sum(x > 0) > ((1/56*4)*length(x)),
  TRUE
)

# Convert to long format for tests
genus_n4 <- psmelt(physeq_genus_n4)

# Aggregate data by genus, depth, and stratigraphic unit (Part2)
df_genus <- aggregate(
  Abundance ~ Sample + Genus + OTU + Unit + Depth_sampling,
  data = genus_n4,
  sum
)

#--------------- Wilcoxon non-parametric analysis

genus.all <- levels(factor(df_genus$OTU))
tmp <- df_genus[-which(df_genus$Abundance==0),] # remove null values
genus.non.null <- table(tmp$OTU)
table(tmp$OTU, tmp$Unit) # non null samples per genus and periods

# We want to test differences between periods III and IV (shovel vs. pristine); II and III (tractor vs. shovel); I and II (Resilience)

df.test <- NULL

for(i in 1:length(genus.all)){
  
  # Select a genius
  genus.i <- df_genus[which(df_genus$OTU == genus.all[i]),]
  
  # define the periods
  genus.i.period.1 <- genus.i[which(genus.i$Unit == "I"),]
  genus.i.period.2 <- genus.i[which(genus.i$Unit == "II"),]
  genus.i.period.3 <- genus.i[which(genus.i$Unit == "III"),]
  genus.i.period.4 <- genus.i[which(genus.i$Unit == "IV"),]
  
  # Proportion of non null samples in each period
  prop.non.null <- tapply(genus.i$Abundance, genus.i$Unit, function(x)length(which(x==0))/length(x))
  
  # Proportion of non null samples in all sample
  prop.non.null.tot <- length(which(genus.i$Abundance==0))/length(genus.i$Abundance)
  
  # Average abundance in each period 
  avg.abd <- tapply(genus.i$Abundance, genus.i$Unit, mean)
  
  #---- test 1 (shovel vs. pristine)
  test1 <- wilcox.test(genus.i.period.3$Abundance, genus.i.period.4$Abundance)
  rel.chge1 <- ((avg.abd[3]+1) - (avg.abd[4]+1))/(avg.abd[4]+1)*100
  pval1 <- test1$p.value
  
  #---- test 2 (tractor vs. shovel)
  test2 <- wilcox.test(genus.i.period.2$Abundance, genus.i.period.3$Abundance)
  rel.chge2 <- ((avg.abd[2]+1) - (avg.abd[3]+1))/(avg.abd[3]+1)*100
  pval2 <- test2$p.value
  
  #---- test 3 (Resilience)
  test3 <- wilcox.test(genus.i.period.1$Abundance, genus.i.period.2$Abundance)
  rel.chge3 <- ((avg.abd[1]+1) - (avg.abd[2]+1))/(avg.abd[2]+1)*100
  pval3 <- test3$p.value
  
  # Store results
  out <- c(genus.all[i], prop.non.null, prop.non.null.tot, avg.abd, rel.chge1, pval1, rel.chge2, pval2, rel.chge3, pval2)
  df.test <- rbind(df.test, out)
  
}

df.test <- as.data.frame(df.test)
colnames(df.test) <- c("genus", 
                       "prop0.I", "prop0.II", "prop0.III", "prop0.IV","prop0.tot",
                       "avgAbd.I", "avgAbd.II", "avgAbd.III", "avgAbd.IV",
                       "rel.chge.test1", "pval.test1", "rel.chge.test2", "pval.test2", "rel.chge.test3", "pval.test3")

df.test[,-1] <- apply(df.test[,-1], 2, function(x)round(as.numeric(x), 3))
rownames(df.test) <- NULL

#--- Some summary statistics
tmp <- as.matrix(df.test[c("rel.chge.test1", "rel.chge.test2", "rel.chge.test3")])
rownames(tmp) <- df.test$genus

prop.signif <- apply(df.test[c("pval.test1", "pval.test2", "pval.test3")], 2, function(x)length(which(x<0.1))/length(x)) # proportion of species with significant effects
prop.decline <- apply(tmp, 2, function(x)length(which(x<0))/length(x)) # proportion of species declining regardless of significance
prop.increase <- apply(tmp, 2, function(x)length(which(x>0))/length(x)) # proportion of species increasing regardless of significance
med.change <- apply(tmp, 2, median, na.rm=T) # median of relative changes
avg.change <- apply(tmp, 2, mean, na.rm=T) # average of relative changes

#--- add Genus 
taxo <- select(genus_n4, "Genus", "OTU")
colnames(df.test)[1] <- "OTU"
df.test2 <- merge(df.test,taxo,  by = "OTU") 
df.test2 <- df.test2 %>% distinct(OTU, .keep_all = TRUE)

# --- Filter significant genera ---
df_filtered <- df.test2 %>%
  filter(
    pval.test1 < 0.05 | pval.test2 < 0.05 | pval.test3 < 0.05 | prop0.tot < 0.25)

selection_genera <- df_filtered$OTU

genera_filtered <- df_genus %>%
  filter(OTU %in% selection_genera)

genera_filtered$Depth_sampling <- as.numeric(genera_filtered$Depth_sampling) # 23 genera

# Visualisation of some of these genera (based on theur temporal distribution and bibliographic information) 
genera_figure <- c(" Desmodesmus"," Amoebogregarina", " Colpodella",
                " Fungi_XXXX", " Selenidium1"," Paracercomonas", " Scrippsiella",
                " Coelastrum"," Paramicrosporidium")

genera_figure <- genera_filtered %>% filter(Genus %in% genera_figure)

genera_figure$Genus <- factor(genera_figure$Genus,
                                 levels = c(" Paracercomonas", " Amoebogregarina"," Colpodella"," Scrippsiella",
                                            " Paramicrosporidium"," Selenidium1"," Desmodesmus"," Coelastrum",
                                            " Fungi_XXXX" ))

col_vector_genus <- c("#FAACEB","#6376DB","#01579E","#7DD9F8",
                      "#824600","#6376DB","#FFA54B","#FFA54B",
                      "#824600")

fig4_genus <- ggplot(genera_figure, aes(x = -Depth_sampling, y = Abundance, fill = Genus)) +
  geom_vline(aes(xintercept = -41), linetype = "dashed", color = "grey", size = 0.5) +
  geom_vline(aes(xintercept = -93), linetype = "dashed", color = "grey", size = 0.5) +
  geom_vline(aes(xintercept = -137), linetype = "dashed", color = "grey", size = 0.5) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = col_vector_genus, guide = guide_legend(ncol = 1)) +
  geom_point(
    data = subset(genera_figure, Abundance == 0),
    aes(x = -Depth_sampling, y = 0),  # y = 0 to highlight absence in the sample
    shape = 21,   
    size = 1,
    fill = "white"
  ) +
  facet_grid(. ~ Genus, scales = "free_x") +
  theme_pubr(base_size = 6) +
  labs(
    x = "Depth (cm)",
    y = "Relative Abundance among resting stage taxa"
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    strip.text.y = element_text(angle = 0, size = 10),
    legend.position = "none",
    panel.grid.minor.y = element_blank()
  ) +
  coord_flip()

#Panel of Nickel concentration along the core 

# Define a helper function to linearly interpolate missing values (to create a continuous visualisation)
interpolate_na <- function(df, variables) {
  df <- df %>% 
    mutate(across(
      all_of(variables),
      ~ approx(
        Depth[!is.na(.)], 
        .[!is.na(.)], 
        xout = Depth, 
        method = "linear", 
        rule = 2
      )$y
    )) %>%
    ungroup()
  
  return(df)
}

# --- Apply interpolation to selected variables ---
variables_to_interpolate <- "Ni" #possible to select several variables
df_interpolated <- interpolate_na(env_df, variables_to_interpolate)

# --- Prepare depth factor for plotting (reverse order for vertical axis) ---
rownames(df_interpolated) <- gsub("-", "", df_interpolated$Codification)
df_interpolated$Depth_fact <- factor(
  as.factor(df_interpolated$Depth),
  levels = rev(unique(df_interpolated$Depth))
)

# --- Plot Nickel concentration as a vertical gradient ---
fig4__nickel <- ggplot(df_interpolated, aes(x = Depth_fact, y = 1, fill = Ni)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("#023436", "#FFC65C", "#FF5416")) +
  labs(x = "Depth (cm)", y = NULL, fill = "Nickel (mg kg -1)") +
  coord_flip() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

fig4_combined <- fig4__nickel + fig4_foram + fig4_genus + plot_layout(ncol = 3,nrow=1, width = c(0.2,2,3))
ggsave(filename = "Outputs/Figure4_genus.pdf", plot = fig4_combined,
       width = 200, height = 300, units = "mm", dpi = 600)


