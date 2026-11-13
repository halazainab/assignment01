###1 - OVERVIEW----

#####Question----
#H ow does the transposable element differ between Trichinella spiralis and Trichinella murrelli in terms of overall repeat abundance, major element class composition, and element length distribution?

####Data Access----
# This code file should run on everyones system assuming they have downloaded the required packages and are in the latest software version. 

# Please make sure to download datasets from this NCBI link for T. Spirallis and T. Murrelli and save with the appropriate name and data/ directory:
#T. Spirallis: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_001447595.1/ 
#T. Murrelli: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_001447425.1/ 

####Installation---- 
#Package installation & load libraries + Data exploration for both species (T. spiralis & T. murrelli)

install.packages("tidyverse")
install.packages("permute")
install.packages("patchwork")

rm(list = ls())

library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
library(ggplot2)
library(dplyr)
library(stringr)
library(patchwork)

###2 - DATA EXPLORATION---- 

# Data wrangling. We will add a species label to each dataset and then combine them into one master table.

# Date of Download: October 23, 2026

#####trichinella spiralis----
spiralis_raw <- read_csv("../data/trichinella_spiralis_te.csv")

# check class
class(spiralis_raw)

# see variables 
names(spiralis_raw)

# see summary
summary(spiralis_raw)

# check dimensions
dim(spiralis_raw)

#####trichinella murelli----
murrelli_raw <- read_csv("../data/trichinella murrelli_te.csv")

# check class
class(murrelli_raw)

# see variables 
names(murrelli_raw)

# see summary
summary(murrelli_raw)

# check dimensions
dim(murrelli_raw)

###3 - DATA CLEANING---- 
#####Combination----

spiralis_clean <- spiralis_raw |>
  select(qry_id, qry_start, qry_end, qry_width, matching_class, perc_div, sw_score) |>
  mutate(species = "T. spiralis")

murrelli_clean <- murrelli_raw |>
  select(qry_id, qry_start, qry_end, qry_width, matching_class, perc_div, sw_score) |>
  mutate(species = "T. murrelli")

df_tstm <- bind_rows(spiralis_clean, murrelli_clean)

####Broad Classification----
# Group detailed repeat annotations into broad, interpretable classes using str_detect() from stingr.

unique(df_tstm$matching_class)

df_tstm <- df_tstm |>
  mutate(
    broad_class = case_when(
      str_detect(matching_class, "DNA") ~ "DNA Transposon",
      str_detect(matching_class, "LINE") ~ "LINE",
      str_detect(matching_class, "SINE") ~ "SINE",
      str_detect(matching_class, "LTR") ~ "LTR Retrotransposon",
      str_detect(matching_class, "Simple_repeat") ~ "Simple Repeat",
      str_detect(matching_class, "simple_repeat") ~ "Simple Repeat",
      str_detect(matching_class, "Satellite") ~ "Satellite",
      str_detect(matching_class, "PLE") ~ "PLE",
      str_detect(matching_class, "Helitron") ~ "RC/Helitron",
      str_detect(matching_class, "Low_complexity") ~ "Low Complexity",
      str_detect(matching_class, "Retroposon") ~ "Other Retrotransposon",
      TRUE ~ "Other/Unknown"))

#check unique groups in broad_class
unique(df_tstm$broad_class)

####Non-TE Filtering----

# Filter out non-TEs (decided to keep satellites as a good non transposable element comparison cohort for discussion)

combined_te <- df_tstm |>
  filter(!broad_class %in% c("Simple Repeat", "Low Complexity", "Other/Unknown")) |>
  mutate(broad_class = fct_drop(broad_class)) |>
  mutate(broad_class = fct_infreq(broad_class))

df_tstmclean <- combined_te |>
  mutate(df_telength = qry_end - qry_start +1)

df_tstmclean |>
  select(qry_start, qry_end, df_telength) |>
  head()

###4 - SUMMARY STATS----

# summary statistics for results slides

####Metric 1----

#Overall repeat statistics per species

overall_te_summary <- combined_te |>
  group_by(species) |>
  summarise(
    total_te_count = n(),
    median_length_bp = median(qry_width),
    mean_length_bp = mean(qry_width),
    total_te_bp = sum(qry_width)
  )

cat("\n--- Transposable Element Summary per Species ---\n")
print(overall_te_summary)

####Metric 2----

#Abundance and percentage per repeat class

te_class_summary <- combined_te |>
  count(species, broad_class) |>
  group_by(species) |>
  mutate(percentage = (n / sum(n)) * 100) |>
  ungroup()

cat("\n--- Transposable Element Class Proportions ---\n")
print(te_class_summary)

###5 - DATA VISUALIZATION----

####QrylengthXspecies----
#figure 1: Element Length Distribution (EDA / Quality Check)

fig1 <- ggplot(combined_te, aes(x = species, y = qry_width, fill = species)) + geom_boxplot(outlier.alpha = 0.15) +
  scale_y_log10() + 
  scale_fill_manual(values = c("T. spiralis" = "#2b8cbe", "T. murrelli" = "#e6550d")) + labs(
    title = "Comparison of Transposable Element Lengths",
    x = "Species",
    y = "TE Length in bp (log10 scale)"
  ) + 
  theme_classic() +
  theme(legend.position = "none")

#print to Plots viewer and save directly to output/
print(fig1)
ggsave("../output/figure1_length_distribution.png", plot = fig1, width = 6, height = 4.5)

####Class composition----
#figure 2: Major Repeat Class Composition (Main Comparative Result)

fig2 <- ggplot(te_class_summary, aes(x = broad_class, y = percentage, fill = species)) + geom_col(position = position_dodge(width = 0.8), width = 0.7) + scale_fill_manual(values = c("T. spiralis" = "#2b8cbe", "T. murrelli" = "#e6550d")) + labs(
  title = "Proportion of Transposable Element Classes",
  x = "Transposable Element Class",
  y = "Percentage of Total TEs (%)",
  fill = "Species"
) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "top"
  )

print(fig2)
ggsave("../output/figure2_class_composition.png", plot = fig2, width = 7, height = 5)

####Contig repeat density----

#step 8: figure 3: Contig-Level Repeat Density & Wilcoxon Test (Novelty Component)
#count number of repeats per contig using qry_id

contig_te_counts <- combined_te |>
  count(species, qry_id)

#perform Wilcoxon rank-sum test

wilcox_te_res <- wilcox.test(n ~ species, data = contig_te_counts)
p_val_label <- paste0("Wilcoxon rank-sum p = ", format.pval(wilcox_te_res$p.value, digits = 3))

fig3 <- ggplot(contig_te_counts, aes(x = species, y = n, fill = species)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.2) +
  scale_y_log10() +
  scale_fill_manual(values = c("T. spiralis" = "#2b8cbe", "T. murrelli" = "#e6550d")) +
  annotate("text", x = 1.5, y = max(contig_te_counts$n) * 0.8, label = p_val_label, fontface = "italic") +
  labs(
    title = "Transposable Element Density per Contig",
    subtitle = "Comparing spatial TE clustering across assembly contigs",
    x = "Species",
    y = "Number of TEs per Contig (log10 scale)"
  ) + 
  theme_classic() +
  theme(legend.position = "none")

print(fig3)
ggsave("../output/figure3_contig_density.png", plot = fig3, width = 6, height = 4.5)

cat("\n--- Statistical Test on TE Contig Density ---\n")
print(wilcox_te_res)

###

df_percdiv <- aggregate(x = df_tstmclean, by = list(MatchClass = df_tstmclean$broad_class), FUN = function(x) mean(x, na.rm = TRUE))
summary(df_percdiv)

####Divergence----
# consider comparing and confirming the evolutionary relationship between the two: At a rate of approximately 5% divergence from consensus (Kimura divergence). Under 5% is considered very low divergence. This will address question into which of the two species is more closely resembling the family consensus. Investigation into divergence will shed light into how alike they are (percentage conservation from the same ancestor/reference), whereas lowest divergence will allow for higher resolution into where the divergence starts to occur. 

dfOverDiv <- combined_te |>
  filter(!is.na(matching_class)) |>
  filter(perc_div >= 5) |>
  group_by(species, matching_class) |>
  summarize(avgPD = mean(perc_div, na.rm = TRUE)) |>
  arrange(desc(avgPD)) |>
  print()

####Conservation----
# now consider the lists of least divergence and look for similarity

dfUnderDiv <- combined_te |>
  filter(!is.na(matching_class)) |>
  filter(perc_div <= 5) |>
  group_by(species, matching_class) |>
  summarize(avgPD = mean(perc_div, na.rm = TRUE)) |> ungroup() |>
  arrange(avgPD) |>
  print()

#####Jitter Plot / Test----

ggplot(dfOverDiv, aes(x = avgPD, y = matching_class)) +
  geom_jitter(aes(color = species), width = 0., height = 0, alpha = 0.7, size = 3) +
  scale_color_manual(values = c("T. spiralis" = "blue", "T. murrelli" = "red")) +
  labs(
    x = "Average Divergence",
    y = "Repeat Class",
    color = "Species"
  ) +
  theme_minimal()

# did not enjoy the visualization (sparse + no real direction)

#####Dumbbell Plots----
# attempting dumbbell for better presentation 

spec_colors <- c(
  "T. spiralis" = "green", 
  "T. murrelli" = "darkgreen"
)

ggAvgDiv <- ggplot(dfOverDiv, aes(x = avgPD,y = matching_class, color = species)) +
  geom_line(
    aes(group = matching_class),
    color = "green",
    linewidth = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(values = spec_colors) +
  labs(
    x = "Average Percent Divergence (%)", 
    y = "Repeat Class", 
    color = "Species"
  ) + 
  theme_minimal() +
  labs(title = "Most Diverged Classes\nT.spiralis vs. T. murelli (Averaged)")

ggAvgDiv
# ggsave(
#   filename = "../output/AvgDiv.png",
#   plot = ggAvgDiv,
#   width = 8,
#   height = 8,
#   units = "in",
#   dpi = 300
# )

ggAvgCons <- ggplot(dfUnderDiv, aes(x = avgPD,y = matching_class, color = species)) +
  geom_line(
    aes(group = matching_class),
    color = "green",
    linewidth = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(values = spec_colors) +
  labs(
    x = "Average Percent Divergence (%)", 
    y = "Repeat Class", 
    color = "Species"
  ) + 
  theme_minimal() +
  labs(title = "Most Highly Conserved Classes\nT.spiralis vs. T. murelli (Averaged)")
ggAvgCons

AvgDivergence_plot <- ggAvgDiv + ggAvgCons +
  plot_layout(ncol = 2)

AvgDivergence_plot

ggsave(
  filename = "../output/fig_2_AvgDiv.png",
  plot = AvgDivergence_plot,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)

#####Box plots----
dfOverDiv2 <- combined_te |>
  filter(!is.na(matching_class)) |>
  filter(perc_div >= 5) |>
  group_by(species, matching_class) |>
  arrange(desc(perc_div)) |>
  print()

dfUnderDiv2 <- combined_te |>
  filter(!is.na(matching_class)) |>
  filter(perc_div <= 5) |>
  group_by(species, matching_class) |>
  arrange(desc(perc_div)) |>
  print()

ggDiv <- ggplot() +
  geom_boxplot(data = dfOverDiv2, aes(x = matching_class, y = perc_div, color = species)) +
  labs(
    title = "Repeat Divergence by class\n in T. spiralis vs. T.murelli",
    x = "Repeat Class",
    y = "Divergence (%)", 
    color = "Species"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggCons <- ggplot() +
  geom_boxplot(data = dfUnderDiv2, aes(x = matching_class, y = perc_div, color = species)) +
  labs(
    title = "Repeat Conservation by class\n in T. spiralis vs. T.murelli",
    x = "Repeat Class",
    y = "Divergence (%)", 
    color = "Species"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

Divergence_plot <- ggDiv + ggCons +
  plot_layout(ncol = 2)

Divergence_plot
ggsave(
  filename = "../output/fig_3_Div.png",
  plot = Divergence_plot,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)

#####Multi-Variable Bubble Plots----
# what if we look at sequence length x percent convergence? 

mv_combined <- ggplot(data = df_tstmclean, aes(x = perc_div, y  = df_telength, size = sw_score, fill = broad_class, shape = species)) +
  geom_point(alpha=0.5, color = "black") +
  scale_shape_manual(values = c(
    "T. spiralis" = 21,
    "T. murrelli" = 24)) +
  
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = 21, 
        color = "black", 
        alpha = 1
      )
    )
  ) +
  labs(
    fill = "Repeat Class",
    shape = "Species",
    size = "SW_score") +
  labs(
    title = "Multivariable Comparison between\n T. spiralis vs. T.murelli",
    x = "Percent Divergence (%)",
    y = "Sequence Length (bp)",
  ) +
  theme_minimal()

mv_edited <- ggplot(data = subset(df_tstmclean, broad_class != "LTR Retrotransposon"), aes(x = perc_div, y  = df_telength, size = sw_score, fill = broad_class, shape = species)) +
  geom_point(alpha=0.5, color = "black") +
  scale_shape_manual(values = c(
    "T. spiralis" = 21,
    "T. murrelli" = 24)) +
  
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = 21, 
        color = "black", 
        alpha = 1
      )
    )
  ) +
  labs(
    fill = "Repeat Class",
    shape = "Species",
    size = "SW_score") +
  labs(
    title = "Multivariable Comparison between\n T. spiralis vs. T.murelli (without LTR)",
    x = "Percent Divergence (%)",
    y = "Sequence Length (bp)",
  ) +
  theme_minimal()

Bubble_plot <- mv_combined + mv_edited +
  plot_layout(ncol = 2)
Bubble_plot

ggsave(
  filename = "../output/fig_45_Comp.png",
  plot = Bubble_plot,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)

ggsave(
  filename = "../output/fig_45_Compb.png",
  plot = Bubble_plot,
  width = 20,
  height = 8,
  units = "in",
  dpi = 300
)

