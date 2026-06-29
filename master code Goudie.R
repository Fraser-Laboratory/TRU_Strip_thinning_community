
library(tidyverse)
library(lme4)
library(lmerTest)
library(car)
library(writexl)
library(patchwork)
library(betapart)
library(vegan)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(multcomp)
library(janitor)
# -----------------------------
# 1. Clean column names
# -----------------------------
# This converts names to lower_case format and avoids duplicated naming issues.

final_table <- final_table %>%
  clean_names()


# -----------------------------
# 2. Basic structure of the dataset
# -----------------------------

dim(final_table)
names(final_table)
glimpse(final_table)
head(final_table, 10)


# -----------------------------
# 3. Identify metadata columns
# -----------------------------
# These are design/environmental/identifier variables.
# Adjust this list only if your dataset has additional metadata columns.

metadata_cols <- c(
  "block",
  "plot",
  "year",
  "trt",
  "treatment",
  "strip",
  "strip_id",
  "quadrat",
  "quad",
  "q",
  "distance",
  "distance_road",
  "elevation",
  "slope",
  "orientation"
)

# Keep only metadata columns that are actually present in final_table
metadata_cols <- intersect(metadata_cols, names(final_table))

# Species columns are all columns that are not metadata
species_cols <- setdiff(names(final_table), metadata_cols)


# -----------------------------
# 4. Check detected columns
# -----------------------------

metadata_cols
species_cols

length(metadata_cols)
length(species_cols)


# -----------------------------
# 5. Check key design variables
# -----------------------------

# Year
if ("year" %in% names(final_table)) {
  table(final_table$year, useNA = "ifany")
}

# Treatment
if ("trt" %in% names(final_table)) {
  table(final_table$trt, useNA = "ifany")
}

if ("treatment" %in% names(final_table)) {
  table(final_table$treatment, useNA = "ifany")
}

# Block
if ("block" %in% names(final_table)) {
  table(final_table$block, useNA = "ifany")
}

# Plot / strip identifiers
if ("plot" %in% names(final_table)) {
  table(final_table$plot, useNA = "ifany")
}

if ("strip_id" %in% names(final_table)) {
  table(final_table$strip_id, useNA = "ifany")
}

if ("quadrat" %in% names(final_table)) {
  table(final_table$quadrat, useNA = "ifany")
}


# -----------------------------
# 6. Summary of numeric metadata variables
# -----------------------------

numeric_metadata_cols <- metadata_cols[
  sapply(final_table[metadata_cols], is.numeric)
]

if (length(numeric_metadata_cols) > 0) {
  summary(final_table[numeric_metadata_cols])
}


# -----------------------------
# 7. Check species cover columns
# -----------------------------
# This assumes species columns contain cover values, usually 0 to 100.

summary(final_table[species_cols])


# Minimum and maximum cover value for each species
species_ranges <- final_table %>%
  summarise(
    across(
      all_of(species_cols),
      list(
        min = ~min(.x, na.rm = TRUE),
        max = ~max(.x, na.rm = TRUE)
      )
    )
  )

species_ranges


# -----------------------------
# 8. Species frequency check
# -----------------------------
# Number of quadrat-year observations where each species has cover > 0.

nonzero_counts <- final_table %>%
  summarise(
    across(
      all_of(species_cols),
      ~sum(.x > 0, na.rm = TRUE)
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "species",
    values_to = "nonzero_observations"
  ) %>%
  arrange(desc(nonzero_observations))

head(nonzero_counts, 20)
tail(nonzero_counts, 20)


# -----------------------------
# 9. Optional: check missing values
# -----------------------------

missing_summary <- final_table %>%
  summarise(
    across(
      everything(),
      ~sum(is.na(.x))
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "n_missing"
  ) %>%
  arrange(desc(n_missing))

missing_summary

############################################################################

# ============================================================
# GOUDIE PAPER - PART 1B
# Create sampling identifiers: strip_id, quadrat_id, and sample_id
# ============================================================

# Make sure treatment, year, and block are factors/clean variables
final_table <- final_table %>%
  mutate(
    block = factor(block),
    year = factor(year, levels = c("2018", "2019", "2020")),
    trt = factor(trt, levels = c("Control", "10 m", "15 m", "20 m")),
    orientation = factor(orientation)
  )

# ------------------------------------------------------------
# Create strip and quadrat identifiers
# ------------------------------------------------------------
# Assumption:
# Within each block × treatment × year combination, rows are ordered as:
# strip 1: quadrats 1–5
# strip 2: quadrats 1–5
#
# This matches the design:
# 3 blocks × 4 treatments × 2 strips × 5 quadrats × 3 years = 360 rows

final_table <- final_table %>%
  group_by(block, trt, year) %>%
  mutate(
    row_within_block_trt_year = row_number(),
    strip_number = ceiling(row_within_block_trt_year / 5),
    quadrat_number = ((row_within_block_trt_year - 1) %% 5) + 1
  ) %>%
  ungroup() %>%
  mutate(
    strip_id = paste(block, trt, strip_number, sep = "_"),
    quadrat_id = paste(block, trt, strip_number, quadrat_number, sep = "_"),
    sample_id = paste(quadrat_id, year, sep = "_")
  )

# Check that structure is correct
final_table %>%
  count(year, block, trt, strip_number)

final_table %>%
  count(year, block, trt, strip_number, quadrat_number)

# Check number of unique IDs
n_distinct(final_table$strip_id)
n_distinct(final_table$quadrat_id)
n_distinct(final_table$sample_id)

################################################## Beta diversity #################################################





############################################################
##   metadata and community matrix


env_cols <- c(
  "block", "year", "trt", "distance", "elevation",
  "slope", "orientation", "label", "sampling"
)

env_cols <- intersect(env_cols, names(vegetation_data))

metadata <- vegetation_data %>%
  dplyr::select(dplyr::all_of(env_cols)) %>%
  dplyr::mutate(
    block       = factor(block),
    year        = factor(year, levels = c("2018", "2019", "2020")),
    trt         = factor(trt),
    orientation = factor(orientation),
    label       = factor(label),
    sampling    = factor(sampling)
  )

species_cols <- setdiff(names(vegetation_data), env_cols)

comm_cover <- vegetation_data %>%
  dplyr::select(dplyr::all_of(species_cols))

comm_pa <- comm_cover %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::everything(),
      ~ ifelse(.x > 0, 1, 0)
    )
  )

comm_pa <- comm_pa[, colSums(comm_pa, na.rm = TRUE) > 0, drop = FALSE]

dim(comm_pa)
head(comm_pa[, 1:10])


## Jaccard beta-diversity partitioning


beta_jac <- betapart::beta.pair(
  comm_pa,
  index.family = "jaccard"
)

jac_total    <- beta_jac$beta.jac
jac_turnover <- beta_jac$beta.jtu
jac_nested   <- beta_jac$beta.jne

mean_jaccard <- data.frame(
  component = c("Total Jaccard", "Turnover", "Nestedness"),
  mean_value = c(
    mean(jac_total, na.rm = TRUE),
    mean(jac_turnover, na.rm = TRUE),
    mean(jac_nested, na.rm = TRUE)
  )
)

mean_jaccard


############################################################
## 4. Full PERMANOVA on Jaccard turnover
############################################################

set.seed(123)

permanova_turnover <- vegan::adonis2(
  jac_turnover ~ year + trt + distance + elevation + slope + orientation,
  data = metadata,
  strata = metadata$block,
  permutations = 999
)

permanova_turnover


############################################################
## 5. Table 1: Full PERMANOVA table
############################################################

table1_permanova <- as.data.frame(permanova_turnover) %>%
  tibble::rownames_to_column("Factor") %>%
  dplyr::rename(
    df = Df,
    R2 = R2,
    F_value = F,
    p_value = `Pr(>F)`
  ) %>%
  dplyr::mutate(
    Factor = dplyr::recode(
      Factor,
      "year" = "Year",
      "trt" = "Thinning",
      "distance" = "Distance",
      "elevation" = "Elevation",
      "slope" = "Slope",
      "orientation" = "Orientation",
      "Residual" = "Residual",
      "Total" = "Total"
    ),
    p_label = dplyr::case_when(
      is.na(p_value) ~ "",
      p_value < 0.001 ~ "<0.001",
      TRUE ~ as.character(round(p_value, 3))
    ),
    sig = dplyr::case_when(
      is.na(p_value) ~ "",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

table1_permanova


############################################################
## 6. Pairwise PERMANOVA by year
############################################################

run_pairwise_year_permanova <- function(year1, year2, metadata, comm_pa) {
  
  keep_rows <- metadata$year %in% c(year1, year2)
  
  meta_sub <- droplevels(metadata[keep_rows, , drop = FALSE])
  comm_sub <- comm_pa[keep_rows, , drop = FALSE]
  
  beta_sub <- betapart::beta.pair(
    comm_sub,
    index.family = "jaccard"
  )
  
  jac_sub <- beta_sub$beta.jtu
  
  res <- vegan::adonis2(
    jac_sub ~ year,
    data = meta_sub,
    strata = meta_sub$block,
    permutations = 999
  )
  
  list(
    years = paste(year1, year2, sep = " vs "),
    permanova = res
  )
}

set.seed(123)

pair_2018_2019 <- run_pairwise_year_permanova("2018", "2019", metadata, comm_pa)
pair_2018_2020 <- run_pairwise_year_permanova("2018", "2020", metadata, comm_pa)
pair_2019_2020 <- run_pairwise_year_permanova("2019", "2020", metadata, comm_pa)

pair_2018_2019$permanova
pair_2018_2020$permanova
pair_2019_2020$permanova


############################################################
## 7. Table 2: Pairwise PERMANOVA table
############################################################

extract_pairwise <- function(pair_object) {
  
  x <- as.data.frame(pair_object$permanova)
  
  data.frame(
    Comparison = pair_object$years,
    R2 = x$R2[1],
    F_value = x$F[1],
    p_value = x$`Pr(>F)`[1]
  )
}

table2_pairwise_permanova <- dplyr::bind_rows(
  extract_pairwise(pair_2018_2019),
  extract_pairwise(pair_2019_2020),
  extract_pairwise(pair_2018_2020)
) %>%
  dplyr::mutate(
    Comparison = factor(
      Comparison,
      levels = c("2018 vs 2019", "2019 vs 2020", "2018 vs 2020")
    ),
    p_label = dplyr::case_when(
      p_value < 0.001 ~ "<0.001",
      TRUE ~ as.character(round(p_value, 3))
    ),
    sig = dplyr::case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  dplyr::arrange(Comparison)

table2_pairwise_permanova


############################################################
## 8. Bootstrap SE for pairwise R2
############################################################

boot_R2_pair <- function(year1, year2, metadata, comm_pa, n_boot = 199) {
  
  keep_rows <- metadata$year %in% c(year1, year2)
  
  meta0 <- droplevels(metadata[keep_rows, , drop = FALSE])
  comm0 <- comm_pa[keep_rows, , drop = FALSE]
  
  idx1 <- which(meta0$year == year1)
  idx2 <- which(meta0$year == year2)
  
  r2_boot <- numeric(n_boot)
  
  for (b in seq_len(n_boot)) {
    
    idx <- c(
      sample(idx1, length(idx1), replace = TRUE),
      sample(idx2, length(idx2), replace = TRUE)
    )
    
    meta_b <- droplevels(meta0[idx, , drop = FALSE])
    comm_b <- comm0[idx, , drop = FALSE]
    
    beta_b <- betapart::beta.pair(
      comm_b,
      index.family = "jaccard"
    )
    
    jac_b <- beta_b$beta.jtu
    
    ad_b <- vegan::adonis2(
      jac_b ~ year,
      data = meta_b,
      strata = meta_b$block,
      permutations = 0
    )
    
    r2_boot[b] <- ad_b$R2[1]
  }
  
  r2_boot <- r2_boot[!is.na(r2_boot)]
  sd(r2_boot) * 100
}


############################################################
## 9. Data for Figure 2A: nestedness vs turnover
############################################################

turn_vals <- as.numeric(jac_turnover)
nest_vals <- as.numeric(jac_nested)

beta_box <- tibble::tibble(
  component = c(
    rep("Nestedness", length(nest_vals)),
    rep("Turnover", length(turn_vals))
  ),
  value = c(nest_vals, turn_vals)
) %>%
  tidyr::drop_na() %>%
  dplyr::mutate(
    component = factor(component, levels = c("Nestedness", "Turnover"))
  )


############################################################
## 10. Figure 2A: nestedness vs turnover
############################################################


fig2A <- ggplot(beta_box, aes(x = component, y = value, fill = component)) +
  geom_boxplot(
    outlier.size = 0.6,
    color = "black",
    width = 0.6
  ) +
  annotate(
    "text",
    x = 0.58,
    y = 0.95,
    label = "(A)",
    family = "Times New Roman",
    fontface = "bold",
    size = 5
  ) +
  scale_fill_manual(
    values = c(
      "Nestedness" = "#d9d9d9",
      "Turnover" = "#4d4d4d"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    expand = c(0, 0)
  ) +
  labs(
    x = "Jaccard β-component",
    y = "Jaccard dissimilarity β"
  ) +
  theme_classic(base_family = "Times New Roman") +
  theme(
    axis.text.x = element_text(size = 13, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 16, color = "black", face = "bold"),
    axis.title.y = element_text(size = 16, color = "black", face = "bold"),
    legend.position = "none",
    plot.margin = margin(5.5, 22, 5.5, 5.5)
  )

fig2A

############################################################
## 11. Data for Figure 2B: pairwise PERMANOVA R2
############################################################

set.seed(123)

pairwise_R2 <- table2_pairwise_permanova %>%
  dplyr::mutate(
    comparison = factor(
      as.character(Comparison),
      levels = c("2018 vs 2019", "2018 vs 2020", "2019 vs 2020")
    ),
    R2_percent = R2 * 100,
    SE_percent = c(
      boot_R2_pair("2018", "2019", metadata, comm_pa),
      boot_R2_pair("2018", "2020", metadata, comm_pa),
      boot_R2_pair("2019", "2020", metadata, comm_pa)
    ),
    fill_col = c("#4d4d4d", "#8c8c8c", "#d9d9d9"),
    star_y = R2_percent + SE_percent + 0.8
  )

pairwise_R2


############################################################
## 12. Figure 2B: pairwise PERMANOVA R2
############################################################

fig2B <- ggplot(pairwise_R2, aes(x = comparison, y = R2_percent, fill = fill_col)) +
  geom_col(
    color = "black",
    width = 0.6
  ) +
  geom_errorbar(
    aes(
      ymin = R2_percent - SE_percent,
      ymax = R2_percent + SE_percent
    ),
    width = 0.15,
    linewidth = 0.7
  ) +
  geom_text(
    aes(y = star_y, label = sig),
    family = "Times New Roman",
    fontface = "bold",
    size = 5
  ) +
  annotate(
    "text",
    x = 0.70,
    y = 19.0,
    label = "(B)",
    family = "Times New Roman",
    fontface = "bold",
    size = 5
  ) +
  scale_fill_identity() +
  scale_y_continuous(
    limits = c(0, 20),
    breaks = seq(0, 20, by = 5),
    expand = c(0, 0)
  ) +
  labs(
    x = "Year comparison",
    y = "Turnover variance explained (R², %)"
  ) +
  theme_classic(base_family = "Times New Roman") +
  theme(
    axis.text.x = element_text(size = 13, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 16, color = "black", face = "bold"),
    axis.title.y = element_text(size = 16, color = "black", face = "bold"),
    legend.position = "none",
    plot.margin = margin(5.5, 5.5, 5.5, 22)
  )

fig2B

############################################################
## 13. Combine Figure 2 with patchwork
############################################################

fig2 <- fig2A + fig2B +
  patchwork::plot_layout(ncol = 2)

fig2


############################################################
## 14. Save outputs
############################################################



ggsave(
  filename = "Figure2_beta_diversity.tiff",
  plot = fig2,
  dpi = 600,
  width = 10,
  height = 4.5,
  units = "in",
  compression = "lzw"
)

################################ Core vs transient richness ################################

############################################################
##  Define metadata and species columns
############################################################

env_cols <- c(
  "block", "year", "trt", "distance",
  "elevation", "slope", "orientation",
  "label", "sampling"
)

env_cols <- intersect(env_cols, names(vegetation_data))

species_cols <- setdiff(names(vegetation_data), env_cols)


############################################################
## 2. Convert species cover to presence-absence
############################################################

comm_pa_core <- vegetation_data %>%
  dplyr::select(dplyr::all_of(species_cols)) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::everything(),
      ~ ifelse(.x > 0, 1, 0)
    )
  )

comm_pa_core <- comm_pa_core[
  ,
  colSums(comm_pa_core, na.rm = TRUE) > 0,
  drop = FALSE
]

species_pa <- colnames(comm_pa_core)


############################################################
## 3. Calculate species temporal occupancy across years
############################################################

df_pa <- dplyr::bind_cols(
  vegetation_data %>%
    dplyr::select(dplyr::any_of(c("year", "label", "sampling"))),
  comm_pa_core
)

species_occupancy <- df_pa %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(species_pa),
      ~ max(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(species_pa),
      ~ sum(.x, na.rm = TRUE)
    )
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "species",
    values_to = "n_years"
  )


############################################################
## 4. Classify species as core or transient
############################################################

species_class <- species_occupancy %>%
  dplyr::mutate(
    class = dplyr::case_when(
      n_years >= 2 ~ "core",
      n_years == 1 ~ "transient",
      TRUE ~ NA_character_
    )
  )

table(species_class$class)


############################################################
## 5. Define core and transient species lists
############################################################

core_species <- species_class %>%
  dplyr::filter(class == "core") %>%
  dplyr::pull(species)

transient_species <- species_class %>%
  dplyr::filter(class == "transient") %>%
  dplyr::pull(species)


############################################################
## 6. Create quadrat-level richness table
############################################################

richness_df <- vegetation_data %>%
  dplyr::mutate(
    total_richness = rowSums(comm_pa_core, na.rm = TRUE),
    core_richness = rowSums(
      comm_pa_core[, core_species, drop = FALSE],
      na.rm = TRUE
    ),
    transient_richness = rowSums(
      comm_pa_core[, transient_species, drop = FALSE],
      na.rm = TRUE
    ),
    transient_prop = dplyr::if_else(
      total_richness > 0,
      transient_richness / total_richness,
      NA_real_
    )
  ) %>%
  dplyr::select(
    dplyr::any_of(c(
      "block", "label", "sampling", "year", "trt",
      "distance", "elevation", "slope", "orientation"
    )),
    total_richness,
    core_richness,
    transient_richness,
    transient_prop
  )


############################################################
## 7. Check richness outputs
############################################################

head(species_class)
head(richness_df)

summary(richness_df$total_richness)
summary(richness_df$core_richness)
summary(richness_df$transient_richness)
summary(richness_df$transient_prop)


############################################################
## 8. Prepare richness data for models
############################################################

richness_df2 <- richness_df %>%
  dplyr::mutate(
    block = factor(block),
    year = factor(year, levels = c("2018", "2019", "2020")),
    trt = factor(trt, levels = c("Control", "10 m", "15 m", "20 m")),
    orientation = factor(orientation)
  )


############################################################
## 9. Fit final richness models
############################################################

mod_total <- glmmTMB(
  total_richness ~ year * trt + (1 | block),
  family = compois,
  data = richness_df2
)

mod_core <- glmmTMB(
  core_richness ~ year * trt + (1 | block),
  family = compois,
  data = richness_df2
)

mod_transient <- glmmTMB(
  transient_richness ~ year * trt + (1 | block),
  family = poisson,
  data = richness_df2
)

mod_transient_prop <- glmmTMB(
  cbind(transient_richness, core_richness) ~ year * trt + (1 | block),
  family = binomial,
  data = richness_df2
)


############################################################
## 10. Type II Wald tests
############################################################

anova_total <- car::Anova(mod_total, type = 2)
anova_core <- car::Anova(mod_core, type = 2)
anova_transient <- car::Anova(mod_transient, type = 2)
anova_transient_prop <- car::Anova(mod_transient_prop, type = 2)

anova_total
anova_core
anova_transient
anova_transient_prop


############################################################
## 11. Table 3: richness model summary table
############################################################

extract_richness_anova <- function(anova_object, response_name) {
  
  anova_df <- as.data.frame(anova_object) %>%
    tibble::rownames_to_column("Effect")
  
  p_col <- grep("Pr\\(", names(anova_df), value = TRUE)
  
  anova_df %>%
    dplyr::transmute(
      Response = response_name,
      Effect = dplyr::recode(
        Effect,
        "year" = "Year",
        "trt" = "Thinning",
        "year:trt" = "Year × Thinning"
      ),
      Chisq = Chisq,
      df = Df,
      p_value = .data[[p_col]],
      p_label = dplyr::case_when(
        p_value < 0.001 ~ "<0.001",
        TRUE ~ as.character(round(p_value, 3))
      ),
      sig = dplyr::case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
}

table3_richness <- dplyr::bind_rows(
  extract_richness_anova(anova_total, "Total richness"),
  extract_richness_anova(anova_core, "Core richness"),
  extract_richness_anova(anova_transient, "Transient richness"),
  extract_richness_anova(anova_transient_prop, "Transient proportion")
)

table3_richness


############################################################
## 12. DHARMa diagnostic table
############################################################

get_dharma_tests <- function(model, model_name) {
  
  sim <- DHARMa::simulateResiduals(model)
  
  disp <- DHARMa::testDispersion(sim)
  zi <- DHARMa::testZeroInflation(sim)
  
  tibble::tibble(
    model = model_name,
    dispersion_statistic = as.numeric(disp$statistic),
    dispersion_p = disp$p.value,
    zero_inflation_statistic = as.numeric(zi$statistic),
    zero_inflation_p = zi$p.value
  )
}

dharma_table <- dplyr::bind_rows(
  get_dharma_tests(mod_total, "Total richness"),
  get_dharma_tests(mod_core, "Core richness"),
  get_dharma_tests(mod_transient, "Transient richness"),
  get_dharma_tests(mod_transient_prop, "Transient proportion")
)

dharma_table




################################ Figure 3 plotting ################################

############################################################
## Figure 3: Core and transient richness structure
## Input object: richness_df
## Output: fig3
############################################################


############################################################
## 1. Prepare data
############################################################

richness_df2 <- richness_df %>%
  dplyr::mutate(
    block = factor(block),
    year = factor(year, levels = c("2018", "2019", "2020")),
    trt = factor(trt, levels = c("10 m", "15 m", "20 m", "Control")),
    orientation = factor(orientation)
  )


############################################################
## 2. Summary by year
############################################################

summary_year <- richness_df2 %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(
    total_mean = mean(total_richness, na.rm = TRUE),
    total_se   = sd(total_richness, na.rm = TRUE) / sqrt(dplyr::n()),
    
    core_mean = mean(core_richness, na.rm = TRUE),
    core_se   = sd(core_richness, na.rm = TRUE) / sqrt(dplyr::n()),
    
    trans_mean = mean(transient_richness, na.rm = TRUE),
    trans_se   = sd(transient_richness, na.rm = TRUE) / sqrt(dplyr::n()),
    
    prop_mean = mean(transient_prop, na.rm = TRUE),
    prop_se   = sd(transient_prop, na.rm = TRUE) / sqrt(dplyr::n()),
    
    .groups = "drop"
  )


############################################################
## 3. Summary by thinning treatment
############################################################

summary_trt <- richness_df2 %>%
  dplyr::group_by(trt) %>%
  dplyr::summarise(
    total_mean = mean(total_richness, na.rm = TRUE),
    total_se   = sd(total_richness, na.rm = TRUE) / sqrt(dplyr::n()),
    
    core_mean = mean(core_richness, na.rm = TRUE),
    core_se   = sd(core_richness, na.rm = TRUE) / sqrt(dplyr::n()),
    
    trans_mean = mean(transient_richness, na.rm = TRUE),
    trans_se   = sd(transient_richness, na.rm = TRUE) / sqrt(dplyr::n()),
    
    prop_mean = mean(transient_prop, na.rm = TRUE),
    prop_se   = sd(transient_prop, na.rm = TRUE) / sqrt(dplyr::n()),
    
    .groups = "drop"
  )


############################################################
## 4. Panel data
############################################################

panel_A_data <- summary_year %>%
  dplyr::select(year, total_mean, total_se, core_mean, core_se) %>%
  tidyr::pivot_longer(
    cols = -year,
    names_to = c("metric", ".value"),
    names_pattern = "(total|core)_(mean|se)"
  ) %>%
  dplyr::mutate(
    metric = factor(
      metric,
      levels = c("total", "core"),
      labels = c("Total richness (S')", "Core richness")
    )
  )

panel_B_data <- summary_year %>%
  dplyr::select(year, trans_mean, trans_se, prop_mean, prop_se) %>%
  tidyr::pivot_longer(
    cols = -year,
    names_to = c("metric", ".value"),
    names_pattern = "(trans|prop)_(mean|se)"
  ) %>%
  dplyr::mutate(
    metric = factor(
      metric,
      levels = c("trans", "prop"),
      labels = c("Transient richness", "Transient proportion")
    )
  )

panel_C_data <- summary_trt %>%
  dplyr::select(trt, total_mean, total_se, core_mean, core_se) %>%
  tidyr::pivot_longer(
    cols = -trt,
    names_to = c("metric", ".value"),
    names_pattern = "(total|core)_(mean|se)"
  ) %>%
  dplyr::mutate(
    metric = factor(
      metric,
      levels = c("total", "core"),
      labels = c("Total richness (S')", "Core richness")
    )
  )

panel_D_data <- summary_trt %>%
  dplyr::select(trt, trans_mean, trans_se, prop_mean, prop_se) %>%
  tidyr::pivot_longer(
    cols = -trt,
    names_to = c("metric", ".value"),
    names_pattern = "(trans|prop)_(mean|se)"
  ) %>%
  dplyr::mutate(
    metric = factor(
      metric,
      levels = c("trans", "prop"),
      labels = c("Transient richness", "Transient proportion")
    )
  )


############################################################
## 5. Colours
############################################################

metric_cols <- c(
  "Total richness (S')" = "black",
  "Core richness" = "#00A83A",
  "Transient richness" = "red",
  "Transient proportion" = "#F4A300"
)


############################################################
## 6. Common theme
############################################################

theme_fig3 <- theme_classic(base_family = "Times New Roman") +
  theme(
    axis.text.x = element_text(size = 13, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 16, color = "black", face = "bold"),
    axis.title.y = element_text(size = 16, color = "black", face = "bold"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA)
  )


############################################################
## 7. Panel A
############################################################


panel_A_data <- summary_year %>%
  dplyr::select(year, total_mean, total_se, core_mean, core_se) %>%
  tidyr::pivot_longer(
    cols = -year,
    names_to = c("metric", ".value"),
    names_pattern = "(total|core)_(mean|se)"
  ) %>%
  dplyr::mutate(
    metric = dplyr::recode(
      metric,
      total = "Total richness (S')",
      core = "Core richness"
    ),
    metric = factor(
      metric,
      levels = c(
        "Total richness (S')",
        "Core richness",
        "Transient richness",
        "Transient proportion"
      )
    )
  )

panel_A_legend_data <- data.frame(
  year = factor(
    c("2018", "2018", "2018", "2018"),
    levels = c("2018", "2019", "2020")
  ),
  mean = c(NA, NA, NA, NA),
  metric = factor(
    c(
      "Total richness (S')",
      "Core richness",
      "Transient richness",
      "Transient proportion"
    ),
    levels = c(
      "Total richness (S')",
      "Core richness",
      "Transient richness",
      "Transient proportion"
    )
  )
)

panel_A <- ggplot(panel_A_data, aes(x = year, y = mean, color = metric, group = metric)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.08,
    linewidth = 0.7
  ) +
  geom_line(
    data = panel_A_legend_data,
    aes(x = year, y = mean, color = metric, group = metric),
    linewidth = 1.1,
    show.legend = TRUE,
    na.rm = TRUE
  ) +
  geom_point(
    data = panel_A_legend_data,
    aes(x = year, y = mean, color = metric),
    size = 3,
    show.legend = TRUE,
    na.rm = TRUE
  ) +
  annotate(
    "text",
    x = 0.65,
    y = 9.05,
    label = "(A)",
    family = "Times New Roman",
    fontface = "bold",
    size = 5
  ) +
  scale_color_manual(
    values = metric_cols,
    breaks = c(
      "Total richness (S')",
      "Core richness",
      "Transient richness",
      "Transient proportion"
    ),
    drop = FALSE
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        color = c("black", "#00A83A", "red", "#F4A300"),
        linewidth = c(1.1, 1.1, 1.1, 1.1),
        size = c(3, 3, 3, 3)
      )
    )
  ) +
  scale_y_continuous(
    limits = c(5.85, 9.2),
    breaks = c(6, 7, 8, 9),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_fig3 +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(
      size = 16,
      color = "black",
      face = "bold",
      vjust = 0.5,
      margin = margin(r = 10)
    ),
    legend.title = element_blank(),
    legend.text = element_text(size = 12, color = "black"),
    legend.position = c(0.66, 0.28),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA)
  )

panel_A

############################################################
## 8. Panel B
############################################################

panel_B <- ggplot(panel_B_data, aes(x = year, y = mean, color = metric, group = metric)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.08,
    linewidth = 0.7
  ) +
  annotate(
    "text",
    x = 0.65,
    y = 0.97,
    label = "(B)",
    family = "Times New Roman",
    fontface = "bold",
    size = 5
  ) +
  scale_color_manual(values = metric_cols, drop = FALSE) +
  scale_y_continuous(
    limits = c(-0.05, 1.05),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    expand = c(0, 0)
  ) +
  labs(
    x = "Years",
    y = NULL
  ) +
  theme_fig3 +
  theme(
    legend.position = "none"
  )


############################################################
## 9. Panel C
############################################################

panel_C <- ggplot(panel_C_data, aes(x = trt, y = mean, color = metric, group = metric)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.08,
    linewidth = 0.7
  ) +
  annotate(
    "text",
    x = 0.78,
    y = 9.05,
    label = "(C)",
    family = "Times New Roman",
    fontface = "bold",
    size = 5
  ) +
  scale_color_manual(values = metric_cols, drop = FALSE) +
  scale_y_continuous(
    limits = c(5.85, 9.2),
    breaks = c(6, 7, 8, 9),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_fig3 +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none"
  )


############################################################
## 10. Panel D
############################################################

panel_D <- ggplot(panel_D_data, aes(x = trt, y = mean, color = metric, group = metric)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.08,
    linewidth = 0.7
  ) +
  annotate(
    "text",
    x = 0.78,
    y = 0.97,
    label = "(D)",
    family = "Times New Roman",
    fontface = "bold",
    size = 5
  ) +
  scale_color_manual(values = metric_cols, drop = FALSE) +
  scale_y_continuous(
    limits = c(-0.05, 1.05),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    expand = c(0, 0)
  ) +
  labs(
    x = "Thinning",
    y = NULL
  ) +
  theme_fig3 +
  theme(
    legend.position = "none"
  )


############################################################
## 11. Combine Figure 3 correctly
############################################################

panel_A <- panel_A + labs(y = NULL)
panel_B <- panel_B + labs(y = NULL)
panel_C <- panel_C + labs(y = NULL)
panel_D <- panel_D + labs(y = NULL)

shared_y <- patchwork::wrap_elements(
  full = grid::textGrob(
    "Richness",
    rot = 90,
    gp = grid::gpar(
      fontfamily = "Times New Roman",
      fontsize = 16,
      fontface = "bold"
    )
  )
)

design_fig3 <- "
LAC
LBD
"

fig3 <- patchwork::wrap_plots(
  L = shared_y,
  A = panel_A,
  C = panel_C,
  B = panel_B,
  D = panel_D,
  design = design_fig3,
  widths = c(0.04, 1, 1),
  heights = c(1, 1)
) &
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

fig3

############################################################
## 12. Save Figure 3
############################################################

ggsave(
  filename = "Figure3_core_transient_richness.tiff",
  plot = fig3,
  dpi = 600,
  width = 10,
  height = 7.5,
  units = "in",
  compression = "lzw",
  bg = "white"
)



################################# Biomass #################################################

############################################################
## Biomass and richness analysis
############################################################

str(Biomass_t)
str(richness_df)




############################################################
## 1. Strip-level richness, 2019–2020 only
############################################################

richness_strip <- richness_df %>%
  mutate(
    year  = as.numeric(as.character(year)),
    block = as.character(block),
    label = as.character(label),
    trt   = as.character(trt)
  ) %>%
  filter(year %in% c(2019, 2020)) %>%
  group_by(block, label, year, trt, distance, elevation, slope, orientation) %>%
  summarise(
    total_richness     = mean(total_richness, na.rm = TRUE),
    core_richness      = mean(core_richness, na.rm = TRUE),
    transient_richness = mean(transient_richness, na.rm = TRUE),
    transient_prop     = mean(transient_prop, na.rm = TRUE),
    .groups = "drop"
  )


############################################################
## 2. Clean biomass table
############################################################

Biomass_t2 <- Biomass_t %>%
  rename(
    block     = blcok,
    AGB_total = AGB
  ) %>%
  mutate(
    block = as.character(block),
    label = as.character(label),
    year  = as.numeric(year)
  ) %>%
  dplyr::select(block, label, year, AGB_total)


############################################################
## 3. Join biomass with strip-level richness
############################################################

biomass_richness <- Biomass_t2 %>%
  inner_join(richness_strip, by = c("block", "label", "year")) %>%
  mutate(
    block = factor(block),
    year  = factor(year),
    trt   = factor(trt)
  )

biomass_richness


############################################################
## 4. Save joined table
############################################################

write_xlsx(
  biomass_richness,
  path = "C:/Users/Research Greenhouse/Desktop/paper goudie/goudie2/data/biomass_richness.xlsx"
)


############################################################
## 5. Log-transformed richness variables
############################################################

biomass_richness <- biomass_richness %>%
  mutate(
    log_transient = log1p(transient_richness),
    log_core      = log1p(core_richness)
  )


############################################################
## 6. Mixed models
############################################################

## Transient richness model
model_log_transient <- lmer(
  AGB_total ~ log_transient * trt * year + (1 | block),
  data = biomass_richness
)

summary(model_log_transient)
Anova(model_log_transient, type = 2)


## Core richness model
model_log_core <- lmer(
  AGB_total ~ log_core * trt * year + (1 | block),
  data = biomass_richness
)

summary(model_log_core)
Anova(model_log_core, type = 2)


## Total biomass model
total_model <- lmer(
  AGB_total ~ trt * year + (1 | block),
  data = biomass_richness
)

summary(total_model)
Anova(total_model, type = 2)


############################################################
## Biomass figure cleaning + patchwork


############################################################
## 1. Clean factors and log variable
############################################################

biomass_richness <- biomass_richness %>%
  mutate(
    year = factor(year, levels = c("2019", "2020")),
    trt  = factor(trt, levels = c("10 m", "15 m", "20 m", "Control")),
    log_transient = log10(transient_richness + 1)
  )

############################################################
## 2. Summary tables for bar plots
############################################################

year_means <- biomass_richness %>%
  group_by(year) %>%
  summarise(
    mean_biomass = mean(AGB_total, na.rm = TRUE),
    se_biomass   = sd(AGB_total, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

trt_means <- biomass_richness %>%
  group_by(trt) %>%
  summarise(
    mean_biomass = mean(AGB_total, na.rm = TRUE),
    se_biomass   = sd(AGB_total, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

############################################################
## 3. Tukey letters from total biomass model
############################################################



## Year letters
year_letters <- multcomp::cld(
  emmeans(total_model, ~ year),
  Letters = letters,
  adjust = "tukey"
) %>%
  as.data.frame() %>%
  mutate(
    year = factor(as.character(year), levels = c("2019", "2020")),
    .group = gsub(" ", "", .group)
  ) %>%
  dplyr::select(year, .group)

## Treatment letters
trt_letters <- multcomp::cld(
  emmeans(total_model, ~ trt),
  Letters = letters,
  adjust = "tukey"
) %>%
  as.data.frame() %>%
  mutate(
    trt = factor(as.character(trt), levels = c("10 m", "15 m", "20 m", "Control")),
    .group = gsub(" ", "", .group)
  ) %>%
  dplyr::select(trt, .group)

## Join letters to summary tables
year_means <- year_means %>%
  left_join(year_letters, by = "year") %>%
  mutate(label_y = mean_biomass + se_biomass + 4)

trt_means <- trt_means %>%
  left_join(trt_letters, by = "trt") %>%
  mutate(label_y = mean_biomass + se_biomass + 4)

############################################################
## 4. Treatment colours
############################################################

trt_colors <- c(
  "10 m"    = "#6A6A6A",
  "15 m"    = "#B8860B",
  "20 m"    = "#377EB8",
  "Control" = "black"
)

############################################################
## 5. Common theme
############################################################

theme_biomass <- theme_classic(
  base_size = 14,
  base_family = "Times New Roman"
) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.7),
    axis.text.x = element_text(
      color = "black",
      size = 16,
      family = "Times New Roman"
    ),
    axis.text.y = element_text(
      color = "black",
      size = 16,
      family = "Times New Roman"
    ),
    axis.title.x = element_text(
      color = "black",
      size = 18,
      family = "Times New Roman",
      face = "bold"
    ),
    axis.title.y = element_text(
      color = "black",
      size = 18,
      family = "Times New Roman",
      face = "bold"
    ),
    legend.position = "none",
    plot.margin = margin(8, 8, 8, 8)
  )

############################################################
## 6. Panel A: biomass by year
############################################################

pA <- ggplot(year_means, aes(x = year, y = mean_biomass)) +
  geom_col(
    fill = "grey30",
    width = 0.7,
    color = "black"
  ) +
  geom_errorbar(
    aes(
      ymin = mean_biomass - se_biomass,
      ymax = mean_biomass + se_biomass
    ),
    width = 0.15,
    color = "black",
    linewidth = 0.6
  ) +
  geom_text(
    aes(y = label_y, label = .group),
    size = 6,
    family = "Times New Roman"
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = "(A)",
    hjust = -0.45,
    vjust = 1.35,
    size = 6,
    fontface = "bold",
    family = "Times New Roman"
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 60),
    breaks = seq(0, 60, by = 10)
  ) +
  labs(
    x = "Years",
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_biomass

############################################################
## 7. Panel B: biomass by treatment
############################################################

pB <- ggplot(trt_means, aes(x = trt, y = mean_biomass, fill = trt)) +
  geom_col(
    width = 0.7,
    color = "black"
  ) +
  geom_errorbar(
    aes(
      ymin = mean_biomass - se_biomass,
      ymax = mean_biomass + se_biomass
    ),
    width = 0.15,
    color = "black",
    linewidth = 0.6
  ) +
  geom_text(
    aes(y = label_y, label = .group),
    size = 6,
    family = "Times New Roman"
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = "(B)",
    hjust = -0.45,
    vjust = 1.35,
    size = 6,
    fontface = "bold",
    family = "Times New Roman"
  ) +
  scale_fill_manual(values = trt_colors) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 65),
    breaks = seq(0, 60, by = 10)
  ) +
  labs(
    x = "Strips-thinning",
    y = "Biomass (g)"
  ) +
  coord_cartesian(clip = "off") +
  theme_biomass

############################################################
## 8. Panel C: biomass vs transient richness
############################################################

biomass_plot <- biomass_richness %>%
  filter(trt %in% c("10 m", "15 m", "20 m")) %>%
  droplevels()

############################################################
## 9. R2 values per treatment
############################################################

m10 <- lm(AGB_total ~ log_transient, data = subset(biomass_plot, trt == "10 m"))
m15 <- lm(AGB_total ~ log_transient, data = subset(biomass_plot, trt == "15 m"))
m20 <- lm(AGB_total ~ log_transient, data = subset(biomass_plot, trt == "20 m"))

r2_10 <- summary(m10)$r.squared
r2_15 <- summary(m15)$r.squared
r2_20 <- summary(m20)$r.squared

r2_labels <- data.frame(
  trt = factor(c("10 m", "15 m", "20 m"), levels = c("10 m", "15 m", "20 m")),
  x = c(0.32, 0.35, 0.52),
  y = c(8, 30, 46),
  label = c(
    paste0("R² = ", round(r2_10, 2)),
    paste0("R² = ", round(r2_15, 2)),
    paste0("R² = ", round(r2_20, 2))
  )
)

pC <- ggplot(
  biomass_plot,
  aes(x = log_transient, y = AGB_total, colour = trt, fill = trt)
) +
  geom_point(
    size = 2,
    alpha = 0.8
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 1.1,
    alpha = 0.18
  ) +
  geom_text(
    data = r2_labels,
    aes(x = x, y = y, label = label, colour = trt),
    inherit.aes = FALSE,
    size = 5,
    family = "Times New Roman",
    fontface = "bold",
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = "(C)",
    hjust = -0.45,
    vjust = 1.35,
    size = 6,
    fontface = "bold",
    family = "Times New Roman"
  ) +
  scale_colour_manual(values = trt_colors) +
  scale_fill_manual(values = trt_colors) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 100),
    breaks = seq(0, 100, by = 25)
  ) +
  labs(
    x = "Log 10 Transient richness",
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_biomass

############################################################
## 10. Patchwork: vertical structure like manuscript
############################################################

fig5 <- pA / pB / pC +
  plot_layout(heights = c(1, 1, 1.25))

fig5

############################################################
## 11. Save final figure
############################################################

ggsave(
  filename = "C:/GitHub/TRU_Strip_thinning_community/Figure_4_AGB_biomass_final.tiff",
  plot = fig5,
  dpi = 600,
  width = 7,
  height = 10,
  units = "in",
  compression = "lzw"
)





############################################################
##  C:N ratio 




############################################################
## 1. DATA PREPARATION
############################################################

str(richness_cn_df)
head(richness_cn_df)

richness_cn_df <- richness_cn_df %>%
  mutate(
    block = factor(block),
    label = factor(label),
    year  = factor(year, levels = c("2018", "2019", "2020")),
    trt   = factor(trt, levels = c("Control", "10 m", "15 m", "20 m")),
    CNH = TCH / TNH,        # C:N ratio, 0–10 cm
    CNL = TCL / TNL,        # C:N ratio, 10–20 cm
    log_transient = log1p(transient_richness),
    log_core = log(core_richness)
  )

summary(richness_cn_df$CNH)
summary(richness_cn_df$CNL)

cor(
  richness_cn_df$core_richness,
  richness_cn_df$transient_richness
)


############################################################
## 2. Stats: carbon and nitrogen models
############################################################

## Carbon, 0–10 cm
m_TCH <- lmer(
  TCH ~ core_richness + transient_richness + trt + year +
    (1 | block) + (1 | label),
  data = richness_cn_df
)

summary(m_TCH)
anova(m_TCH, type = 3)


## Nitrogen, 0–10 cm
m_TNH <- lmer(
  TNH ~ core_richness + transient_richness * trt * year +
    (1 | block) + (1 | label),
  data = richness_cn_df
)

summary(m_TNH)
anova(m_TNH, type = 3)


## Carbon, 10–20 cm
m_TCL <- lmer(
  TCL ~ core_richness + transient_richness + trt + year +
    (1 | block) + (1 | label),
  data = richness_cn_df
)

summary(m_TCL)
anova(m_TCL, type = 3)


## Nitrogen, 10–20 cm
m_TNL <- lmer(
  TNL ~ core_richness + transient_richness + trt + year +
    (1 | block) + (1 | label),
  data = richness_cn_df
)

summary(m_TNL)
anova(m_TNL, type = 3)


############################################################
## 3. Stats: C:N ratio models
############################################################

## C:N ratio, 0–10 cm
m_CNH <- lmer(
  CNH ~ core_richness + transient_richness + trt + year +
    (1 | block / label),
  data = richness_cn_df
)

summary(m_CNH)
anova(m_CNH, type = 3)


## C:N ratio, 10–20 cm
m_CNL <- lmer(
  CNL ~ core_richness + transient_richness + trt + year +
    (1 | block / label),
  data = richness_cn_df
)

summary(m_CNL)
anova(m_CNL, type = 3)


############################################################
## 4. Stats for Panel A: CNH by year
############################################################

m_CNH_year <- lmer(
  CNH ~ year + (1 | block / label),
  data = richness_cn_df
)

summary(m_CNH_year)
Anova(m_CNH_year, type = 2)

CNH_year_emm <- emmeans(m_CNH_year, ~ year)

CNH_year_letters <- cld(
  CNH_year_emm,
  Letters = letters,
  adjust = "tukey"
)

CNH_year_letters


############################################################
## 5. Stats for Panel C: CNL by year
############################################################

m_CNL_year <- lmer(
  CNL ~ year + (1 | block / label),
  data = richness_cn_df
)

summary(m_CNL_year)
Anova(m_CNL_year, type = 2)

CNL_year_emm <- emmeans(m_CNL_year, ~ year)

CNL_year_letters <- cld(
  CNL_year_emm,
  Letters = letters,
  adjust = "tukey"
)

CNL_year_letters


############################################################
## 6. Stats for Panel B: CNH vs transient richness
############################################################

m_CNH_trans_fig <- lm(
  CNH ~ log_transient,
  data = richness_cn_df
)

summary(m_CNH_trans_fig)

CNH_trans_R2 <- summary(m_CNH_trans_fig)$r.squared
CNH_trans_p  <- summary(m_CNH_trans_fig)$coefficients["log_transient", "Pr(>|t|)"]

CNH_trans_label <- paste0(
  "R² = ", round(CNH_trans_R2, 2), " ",
  ifelse(CNH_trans_p < 0.001, "***",
         ifelse(CNH_trans_p < 0.01, "**",
                ifelse(CNH_trans_p < 0.05, "*", "ns")))
)

CNH_trans_stats <- data.frame(
  response = "CNH",
  predictor = "log1p(transient_richness)",
  R2 = CNH_trans_R2,
  p_value = CNH_trans_p,
  label = CNH_trans_label
)

CNH_trans_stats


############################################################
## 7. Stats for Panel D: CNL vs transient richness by treatment
############################################################

CNL_trans_stats_by_trt <- richness_cn_df %>%
  filter(trt %in% c("10 m", "15 m", "20 m")) %>%
  group_by(trt) %>%
  do({
    m <- lm(CNL ~ log_transient, data = .)
    data.frame(
      slope = coef(m)[["log_transient"]],
      R2 = summary(m)$r.squared,
      p_value = summary(m)$coefficients["log_transient", "Pr(>|t|)"]
    )
  }) %>%
  ungroup() %>%
  mutate(
    sig = ifelse(
      p_value < 0.001, "***",
      ifelse(p_value < 0.01, "**",
             ifelse(p_value < 0.05, "*", "ns"))
    ),
    label = paste0("R² = ", round(R2, 3), " ", sig)
  )

CNL_trans_stats_by_trt


############################################################
## 8. Plot data
############################################################

richness_cn_plot <- richness_cn_df

richness_cn_D <- richness_cn_plot %>%
  filter(trt %in% c("10 m", "15 m", "20 m"))


############################################################
## 9. Panel A: CNH by year
############################################################

panel_A <- ggplot(richness_cn_plot, aes(x = year, y = CNH)) +
  geom_violin(
    fill = "#CFE8FA",
    color = "black",
    linewidth = 0.7,
    trim = FALSE
  ) +
  geom_boxplot(
    width = 0.18,
    fill = "white",
    color = "black",
    linewidth = 0.6,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.08,
    size = 1.4,
    alpha = 0.45,
    color = "grey30"
  ) +
  annotate("text", x = 1, y = 42, label = "a",
           family = "Times New Roman", fontface = "bold", size = 5) +
  annotate("text", x = 2, y = 42, label = "b",
           family = "Times New Roman", fontface = "bold", size = 5) +
  annotate("text", x = 3, y = 42, label = "c",
           family = "Times New Roman", fontface = "bold", size = 5) +
  annotate("text", x = 0.45, y = 42.8, label = "(A)",
           family = "Times New Roman", fontface = "bold", size = 5) +
  scale_y_continuous(
    limits = c(0, 44),
    breaks = seq(0, 40, by = 10),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 14, base_family = "Times New Roman") +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_blank(),
    legend.position = "none"
  )


############################################################
## 10. Panel C: CNL by year
############################################################

panel_C <- ggplot(richness_cn_plot, aes(x = year, y = CNL)) +
  geom_violin(
    fill = "#D8B7F0",
    color = "black",
    linewidth = 0.7,
    trim = FALSE
  ) +
  geom_boxplot(
    width = 0.18,
    fill = "white",
    color = "black",
    linewidth = 0.6,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.08,
    size = 1.4,
    alpha = 0.45,
    color = "grey30"
  ) +
  annotate("text", x = 1, y = 43, label = "a",
           family = "Times New Roman", fontface = "bold", size = 5) +
  annotate("text", x = 2, y = 31, label = "b",
           family = "Times New Roman", fontface = "bold", size = 5) +
  annotate("text", x = 3, y = 31, label = "b",
           family = "Times New Roman", fontface = "bold", size = 5) +
  annotate("text", x = 0.45, y = 44, label = "(C)",
           family = "Times New Roman", fontface = "bold", size = 5) +
  scale_y_continuous(
    limits = c(0, 45),
    breaks = seq(0, 40, by = 10),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 14, base_family = "Times New Roman") +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_blank(),
    legend.position = "none"
  )


############################################################
## 11. Panel B: CNH vs transient richness
############################################################

panel_B <- ggplot(richness_cn_plot, aes(x = log_transient, y = CNH)) +
  geom_point(
    size = 1.6,
    alpha = 0.55,
    color = "grey30"
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black",
    fill = "#CFE8FA",
    linewidth = 0.8
  ) +
  annotate("text", x = 1.1, y = 28, label = CNH_trans_label,
           family = "Times New Roman", fontface = "bold", size = 5) +
  annotate("text", x = -0.08, y = 42.8, label = "(B)",
           family = "Times New Roman", fontface = "bold", size = 5) +
  scale_y_continuous(
    limits = c(8, 44),
    breaks = seq(10, 40, by = 10),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 14, base_family = "Times New Roman") +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_blank(),
    legend.position = "none"
  )


############################################################
## 12. Panel D: CNL vs transient richness by treatment
############################################################

panel_D <- ggplot(richness_cn_D, aes(x = log_transient, y = CNL)) +
  geom_point(
    size = 1.6,
    alpha = 0.55,
    color = "grey30"
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black",
    fill = "#D8B7F0",
    linewidth = 0.8
  ) +
  geom_text(
    data = CNL_trans_stats_by_trt,
    aes(x = 0.35, y = 38, label = label),
    inherit.aes = FALSE,
    family = "Times New Roman",
    fontface = "bold",
    size = 4.5
  ) +
  annotate("text", x = -0.08, y = 42.8, label = "(D)",
           family = "Times New Roman", fontface = "bold", size = 5) +
  facet_wrap(~ trt, nrow = 1) +
  scale_y_continuous(
    limits = c(6, 44),
    breaks = seq(10, 40, by = 10),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 14, base_family = "Times New Roman") +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.text = element_text(color = "black", size = 13),
    axis.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(0.12, "lines"),
    legend.position = "none"
  )


############################################################
## 13. Combine final figure
############################################################

fig_CN <- patchwork::wrap_plots(
  panel_A, panel_C,
  panel_B, panel_D,
  ncol = 2
)

fig_CN


############################################################
## 14. Save final figure
############################################################

ggsave(
  filename = "Figure_CN_ratio_final.tiff",
  plot = fig_CN,
  path = "C:/Users/Research Greenhouse/Desktop/paper goudie/goudie2/data",
  dpi = 600,
  width = 10,
  height = 8,
  units = "in",
  compression = "lzw"
)