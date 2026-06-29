
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


## ---------------------------------------------------------------
## BLOCK A1 — Prepare metadata and community matrices
## ---------------------------------------------------------------

env_cols <- c("block", "year", "trt", "distance", "elevation", "slope", "orientation")

metadata <- final_table %>%
  select(all_of(env_cols)) %>%
  mutate(
    block       = factor(block),
    year        = factor(year, levels = c("2018", "2019", "2020")),
    trt         = factor(trt),
    orientation = factor(orientation)
  )

str(metadata)

species_cols <- setdiff(names(final_table), env_cols)

comm_cover <- final_table %>%
  select(all_of(species_cols))

comm_pa <- comm_cover %>%
  mutate(across(everything(), ~ ifelse(.x > 0, 1, 0)))

comm_pa <- comm_pa[, colSums(comm_pa) > 0]

dim(comm_pa)
head(comm_pa[, 1:10])


## ---------------------------------------------------------------
## BLOCK A2 — Jaccard partitioning with betapart
## ---------------------------------------------------------------

beta_jac <- betapart::beta.pair(comm_pa, index.family = "jaccard")

jac_total    <- beta_jac$beta.jac
jac_turnover <- beta_jac$beta.jtu
jac_nested   <- beta_jac$beta.jne

mean(jac_turnover)
mean(jac_nested)
mean(jac_total)


## ---------------------------------------------------------------
## BLOCK A3 — Figure 2A: Turnover vs nestedness
## ---------------------------------------------------------------

turn_vals <- as.numeric(jac_turnover)
nest_vals <- as.numeric(jac_nested)

beta_box <- tibble(
  component = c(
    rep("Nestedness", length(nest_vals)),
    rep("Turnover", length(turn_vals))
  ),
  value = c(nest_vals, turn_vals)
) %>%
  drop_na()

light_grey <- "#d9d9d9"
dark_grey  <- "#4d4d4d"

fig2A <- ggplot(beta_box, aes(x = component, y = value, fill = component)) +
  geom_boxplot(outlier.size = 0.6, color = "black", width = 0.6) +
  scale_fill_manual(values = c(
    "Nestedness" = light_grey,
    "Turnover"   = dark_grey
  )) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(
    axis.text.x = element_text(family = "Times New Roman", size = 16, color = "black"),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color = "black"),
    legend.position = "none"
  )

fig2A


## ---------------------------------------------------------------
## BLOCK A4 — PERMANOVA on Jaccard turnover
## ---------------------------------------------------------------

set.seed(123)

permanova_turnover <- vegan::adonis2(
  jac_turnover ~ year + trt + distance + elevation + slope + orientation,
  data   = metadata,
  strata = metadata$block
)

permanova_turnover


## ---------------------------------------------------------------
## BLOCK A5 — Pairwise PERMANOVA: turnover by year
## ---------------------------------------------------------------

run_pairwise_year_permanova <- function(year1, year2, metadata, comm_pa) {
  
  keep_rows <- metadata$year %in% c(year1, year2)
  
  meta_sub <- metadata[keep_rows, , drop = FALSE]
  comm_sub <- comm_pa[keep_rows, , drop = FALSE]
  
  beta_sub <- betapart::beta.pair(comm_sub, index.family = "jaccard")
  jac_sub  <- beta_sub$beta.jtu
  
  res <- vegan::adonis2(
    jac_sub ~ year,
    data   = meta_sub,
    strata = meta_sub$block
  )
  
  list(
    years = paste(year1, year2, sep = " vs "),
    permanova = res
  )
}

pair_2018_2019 <- run_pairwise_year_permanova("2018", "2019", metadata, comm_pa)
pair_2019_2020 <- run_pairwise_year_permanova("2019", "2020", metadata, comm_pa)
pair_2018_2020 <- run_pairwise_year_permanova("2018", "2020", metadata, comm_pa)

pair_2018_2019$permanova
pair_2019_2020$permanova
pair_2018_2020$permanova


## ---------------------------------------------------------------
## BLOCK A6 — Bootstrap SE of pairwise R2 values
## ---------------------------------------------------------------

boot_R2_pair <- function(year1, year2, metadata, comm_pa, n_boot = 199) {
  
  keep <- metadata$year %in% c(year1, year2)
  
  meta0 <- droplevels(metadata[keep, ])
  comm0 <- comm_pa[keep, , drop = FALSE]
  
  idx1 <- which(meta0$year == year1)
  idx2 <- which(meta0$year == year2)
  
  r2_boot <- numeric(n_boot)
  
  for (b in seq_len(n_boot)) {
    
    idx <- c(
      sample(idx1, length(idx1), replace = TRUE),
      sample(idx2, length(idx2), replace = TRUE)
    )
    
    meta_b <- meta0[idx, , drop = FALSE]
    comm_b <- comm0[idx, , drop = FALSE]
    
    beta_b <- betapart::beta.pair(comm_b, index.family = "jaccard")
    jac_b  <- beta_b$beta.jtu
    
    ad_b <- vegan::adonis2(
      jac_b ~ year,
      data         = meta_b,
      strata       = meta_b$block,
      permutations = 0
    )
    
    r2_boot[b] <- ad_b$R2[1]
  }
  
  r2_boot <- r2_boot[!is.na(r2_boot)]
  sd(r2_boot) * 100
}


## ---------------------------------------------------------------
## BLOCK A7 — Figure 2B: Pairwise PERMANOVA R2
## ---------------------------------------------------------------

set.seed(123)

pairwise_R2 <- tibble(
  comparison = factor(
    c("2018 vs 2019", "2018 vs 2020", "2019 vs 2020"),
    levels = c("2018 vs 2019", "2018 vs 2020", "2019 vs 2020")
  ),
  R2_percent = c(
    pair_2018_2019$permanova$R2[1],
    pair_2018_2020$permanova$R2[1],
    pair_2019_2020$permanova$R2[1]
  ) * 100,
  SE_percent = c(
    boot_R2_pair("2018", "2019", metadata, comm_pa),
    boot_R2_pair("2018", "2020", metadata, comm_pa),
    boot_R2_pair("2019", "2020", metadata, comm_pa)
  )
) %>%
  mutate(fill_col = c("#4d4d4d", "#8c8c8c", "#d9d9d9"))

pairwise_R2

fig2B <- ggplot(pairwise_R2, aes(x = comparison, y = R2_percent, fill = fill_col)) +
  geom_col(color = "black", width = 0.6) +
  geom_errorbar(
    aes(
      ymin = R2_percent - SE_percent,
      ymax = R2_percent + SE_percent
    ),
    width = 0.15,
    linewidth = 0.7
  ) +
  scale_fill_identity() +
  scale_y_continuous(limits = c(0, 20), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(
    axis.text.x = element_text(family = "Times New Roman", size = 16, color = "black"),
    axis.text.y = element_text(family = "Times New Roman", size = 16, color = "black"),
    legend.position = "none"
  )

fig2B


## ---------------------------------------------------------------
## BLOCK A8 — Save figures
## ---------------------------------------------------------------

ggsave(
  filename = "Fig2A_boxplot.png",
  plot = fig2A,
  path = "C:/Users/Research Greenhouse/Desktop/paper goudie/goudie2/data",
  width = 6,
  height = 5,
  dpi = 600
)

ggsave(
  filename = "Fig2B_clean.png",
  plot = fig2B,
  path = "C:/Users/Research Greenhouse/Desktop/paper goudie/goudie2/data",
  width = 7,
  height = 5,
  dpi = 600
)


################################ Core vs transient richness ################################

library(tidyverse)

## 1. Define metadata and species columns

env_cols <- c(
  "block", "year", "trt", "distance",
  "elevation", "slope", "orientation",
  "label", "sampling"
)

env_cols <- intersect(env_cols, names(final_div_table1))

species_cols <- setdiff(names(final_div_table1), env_cols)


## 2. Convert species cover to presence-absence

comm_pa <- final_div_table1 %>%
  select(all_of(species_cols)) %>%
  mutate(across(everything(), ~ ifelse(.x > 0, 1, 0)))

comm_pa <- comm_pa[, colSums(comm_pa, na.rm = TRUE) > 0, drop = FALSE]

species_pa <- colnames(comm_pa)


## 3. Calculate species temporal occupancy across years

df_pa <- bind_cols(
  final_div_table1 %>% select(any_of(c("year", "label"))),
  comm_pa
)

species_occupancy <- df_pa %>%
  group_by(year) %>%
  summarise(
    across(all_of(species_pa), ~ max(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  summarise(
    across(all_of(species_pa), ~ sum(.x, na.rm = TRUE))
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "species",
    values_to = "n_years"
  )


## 4. Classify species: core vs transient only

species_class <- species_occupancy %>%
  mutate(
    class = case_when(
      n_years >= 2 ~ "core",
      n_years == 1 ~ "transient",
      TRUE ~ NA_character_
    )
  )

table(species_class$class)


## 5. Define core and transient species lists

core_species <- species_class %>%
  filter(class == "core") %>%
  pull(species)

transient_species <- species_class %>%
  filter(class == "transient") %>%
  pull(species)


## 6. Create richness_df

richness_df <- final_div_table1 %>%
  mutate(
    total_richness = rowSums(comm_pa, na.rm = TRUE),
    core_richness = rowSums(comm_pa[, core_species, drop = FALSE], na.rm = TRUE),
    transient_richness = rowSums(comm_pa[, transient_species, drop = FALSE], na.rm = TRUE),
    transient_prop = ifelse(
      total_richness > 0,
      transient_richness / total_richness,
      NA
    )
  ) %>%
  select(
    any_of(c(
      "block", "label", "sampling", "year", "trt",
      "distance", "elevation", "slope", "orientation"
    )),
    total_richness,
    core_richness,
    transient_richness,
    transient_prop
  )


## 7. Check outputs

head(species_class)
head(richness_df)

summary(richness_df$total_richness)
summary(richness_df$core_richness)
summary(richness_df$transient_richness)
summary(richness_df$transient_prop)

################################ stats Core vs transient  ################################


## ---------------------------------------------------------------
## 1. Prepare data
## ---------------------------------------------------------------

richness_df2 <- richness_df %>%
  mutate(
    block = factor(block),
    year = factor(year, levels = c("2018", "2019", "2020")),
    trt = factor(trt, levels = c("Control", "10 m", "15 m", "20 m")),
    orientation = factor(orientation)
  )


## ---------------------------------------------------------------
## 2. Fit final models
## ---------------------------------------------------------------

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


## ---------------------------------------------------------------
## 3. Type II Wald tests
## ---------------------------------------------------------------

anova_total <- car::Anova(mod_total, type = 2)
anova_core <- car::Anova(mod_core, type = 2)
anova_transient <- car::Anova(mod_transient, type = 2)
anova_transient_prop <- car::Anova(mod_transient_prop, type = 2)

anova_total
anova_core
anova_transient
anova_transient_prop




## ---------------------------------------------------------------
## 4. DHARMa diagnostic table
## ---------------------------------------------------------------

get_dharma_tests <- function(model, model_name) {
  
  sim <- simulateResiduals(model)
  
  disp <- testDispersion(sim)
  zi <- testZeroInflation(sim)
  
  tibble(
    model = model_name,
    dispersion_statistic = as.numeric(disp$statistic),
    dispersion_p = disp$p.value,
    zero_inflation_statistic = as.numeric(zi$statistic),
    zero_inflation_p = zi$p.value
  )
}

dharma_table <- bind_rows(
  get_dharma_tests(mod_total, "Total richness"),
  get_dharma_tests(mod_core, "Core richness"),
  get_dharma_tests(mod_transient, "Transient richness"),
  get_dharma_tests(mod_transient_prop, "Transient proportion")
)

dharma_table



################################ Figure 3 plotting ################################




##  data preparation


richness_df2 <- richness_df %>%
  mutate(
    block = factor(block),
    year = factor(year, levels = c("2018", "2019", "2020")),
    trt = factor(trt, levels = c("10 m", "15 m", "20 m", "Control")),
    orientation = factor(orientation)
  )

summary_year <- richness_df2 %>%
  group_by(year) %>%
  summarise(
    total_mean = mean(total_richness, na.rm = TRUE),
    total_se   = sd(total_richness, na.rm = TRUE) / sqrt(n()),
    core_mean  = mean(core_richness, na.rm = TRUE),
    core_se    = sd(core_richness, na.rm = TRUE) / sqrt(n()),
    trans_mean = mean(transient_richness, na.rm = TRUE),
    trans_se   = sd(transient_richness, na.rm = TRUE) / sqrt(n()),
    prop_mean  = mean(transient_prop, na.rm = TRUE),
    prop_se    = sd(transient_prop, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

summary_trt <- richness_df2 %>%
  group_by(trt) %>%
  summarise(
    total_mean = mean(total_richness, na.rm = TRUE),
    total_se   = sd(total_richness, na.rm = TRUE) / sqrt(n()),
    core_mean  = mean(core_richness, na.rm = TRUE) / sqrt(1) * 0 + mean(core_richness, na.rm = TRUE), # keep style simple
    core_se    = sd(core_richness, na.rm = TRUE) / sqrt(n()),
    trans_mean = mean(transient_richness, na.rm = TRUE),
    trans_se   = sd(transient_richness, na.rm = TRUE) / sqrt(n()),
    prop_mean  = mean(transient_prop, na.rm = TRUE),
    prop_se    = sd(transient_prop, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )


summary_trt <- richness_df2 %>%
  group_by(trt) %>%
  summarise(
    total_mean = mean(total_richness, na.rm = TRUE),
    total_se   = sd(total_richness, na.rm = TRUE) / sqrt(n()),
    core_mean  = mean(core_richness, na.rm = TRUE),
    core_se    = sd(core_richness, na.rm = TRUE) / sqrt(n()),
    trans_mean = mean(transient_richness, na.rm = TRUE),
    trans_se   = sd(transient_richness, na.rm = TRUE) / sqrt(n()),
    prop_mean  = mean(transient_prop, na.rm = TRUE),
    prop_se    = sd(transient_prop, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )


## Panel A data


panel_A_data <- summary_year %>%
  select(year, total_mean, total_se, core_mean, core_se) %>%
  pivot_longer(
    cols = -year,
    names_to = c("metric", ".value"),
    names_pattern = "(total|core)_(mean|se)"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("total", "core", "trans", "prop"),
      labels = c("Total richness (S')", "Core richness",
                 "Transient richness", "Transient proportion")
    )
  )

## add dummy rows so panel A legend shows all 4 series
panel_A_dummy <- tibble(
  year = factor(c("2018", "2018"), levels = c("2018", "2019", "2020")),
  metric = factor(c("Transient richness", "Transient proportion"),
                  levels = c("Total richness (S')", "Core richness",
                             "Transient richness", "Transient proportion")),
  mean = NA_real_,
  se = NA_real_
)

panel_A_data <- bind_rows(panel_A_data, panel_A_dummy)


## Panel B data


panel_B_data <- summary_year %>%
  select(year, trans_mean, trans_se, prop_mean, prop_se) %>%
  pivot_longer(
    cols = -year,
    names_to = c("metric", ".value"),
    names_pattern = "(trans|prop)_(mean|se)"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("trans", "prop"),
      labels = c("Transient richness", "Transient proportion")
    )
  )


## Panel C data


panel_C_data <- summary_trt %>%
  select(trt, total_mean, total_se, core_mean, core_se) %>%
  pivot_longer(
    cols = -trt,
    names_to = c("metric", ".value"),
    names_pattern = "(total|core)_(mean|se)"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("total", "core"),
      labels = c("Total richness (S')", "Core richness")
    )
  )


## Panel D data


panel_D_data <- summary_trt %>%
  select(trt, trans_mean, trans_se, prop_mean, prop_se) %>%
  pivot_longer(
    cols = -trt,
    names_to = c("metric", ".value"),
    names_pattern = "(trans|prop)_(mean|se)"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("trans", "prop"),
      labels = c("Transient richness", "Transient proportion")
    )
  )



metric_cols <- c(
  "Total richness (S')" = "black",
  "Core richness" = "#00A83A",
  "Transient richness" = "red",
  "Transient proportion" = "#F4A300"
)


## Panel A plot


panel_A <- ggplot(panel_A_data, aes(x = year, y = mean, color = metric, group = metric)) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 2.2, na.rm = TRUE) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.08, linewidth = 0.6, na.rm = TRUE) +
  scale_color_manual(
    values = metric_cols,
    breaks = c("Total richness (S')", "Core richness",
               "Transient richness", "Transient proportion"),
    drop = FALSE
  ) +
  labs(x = NULL, y = "Richness", color = NULL) +
  annotate("text", x = 0.88, y = 8.95, label = "(A)",
           family = "Times New Roman", fontface = "bold", size = 5) +
  theme_classic() +
  theme(
    text = element_text(family = "Times New Roman"),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = c(0.68, 0.37),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 12)
  )


## Panel B plot


panel_B <- ggplot(panel_B_data, aes(x = year, y = mean, color = metric, group = metric)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.08, linewidth = 0.6) +
  scale_color_manual(values = metric_cols[c("Transient richness", "Transient proportion")]) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0.02)) +
  labs(x = "Years", y = NULL, color = NULL) +
  annotate("text", x = 0.88, y = 0.98, label = "(B)",
           family = "Times New Roman", fontface = "bold", size = 5) +
  theme_classic() +
  theme(
    text = element_text(family = "Times New Roman"),
    axis.text = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 16, face = "bold"),
    legend.position = "none"
  )


## Panel C plot


panel_C <- ggplot(panel_C_data, aes(x = trt, y = mean, color = metric, group = metric)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.08, linewidth = 0.6) +
  scale_color_manual(values = metric_cols[c("Total richness (S')", "Core richness")]) +
  labs(x = NULL, y = NULL, color = NULL) +
  annotate("text", x = 0.78, y = 9.0, label = "(C)",
           family = "Times New Roman", fontface = "bold", size = 5) +
  theme_classic() +
  theme(
    text = element_text(family = "Times New Roman"),
    axis.text = element_text(size = 14, color = "black"),
    legend.position = "none"
  )


## Panel D plot


panel_D <- ggplot(panel_D_data, aes(x = trt, y = mean, color = metric, group = metric)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.08, linewidth = 0.6) +
  scale_color_manual(values = metric_cols[c("Transient richness", "Transient proportion")]) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0.02)) +
  labs(x = "Thinning", y = NULL, color = NULL) +
  annotate("text", x = 0.78, y = 0.98, label = "(D)",
           family = "Times New Roman", fontface = "bold", size = 5) +
  theme_classic() +
  theme(
    text = element_text(family = "Times New Roman"),
    axis.text = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 16, face = "bold"),
    legend.position = "none"
  )


## Combine


fig3 <- patchwork::wrap_plots(
  panel_A, panel_C,
  panel_B, panel_D,
  ncol = 2
)

fig3

ggsave(
  filename = "Fig3_richness_components_updated.png",
  plot = fig3,
  path = "C:/Users/Research Greenhouse/Desktop/paper goudie/goudie2/data",
  width = 10,
  height = 7,
  dpi = 600
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
## 7. Plot: total biomass by year
############################################################

ggplot(year_means, aes(x = year, y = mean_biomass)) +
  geom_col(fill = "grey30", width = 0.7) +
  geom_errorbar(
    aes(
      ymin = mean_biomass - se_biomass,
      ymax = mean_biomass + se_biomass
    ),
    width = 0.15,
    color = "black"
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 14, base_family = "Times New Roman") +
  theme(
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
    legend.position = "none"
  )


############################################################
## 8. Plot: total biomass by treatment
############################################################

trt_colors <- c(
  "10 m"    = "#6A6A6A",
  "15 m"    = "#B8860B",
  "20 m"    = "#377EB8",
  "Control" = "black"
)

ggplot(trt_means, aes(x = trt, y = mean_biomass, fill = trt)) +
  geom_col(width = 0.7, color = "black") +
  geom_errorbar(
    aes(
      ymin = mean_biomass - se_biomass,
      ymax = mean_biomass + se_biomass
    ),
    width = 0.15,
    color = "black"
  ) +
  scale_fill_manual(values = trt_colors) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 65),
    breaks = seq(0, 60, by = 10)
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 14, base_family = "Times New Roman") +
  theme(
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
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )


############################################################
## 9. Plot: biomass vs transient richness
############################################################

trt_colors <- c(
  "10 m" = "#6A6A6A",
  "15 m" = "#B8860B",
  "20 m" = "#377EB8"
)

y_max <- max(biomass_plot$AGB_total, pred_grid$upper, na.rm = TRUE)

pp <- ggplot() +
  geom_point(
    data = biomass_plot,
    aes(x = log_transient, y = AGB_total, colour = trt),
    size = 2,
    alpha = 0.8
  ) +
  geom_ribbon(
    data = pred_grid,
    aes(x = log_transient, ymin = lower, ymax = upper, fill = trt),
    alpha = 0.18,
    colour = NA
  ) +
  geom_line(
    data = pred_grid,
    aes(x = log_transient, y = pred, colour = trt),
    linewidth = 1.1
  ) +
  scale_colour_manual(values = trt_colors) +
  scale_fill_manual(values = trt_colors) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(ylim = c(0, y_max * 1.05)) +
  theme_classic(base_size = 14, base_family = "Times New Roman") +
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
    axis.title = element_blank(),
    legend.position = "none"
  )

pp

ggsave(
  filename = "transient_biomass_plot_clean1.tiff",
  plot = pp,
  path = "C:/Users/Research Greenhouse/Desktop/paper goudie/goudie2/data",
  dpi = 600,
  width = 7,
  height = 5,
  units = "in",
  compression = "lzw"
)


############################################################
## 10. R2 per treatment curve
############################################################

df_10 <- subset(biomass_richness, trt == "10 m")
df_15 <- subset(biomass_richness, trt == "15 m")
df_20 <- subset(biomass_richness, trt == "20 m")

m10 <- lm(AGB_total ~ log_transient, data = df_10)
m15 <- lm(AGB_total ~ log_transient, data = df_15)
m20 <- lm(AGB_total ~ log_transient, data = df_20)

r2_10 <- summary(m10)$r.squared
r2_15 <- summary(m15)$r.squared
r2_20 <- summary(m20)$r.squared

r2_table <- data.frame(
  Treatment = c("10 m", "15 m", "20 m"),
  R2 = c(r2_10, r2_15, r2_20)
)

r2_table

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