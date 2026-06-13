library(tidyverse)
library(maps)

# Install maps if needed: install.packages("maps")

df <- read_csv("combined_edi_bg_results_20260422.csv",
               na = c("NA", ""))

# ── 1. Global map of dataset locations ────────────────────────────────────────

map_df <- df %>%
  filter(!is.na(bb_west), !is.na(bb_east), !is.na(bb_south), !is.na(bb_north)) %>%
  mutate(
    lon = (bb_west + bb_east) / 2,
    lat = (bb_south + bb_north) / 2,
    ecosystem_type = replace_na(ecosystem_type, "unknown") %>% str_to_title()
  )

ecosystem_colors <- c(
  "Forest"     = "#4daf4a", "Freshwater" = "#377eb8", "Wetlands"   = "#984ea3",
  "Tundra"     = "#a6cee3", "Grassland"  = "#ff7f00", "Ocean"      = "#1f78b4",
  "Mountain"   = "#b15928", "Desert"     = "#e6ab02", "Savannah"   = "#cab2d6",
  "Other"      = "#999999", "Unknown"    = "#cccccc"
)

world <- map_data("world")

p_map <- ggplot() +
  geom_polygon(data = world,
               aes(x = long, y = lat, group = group),
               fill = "gray90", color = "gray70", linewidth = 0.15) +
  geom_point(data = map_df,
             aes(x = lon, y = lat, color = ecosystem_type),
             size = 2, alpha = 0.75) +
  scale_color_manual(values = ecosystem_colors, name = "Ecosystem type") +
  coord_fixed(ratio = 1.3, xlim = c(-180, 180), ylim = c(-90, 90)) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_line(color = "gray95"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  ) +
  guides(color = guide_legend(nrow = 3, override.aes = list(size = 3))) +
  labs(
    title = "Global distribution of SABRE belowground datasets",
    subtitle = paste0(nrow(map_df), " of ", nrow(df), " datasets with geographic coordinates"),
    x = NULL, y = NULL
  )

# ── 2. Histogram of years of data available ───────────────────────────────────

years_df <- df %>%
  filter(!is.na(year_start), !is.na(year_end)) %>%
  mutate(years_of_data = pmax(year_end - year_start + 1, 1))

p_years <- ggplot(years_df, aes(x = years_of_data)) +
  geom_histogram(binwidth = 2, fill = "#2b7bba", color = "white", linewidth = 0.3) +
  scale_x_continuous(breaks = seq(0, max(years_df$years_of_data, na.rm = TRUE) + 5, by = 5)) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank()) +
  labs(
    title = "Temporal span of SABRE belowground datasets",
    subtitle = paste0(nrow(years_df), " datasets with start and end years recorded"),
    x = "Years of data available",
    y = "Number of datasets"
  )

# ── 3. Bar chart of belowground variables ─────────────────────────────────────

# Shorten a few long names for readability
shorten <- c(
  "microbial community composition (amplicon)"       = "Microbial community (amplicon)",
  "microbial community composition (ARISA)"          = "Microbial community (ARISA)",
  "microbial community composition (amplicon/metagenomics)" = "Microbial community (amplicon/metagenomics)",
  "microbial community diversity"                    = "Microbial community diversity",
  "microbial community composition"                  = "Microbial community composition",
  "PLFA (microbial community composition)"           = "PLFA (community composition)",
  "arbuscular mycorrhizal fungi biomass"             = "Arbuscular mycorrhizal fungi",
  "aquatic microbial/phytoplankton biomass"          = "Aquatic microbial/phytoplankton",
  "microbial mat biomass (AFDM)"                     = "Microbial mat biomass (AFDM)",
  "denitrification enzyme activity"                  = "Denitrification enzyme activity",
  "soil nitrogen mineralization"                     = "Soil N mineralization",
  "microbial biomass carbon"                         = "Microbial biomass C",
  "microbial biomass nitrogen"                       = "Microbial biomass N",
  "fungal:bacterial ratio (PLFA)"                    = "Fungal:bacterial ratio (PLFA)"
)

bg_vars <- df %>%
  filter(!is.na(belowground_vars)) %>%
  separate_rows(belowground_vars, sep = ";\\s*") %>%
  mutate(belowground_vars = str_trim(belowground_vars)) %>%
  filter(belowground_vars != "") %>%
  mutate(belowground_vars = recode(belowground_vars, !!!shorten)) %>%
  mutate(belowground_vars = str_to_sentence(belowground_vars)) %>%
  count(belowground_vars, sort = TRUE) %>%
  slice_head(n = 20)

p_bgvars <- ggplot(bg_vars, aes(x = n, y = reorder(belowground_vars, n))) +
  geom_col(fill = "#8c4a2f", color = NA) +
  geom_text(aes(label = n), hjust = -0.2, size = 3.2) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(size = 9)
  ) +
  labs(
    title = "Most common belowground variables in SABRE datasets",
    subtitle = paste0("Top 20 variables across ",
                      df %>% filter(!is.na(belowground_vars)) %>% nrow(),
                      " datasets with identified belowground variables"),
    x = "Number of datasets",
    y = NULL
  )

# ── 4. Bar chart of study types ───────────────────────────────────────────────

study_type <- df %>%
  filter(!is.na(study_type)) %>%
  mutate(study_type = str_trim(study_type)) %>%
  filter(study_type != "") %>%
  count(study_type, sort = TRUE)

p_studytype <- ggplot(study_type, aes(x = n, y = reorder(study_type, n))) +
  geom_col(fill = "#6a3d9a", color = NA) +
  geom_text(aes(label = n), hjust = -0.2, size = 3.2) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(size = 9)
  ) +
  labs(
    title = "Study types in SABRE datasets",
    subtitle = paste0(nrow(study_type), " unique study types across ",
                      df %>% filter(!is.na(study_type), str_trim(study_type) != "") %>% nrow(),
                      " datasets"),
    x = "Number of datasets",
    y = NULL
  )

# ── Save plots ────────────────────────────────────────────────────────────────

ggsave("sabre_map.png",      p_map,      width = 10, height = 6,   dpi = 150)
ggsave("sabre_years.png",    p_years,    width = 7,  height = 4.5, dpi = 150)
ggsave("sabre_bgvars.png",   p_bgvars,   width = 8,  height = 6,   dpi = 150)
ggsave("sabre_studytype.png", p_studytype, width = 8, height = 6,   dpi = 150)

message("Plots saved: sabre_map.png, sabre_years.png, sabre_bgvars.png, sabre_studytype.png")
