suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(here)
})

# Input Paths
drivers_path   <- here("export", "01_top_hepatocyte_lineage_drivers.csv")
occupancy_path <- here("export", "03_motif_occupancy_matrix_p1_p8.csv")

# Load
drivers   <- read_csv(drivers_path, show_col_types = FALSE)
occupancy <- read_csv(occupancy_path, show_col_types = FALSE)

# Synthesis
final_synthesis <- drivers %>%
  inner_join(
    occupancy %>% pivot_longer(-Peak_ID, names_to = "TF", values_to = "Has_Motif") %>% filter(Has_Motif == 1),
    by = "TF"
  ) %>%
  select(TF, Peak_ID, combined_score, chromvar_z, dorothea_z) %>%
  arrange(Peak_ID, desc(combined_score))

write_csv(final_synthesis,
          here("export", "05_final_tf_peak_binding_synthesis.csv"))