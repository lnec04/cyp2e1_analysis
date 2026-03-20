suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(here)
})

# Input Paths
drivers_path   <- here("github", "cyp2e1", "01_tf_activity", "tables", "top_hepatocyte_lineage_drivers.csv")
occupancy_path <- here("github", "cyp2e1", "03_motif_analysis", "tables", "motif_occupancy_matrix_p1_p8.csv")

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

write_csv(final_synthesis, here("github", "cyp2e1", "final_tf_peak_binding_synthesis.csv"))