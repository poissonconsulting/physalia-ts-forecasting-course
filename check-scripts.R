# check for tasks to do
findR::findRscript(pattern = "TODO", path = ".", case.sensitive = FALSE,
                   comments = TRUE) %>%
  filter(path_to_file != "./check-scripts.R")

# check for place markers
findR::findRscript(pattern = "HERE", path = ".", case.sensitive = TRUE,
                   comments = TRUE) %>%
  filter(path_to_file != "./check-scripts.R")

findR::findRscript(pattern = "DELETE", path = ".", case.sensitive = TRUE,
                   comments = TRUE) %>%
  filter(path_to_file != "./check-scripts.R")
