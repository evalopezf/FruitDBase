args <- commandArgs(trailingOnly = TRUE)

json_files <- args[1:(length(args)-1)]
output_rds <- args[length(args)]

library(jsonlite)
library(dplyr)

lista_json <- list()

for (json_file in json_files) {

  if (file.exists(json_file)) {

    info <- fromJSON(json_file)

    sample <- basename(dirname(json_file))

    lista_json[[sample]] <- data.frame(
      sample = sample,
      n_processed = info$n_processed,
      n_pseudoaligned = info$n_pseudoaligned,
      p_pseudoaligned = info$p_pseudoaligned
    )
  }
}

report <- bind_rows(lista_json)

saveRDS(report, output_rds)
