# main_template.R

cat("Working directory:", getwd(), "\n")

output_dir <- file.path(getwd(), "outputs/script18")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cat("Created directory:", output_dir, "\n")

source("DT_Algorithm.R")

the_ratio <- c(0, 0.25, 0.5, 0.75, 1)

for(qq in the_ratio){
  
  results_DT <- data.frame(
    Run = integer(),
    Length = numeric(),
    Cost = numeric(),
    NumDisambigs = integer(),
    Path = integer()
  )
  
  for (i in 1:100) {
    
    obs_gen_para <- read.csv(paste0("pattern/CSR", qq, "/obs_info_all_", qq, "_", i, ".csv"))
    
    result <- DT_fixed_Alg(obs_gen_para, 18)
    
    results_DT[i, ] <- list(
      Run = i,
      Length = result$Length_total,
      Cost = result$Cost_total,
      NumDisambigs = length(result$Disambiguate_state),
      Path = list(result$Optimal_path)
    )
  }
  
  results <- results_DT
  
  results_out <- data.frame(
    Index = paste0('"', 1:nrow(results), '"'),
    results[, c("Length", "Cost", "NumDisambigs")]
  )
  
  header <- '"length" "cost" "number_of_disambiguations"'
  
  txt_path <- file.path(output_dir, paste0("results_DT18_CSR_", qq, ".txt"))
  
  writeLines(header, txt_path)
  
  write.table(
    results_out,
    file = txt_path,
    append = TRUE,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    sep = " "
  )
  
  file_name <- file.path(output_dir, paste0("results_DT18_CSR_", qq, ".rds"))
  saveRDS(results[, "Path"], file = file_name)
  
  
}
