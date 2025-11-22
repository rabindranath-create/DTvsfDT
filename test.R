# main_template.R

cat("Working directory:", getwd(), "\n")

output_dir <- file.path(getwd(), "outputs/script0.5")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cat("Created directory:", output_dir, "\n")

#########################

rm(list = ls())                 # Clear variables
while (!is.null(dev.list())) dev.off()  # Close all plots

#rm(list = ls())                 # Clear variables
#while (!is.null(dev.list())) dev.off()  # Close all plots
library("rstudioapi") # Load rstudioapi package
setwd(dirname(getActiveDocumentContext()$path)) # Set working directory to source file location
getwd() # Check updated working directory






#######################




source("DT_Algorithm.R")

the_ratio <- c(0)

for(qq in the_ratio){
  
  results_DT <- data.frame(
    Run = integer(),
    Length = numeric(),
    Cost = numeric(),
    NumDisambigs = integer(), 
    Path = integer()
  )
  
  for (i in 1:2) {
    
    obs_gen_para <- read.csv(paste0("pattern/CSR", qq, "/obs_info_all_", qq, "_", i, ".csv"))
    
    result <- DT_fixed_Alg(obs_gen_para, 0.5)
    
    results_DT[i, ] <- list(
      Run = i,
      Length = result$Length_total,
      Cost = result$Cost_total,
      NumDisambigs = length(result$Disambiguate_state),
      Path = list(result$Optimal_path)
    )
  }
  
  results <- results_DT
  
 
}



results_out <- data.frame(
  Index = paste0('"', 1:nrow(results), '"'),
  results[, c("Length", "Cost", "NumDisambigs")]
)

header <- '"length" "cost" "number_of_disambiguations"'

txt_path <- file.path(output_dir, paste0("results_DT0.5_CSR_", qq, ".txt"))

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

# Build the file path
file_name <- file.path(output_dir, paste0("results_DT0.5_CSR_", qq, ".rds"))

# Save the results
saveRDS(results[, "Path"], file = file_name)



# Specify the path to your .rds file
file_name <- file.path(output_dir, paste0("results_DT0.5_CSR_", qq, ".rds"))

# Load the RDS file
results_loaded <- readRDS(file_name)

# Check
str(results_loaded)




##########

######

how_many <- 0

n <- min(length(results_loaded[[1]]), length(results_loaded[[2]]))

for(i in 1:(n-1)){
  if( (results_loaded[[1]][i] == results_loaded[[2]][i]) && (results_loaded[[1]][i+1] == results_loaded[[2]][i+1])){
    how_many = how_many +  1
  }
}
