library(PTXQC)
require(PTXQC)
require(methods)

txt_folder = "C:/Users/kqjf682/Omics_pipelines/phospho/MQ/txt"

r = createReport(txt_folder)

cat(paste0("\nReport generated as '", r$report_file, "'\n\n"))
