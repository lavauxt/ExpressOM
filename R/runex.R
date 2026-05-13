#run_bulk_pipeline(
#  ensembl_package_name = "EnsDb.Hsapiens.v107",
#  count_type = "salmon",
#  data_dir = "./data/counts",  
#  out_dir = "./results",
#  sample_table = "./data/sample_table.csv",
#  level = "IFNa",
#  base = "Unstim",
#  model = "~ condition"
#)

#run_bulk_pipeline(
#  count_type = "matrix",
#  matrix_file = "./data/HELIOS_raw_counts.csv",
#  sample_table = "./data/sample_table.csv",
  # ... other params ...
#)