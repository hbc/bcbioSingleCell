library(bcbioSinglecell)
bcb <- load_run(
    upload_dir = file.path("data", "final"),
    sample_metadata_file = file.path("meta", "sample_metadata.xlsx"),
    interesting_groups = c("genotype", "treatment"),
    experiment_name = "",
    principal_investigator = "",
    email = "")
save_data(bcb)
