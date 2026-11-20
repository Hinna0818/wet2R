qpcr_measured_data <- function() {
  utils::read.csv(system.file("extdata", "qpcr_hsp12.csv", package = "wet2R", mustWork = TRUE),
                   colClasses = c(group = "character", bio_rep = "character", tech_rep = "character"))
}
