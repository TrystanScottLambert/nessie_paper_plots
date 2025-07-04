remove.packages('Nessie')
devtools::install_local('~/Desktop/my_tools/Nessie')

library(Nessie)
library(data.table)

cosmo <- Nessie::FlatCosmology$new(0.7, 0.3)

random_zs <- g09_randoms <- as.data.table(read.csv("~/Desktop/GAMA_paper_plotter/gama_g09_randoms.txt"))
func <- Nessie::create_density_function(g09_randoms$z, length(g09_randoms$z)/400, 0.001453924, cosmo)

calibration_data <- as.data.frame(arrow::read_parquet("~/Desktop/GAMA_paper_plotter/mocks/galform_gals_for_R.parquet"))
data <- calibration_data[calibration_data$Volume == 1, ]
print('data read in')

red_cat <- Nessie::RedshiftCatalog$new(data$RA, data$DEC, data$Zspec, func, cosmo)
red_cat$run_fof(b0 = 0.05, r0 = 18)

red_cat$completeness <- rep(0.95, length(red_cat$ra_array))
red_cat$mock_group_ids <- data$GroupID

score = red_cat$compare_to_mock(min_group_size = 5)
print(score)

another_score = Nessie::calculate_s_score(red_cat$group_ids, red_cat$mock_group_ids, 5)
print(another_score)
