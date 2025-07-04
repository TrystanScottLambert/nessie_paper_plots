#remove.packages('Nessie')
#devtools::install_local('~/Desktop/my_tools/Nessie')

library(Nessie)
library(celestial)
library(data.table)
library(ggplot2)

cosmo <- Nessie::FlatCosmology$new(0.7, 0.3)

# Read in SDSS data
gama_data <- as.data.frame(read.csv("gama_gals.csv"))
colnames(gama_data) = c("ra", "dec", "redshift", "group_id", "ap_mag")
gama_data$group_id[gama_data$group_id == 0] <- -1
gama_data <- gama_data[gama_data$ra > 100, ] # Ignore gama02 region
abs_mag <- gama_data$ap_mag - cosmo$dist_mod(gama_data$redshift)
velocity_errors <- rep(50, length(abs_mag))
completeness <- rep(0.98, length(abs_mag))

g09_randoms <- as.data.table(read.csv("~/Desktop/GAMA_paper_plotter/gama_g09_randoms.txt"))
g12_randoms <- as.data.table(read.csv("~/Desktop/GAMA_paper_plotter/gama_g12_randoms.txt"))
g15_randoms <- as.data.table(read.csv("~/Desktop/GAMA_paper_plotter/gama_g15_randoms.txt"))

func <- Nessie::create_density_function(g09_randoms$z, (length(g09_randoms$z)/400) * 3, 0.004361773, cosmo)

red_cat <- Nessie::RedshiftCatalog$new(gama_data$ra, gama_data$dec, gama_data$redshift, func, cosmo)
red_cat$completeness <- completeness
start.now <- Sys.time()
red_cat$run_fof(b0 = 0.05, r0 = 18)
start.end <- Sys.time()
print(start.end - start.now)

red_cat$mock_group_ids <- gama_data$group_id


score = red_cat$compare_to_mock(min_group_size = 5)
print(score)

new_group_catalog <- red_cat$calculate_group_table(abs_mag, velocity_errors)
old_cat <- Nessie::RedshiftCatalog$new(gama_data$ra, gama_data$dec, gama_data$redshift, func, cosmo)
old_cat$group_ids <- as.integer(gama_data$group_id)
old_group_catalog <- old_cat$calculate_group_table(abs_mag, velocity_errors)

# Write the two catalogs to file.
write.csv(new_group_catalog, "nessie_on_gama.csv", row.names = F)
write.csv(old_group_catalog, "g_v10.csv", row.names = F)






