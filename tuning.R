remove.packages('Nessie')
devtools::install_local('~/Desktop/my_tools/Nessie')

library(Nessie)
library(psych)
library(Highlander)
library(celestial)
library(data.table)
library(arrow)
# Setting up cosmology
cosmo_shark <- Nessie::FlatCosmology$new(h = 0.6751, omega_matter = 0.3)
cosmo <- Nessie::FlatCosmology$new(h = 0.7, omega_matter = 0.3)

# Calculating rho mean from random catalogues.
g09_area <- skyarea(c(129, 141), c(-2, 3))
g12_area <- skyarea(c(174, 186), c(-3, 2))
g15_area <- skyarea(c(211.5, 223.5), c(-2, 3))
g23_area <- skyarea(c(339, 351), c(-35, -30))
combined_area = g09_area['areafrac'] + g12_area['areafrac'] + g15_area['areafrac'] + g23_area['areafrac']

g09_randoms <- as.data.table(read.csv("~/Desktop/GAMA_paper_plotter/gama_g09_randoms.txt"))
g12_randoms <- as.data.table(read.csv("~/Desktop/GAMA_paper_plotter/gama_g12_randoms.txt"))
g15_randoms <- as.data.table(read.csv("~/Desktop/GAMA_paper_plotter/gama_g15_randoms.txt"))
g23_randoms <- as.data.table(read.csv("~/Desktop/GAMA_paper_plotter/gama_g23_randoms.txt"))
gama_combined <- as.data.table(read.csv("~/Desktop/GAMA_paper_plotter/gama_combined_randoms.txt"))

g09_rho_mean <- Nessie::create_density_function(g09_randoms$z, length(g09_randoms$z)/400, g09_area["areafrac"], cosmo)
g12_rho_mean <- Nessie::create_density_function(g12_randoms$z, length(g12_randoms$z)/400, g12_area["areafrac"], cosmo)
g15_rho_mean <- Nessie::create_density_function(g15_randoms$z, length(g15_randoms$z)/400, g15_area["areafrac"], cosmo)
g23_rho_mean <- Nessie::create_density_function(g23_randoms$z, length(g23_randoms$z)/400, g23_area["areafrac"], cosmo)
gama_rho_mean <- Nessie::create_density_function(gama_combined$z, length(gama_combined$z)/400, combined_area, cosmo)
rho_means <- list(g09 = g09_rho_mean, g12 = g12_rho_mean, g15 = g15_rho_mean, g23 = g23_rho_mean)

# Setting up Redshift catalogues
rlim <- 19.65

calibration_data <- as.data.frame(arrow::read_parquet("~/Desktop/GAMA_paper_plotter/mocks/gama_gals_for_R.parquet"))


lightcone_numbers = c(0, 3, 4, 5, 7, 8, 9, 10)
gama_fields = c('g09', 'g12', 'g15', 'g23')
combinations <- expand.grid(lightcone_numbers, gama_fields)
all_fields <- paste0(combinations$Var1, combinations$Var2)

redshift_catalogues <- list()
for (g_field in all_fields) {
    gama_field <- substr(g_field, nchar(g_field) - 2, nchar(g_field))
    local_catalogue <- calibration_data[calibration_data['lightcone_gamafield'] == g_field, ]
    red_cat <- RedshiftCatalog$new(local_catalogue$ra, local_catalogue$dec, local_catalogue$zobs, rho_means[[gama_field]], cosmo)
    red_cat$mock_group_ids <- local_catalogue$GroupID
    red_cat$completeness <- rep(0.95, length(red_cat$ra_array))

    # setting this values to -1
    counts <- table(red_cat$mock_group_ids)
    singleton_ids <- names(counts[counts == 1])
    red_cat$mock_group_ids <- ifelse(red_cat$mock_group_ids %in% singleton_ids, -1, red_cat$mock_group_ids)

    redshift_catalogues[[length(redshift_catalogues) + 1]] <- red_cat
}



redshift_catalogues_shark <- list()
for (lc in lightcone_numbers) {
    shark_data <- as.data.frame(arrow::read_parquet(paste("~/Desktop/GAMA_paper_plotter/mocks/shark_for_R_",lc,".parquet", sep = '')))

    red_cat <- RedshiftCatalog$new(shark_data$ra, shark_data$dec, shark_data$zobs, gama_rho_mean, cosmo_shark)

    #testing different IDS things.
    #mock_ids <- shark_data$id_group_sky
    #counts <- table(mock_ids)
    #new_mock_ids <- ifelse(counts[as.character(mock_ids)] == 1, -1, mock_ids)
    #red_cat$mock_group_ids <- as.integer(new_mock_ids)

    red_cat$mock_group_ids <- shark_data$id_group_sky
    red_cat$completeness <- rep(0.95, length(red_cat$ra_array))

    redshift_catalogues_shark[[length(redshift_catalogues_shark) + 1]] <- red_cat
}


lc_numbers <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)
redshift_catalogues_galform <- list()
for (lc in lc_numbers) {
    galform_data <- as.data.frame(arrow::read_parquet("~/Desktop/GAMA_paper_plotter/mocks/galform_gals_for_R.parquet"))
    local_catalogue <- galform_data[galform_data$Volume == lc, ]

    red_cat <- RedshiftCatalog$new(local_catalogue$RA, local_catalogue$DEC, local_catalogue$Zspec, gama_rho_mean, cosmo)
    red_cat$mock_group_ids <- local_catalogue$GroupID
    red_cat$completeness <- rep(0.95, length(red_cat$ra_array))

    redshift_catalogues_galform[[length(redshift_catalogues_galform) + 1]] <- red_cat
}

red_cat_galform <- redshift_catalogues_galform[[1]]
red_cat_shark <- redshift_catalogues_shark[[1]]

red_cat_galform$run_fof(b0 = 0.05, r0 = 30)
red_cat_shark$run_fof(b0 = 0.05, r0 = 30)



# Trying out the broken abundance matched stuff.
broken <- arrow::read_parquet("gama_gals_for_R.parquet")
broken <- broken[broken$lightcone_gamafield == '2g15', ]
red_cat_broken <- RedshiftCatalog$new(broken$ra, broken$dec, broken$zobs, gama_rho_mean, cosmo)
red_cat_broken$run_fof(b0 = 0.05, r0 = 30)
broken_mock_ids <- broken$GroupID
counts <- table(broken_mock_ids)
new_broken_mock_ids <- ifelse(counts[as.character(broken_mock_ids)] == 1, -1, broken_mock_ids)
red_cat_broken$mock_group_ids <- as.integer(new_broken_mock_ids)
Nessie::calculate_s_score(red_cat_shark$group_ids, red_cat_shark$mock_group_ids, 5)
Nessie::calculate_s_score(red_cat_galform$group_ids, red_cat_galform$mock_group_ids, 5)
Nessie::calculate_s_score(red_cat_broken$group_ids, red_cat_broken$mock_group_ids, 5)

# Running the tuning
#Nessie::tune_group_finder(
#    redshift_catalogues_shark,
#    minimum_group_size = 5,
#    b0_estimate = 0.06,
#    r0_estimate = 36,
#    b0_bounds = c(0.04, 0.1),
#    r0_bounds = c(15, 50),
#    scaling = c(100, 1))

