# Script to build the final group catalog.
#remove.packages('Nessie')
#devtools::install_github('TrystanScottLambert/Nessie')
library(Nessie)
library(data.table)
library(celestial)
library(FoF)
library(ggplot2)
library(dplyr)

calculate_k_corrections <- function(redshift_array) {
    k_corrections <- {}
    for (red in redshift_array) {
        k_corrections <-  c(k_corrections, FoF::KEcorr(red)[2])
    }
    return(k_corrections)
}

calculate_ab_mags <- function(redshift_array, apparent_mags, cosmology) {
    ke_corrs <- calculate_k_corrections(redshift_array)
    distance_moduli <- cosmology$dist_mod(redshift_array)
    return(apparent_mags - distance_moduli - ke_corrs)
}

# Best values from tuning.
b0_tuned <- 0.05
r0_tuned <- 36

b0_aaron <- 0.06
r0_aaron <- 18
cosmo <- Nessie::FlatCosmology$new(0.7, 0.3)


g09_area <- skyarea(c(129, 141), c(-2, 3))
g12_area <- skyarea(c(174, 186), c(-3, 2))
g15_area <- skyarea(c(211.5, 223.5), c(-2, 3))
g23_area <- skyarea(c(339, 351), c(-35, -30))
areas <- list(g09 = g09_area['areafrac'], g12 = g12_area['areafrac'], g15 = g15_area['areafrac'], g23 = g23_area['areafrac'])

create_group_catalog <- function(field, b0, r0) {
    random_file_name <- paste("~/Desktop/GAMA_paper_plotter/gama_",field,"_randoms.txt", sep = '')
    galaxy_data_file_name <- paste("~/Desktop/GAMA_paper_plotter/gama_galaxy_catalogs/",field,"_galaxies.dat", sep = '')

    random_catalog <- as.data.table(read.csv(random_file_name))
    rho_mean <- Nessie::create_density_function(random_catalog$z, length(random_catalog$z)/400, as.numeric(areas[field]), cosmo)
    galaxies <- as.data.frame(read.csv("~/Desktop/GAMA_paper_plotter/gama_galaxy_catalogs/g09_galaxies.dat", sep = ' '))
    red_cat <- RedshiftCatalog$new(galaxies$RA, galaxies$DEC, galaxies$Z, rho_mean, cosmo, completeness = rep(0.95, length(galaxies$RA)))

    absolute_mags <- calculate_ab_mags(galaxies$Z, galaxies$Rpetro, cosmo)
    vel_errs <- rep(50., length(galaxies$Z))

    red_cat$run_fof(b0, r0)
    galaxies['group_id'] <- red_cat$group_ids
    write.csv(galaxies, paste("group_catalogs/group_galaxies_", field,".dat", sep = ""), row.names = FALSE)

    group_catalog <- red_cat$calculate_group_table(absolute_mags, vel_errs)

    return(group_catalog)
}

save_file <- function(group_catalog, field, type) {
    outfile <- paste("group_catalogs/group_catalog_",type,"_", field, ".dat", sep = "")
    write.csv(group_catalog, outfile, row.names = FALSE)
}


fields <- c("g09", "g12", "g15", "g23")
for (field in fields) {
    group_cat <- create_group_catalog(field, b0_tuned, r0_tuned)
    save_file(group_cat, field, "tuned")
}

for (field in fields) {
    group_cat <- create_group_catalog(field, b0_aaron, r0_aaron)
    save_file(group_cat, field, "old")
}



