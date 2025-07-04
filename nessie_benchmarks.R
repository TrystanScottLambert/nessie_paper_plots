remove.packages('Nessie')
#devtools::install_github('TrystanScottLambert/Nessie')
devtools::install_local('~/Desktop/my_tools/Nessie')
library(Nessie)
library(celestial)
library(ggplot2)
library(reshape2)

infile_randoms <- '~/Desktop/GAMA_paper_plotter/gama_g09_randoms.txt'
randoms <- read.csv(infile_randoms)
random_z <- randoms$z

read_gama_field <- function(gama_field) {
    stub <- '~/Desktop/GAMA_paper_plotter/gama_galaxy_catalogs/g'
    infile <- paste(stub, gama_field, "_galaxies.dat", sep = '')
    gama <- read.csv(infile, sep = ' ')
    gama <- gama[gama$Z < 0.5, ]
    gama <- gama[gama$Rpetro < 19.65, ]
    gama['ab_mag'] <- gama$Rpetro - cosdistDistMod(gama$Z)
    gama['vel_err'] <- rep(50., length(gama$RA))
    return(gama)
}

# Function to shift coordinates
shift_coordinates <- function(gama_data, ra_shift = 0, dec_shift = 0) {
    shifted <- gama_data
    shifted$RA <- shifted$RA + ra_shift
    shifted$DEC <- shifted$DEC + dec_shift

    # Handle RA wrapping (0-360 degrees)
    shifted$RA <- shifted$RA %% 360

    # Handle DEC bounds (-90 to +90 degrees)
    shifted$DEC <- pmax(-90, pmin(90, shifted$DEC))

    return(shifted)
}

# Read original fields
gama_09 <- read_gama_field("09")
gama_12 <- read_gama_field("12")
gama_15 <- read_gama_field("15")
gama_23 <- read_gama_field("23")

# Create shifted versions
# Shift up by 30 degrees in DEC
gama_09_up <- shift_coordinates(gama_09, dec_shift = 30)
gama_12_up <- shift_coordinates(gama_12, dec_shift = 30)
gama_15_up <- shift_coordinates(gama_15, dec_shift = 30)
gama_23_up <- shift_coordinates(gama_23, dec_shift = 30)

# Shift around by 180 degrees in RA
gama_09_around <- shift_coordinates(gama_09, ra_shift = 180)
gama_12_around <- shift_coordinates(gama_12, ra_shift = 180)
gama_15_around <- shift_coordinates(gama_15, ra_shift = 180)
gama_23_around <- shift_coordinates(gama_23, ra_shift = 180)

# Shift both up and around
gama_09_up_around <- shift_coordinates(gama_09, ra_shift = 180, dec_shift = 30)
gama_12_up_around <- shift_coordinates(gama_12, ra_shift = 180, dec_shift = 30)
gama_15_up_around <- shift_coordinates(gama_15, ra_shift = 180, dec_shift = 30)
gama_23_up_around <- shift_coordinates(gama_23, ra_shift = 180, dec_shift = 30)

# Create extended combinations
gamas_extended <- list(
    # Original combinations
    gama_1 = gama_09,
    gama_2 = rbind(gama_09, gama_12),
    gama_3 = rbind(gama_09, gama_12, gama_15),
    gama_4 = rbind(gama_09, gama_12, gama_15, gama_23),

    # Add shifted up versions
    gama_5 = rbind(gama_09, gama_12, gama_15, gama_23, gama_09_up),
    gama_6 = rbind(gama_09, gama_12, gama_15, gama_23, gama_09_up, gama_12_up),
    gama_7 = rbind(gama_09, gama_12, gama_15, gama_23, gama_09_up, gama_12_up, gama_15_up),
    gama_8 = rbind(gama_09, gama_12, gama_15, gama_23, gama_09_up, gama_12_up, gama_15_up, gama_23_up),

    # Add shifted around versions
    gama_9 = rbind(gama_09, gama_12, gama_15, gama_23, gama_09_up, gama_12_up, gama_15_up, gama_23_up,
                   gama_09_around),
    gama_10 = rbind(gama_09, gama_12, gama_15, gama_23, gama_09_up, gama_12_up, gama_15_up, gama_23_up,
                    gama_09_around, gama_12_around),
    gama_11 = rbind(gama_09, gama_12, gama_15, gama_23, gama_09_up, gama_12_up, gama_15_up, gama_23_up,
                    gama_09_around, gama_12_around, gama_15_around),
    gama_12 = rbind(gama_09, gama_12, gama_15, gama_23, gama_09_up, gama_12_up, gama_15_up, gama_23_up,
                    gama_09_around, gama_12_around, gama_15_around, gama_23_around),

    # Add all shifted versions (up + around)
    gama_13 = rbind(gama_09, gama_12, gama_15, gama_23, gama_09_up, gama_12_up, gama_15_up, gama_23_up,
                    gama_09_around, gama_12_around, gama_15_around, gama_23_around, gama_09_up_around),
    gama_14 = rbind(gama_09, gama_12, gama_15, gama_23, gama_09_up, gama_12_up, gama_15_up, gama_23_up,
                    gama_09_around, gama_12_around, gama_15_around, gama_23_around,
                    gama_09_up_around, gama_12_up_around),
    gama_15 = rbind(gama_09, gama_12, gama_15, gama_23, gama_09_up, gama_12_up, gama_15_up, gama_23_up,
                    gama_09_around, gama_12_around, gama_15_around, gama_23_around,
                    gama_09_up_around, gama_12_up_around, gama_15_up_around),
    gama_16 = rbind(gama_09, gama_12, gama_15, gama_23, gama_09_up, gama_12_up, gama_15_up, gama_23_up,
                    gama_09_around, gama_12_around, gama_15_around, gama_23_around,
                    gama_09_up_around, gama_12_up_around, gama_15_up_around, gama_23_up_around)
)

# Set up cosmology and density function
cosmo <- FlatCosmology$new(0.7, 0.3)
rho_mean <- create_density_function(random_z, length(gama_09$RA), 0.001453924, cosmo)

# Run timing analysis
timings <- numeric(length(gamas_extended))
names(timings) <- names(gamas_extended)

for (i in seq_along(gamas_extended)) {
    cat("Running FOF on", names(gamas_extended)[i], "with", nrow(gamas_extended[[i]]), "galaxies\n")
    g <- gamas_extended[[i]]
    start.now <- Sys.time()
    catalog <- RedshiftCatalog$new(g$RA, g$DEC, g$Z, rho_mean, cosmo)
    catalog$get_raw_groups(0.06, 18)
    end.now <- Sys.time()
    timings[i] <- as.numeric(difftime(end.now, start.now, units = "secs"))
}

# Construct data frame
df <- data.frame(
    galaxies = sapply(gamas_extended, nrow),
    time = timings
)

# Add theoretical scaling curves - normalize to start at the first data point
first_galaxies <- df$galaxies[1]
first_time <- df$time[1]

# Calculate scaling factors based on first data point
n2_scale <- first_time / (first_galaxies^2)
nlogn_scale <- first_time / (first_galaxies * log(first_galaxies))
n_scale <- first_time / first_galaxies

df$n2 <- df$galaxies^2 * n2_scale
df$nlogn <- df$galaxies * log(df$galaxies) * nlogn_scale
df$n <- df$galaxies * n_scale

# Melt for plotting
df_melt <- melt(df, id.vars = "galaxies", measure.vars = c("time", "n2", "nlogn", "n"))

# Create the plot
p1 <- ggplot(df_melt, aes(x = galaxies, y = value, color = variable, linetype = variable)) +
    geom_line(size = 1.2) +
    geom_point(data = df_melt[df_melt$variable == "time", ], size = 3) +
    scale_color_manual(values = c("time" = "#D55E00", "n2" = "#0072B2", "nlogn" = "#009E73", "n" = "#CC79A7"),
                       labels = c("Actual FOF Time", expression(n^2), expression(n %.% log(n)), expression(n))) +
    scale_linetype_manual(values = c("time" = "solid", "n2" = "dashed", "nlogn" = "dotted", "n" = "dotdash"),
                          labels = c("Actual FOF Time", expression(n^2), expression(n %.% log(n)), expression(n))) +
    labs(
        title = "Extended FOF Runtime vs Galaxy Count with Complexity Curves",
        x = "Number of Galaxies",
        y = "Time (seconds)",
        color = "Legend",
        linetype = "Legend"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom")

# Create log-log plot for better visualization of scaling
p2 <- ggplot(df_melt, aes(x = galaxies, y = value, color = variable, linetype = variable)) +
    geom_line(size = 1.2) +
    geom_point(data = df_melt[df_melt$variable == "time", ], size = 3) +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values = c("time" = "#D55E00", "n2" = "#0072B2", "nlogn" = "#009E73", "n" = "#CC79A7"),
                       labels = c("Actual FOF Time", expression(n^2), expression(n %.% log(n)), expression(n))) +
    scale_linetype_manual(values = c("time" = "solid", "n2" = "dashed", "nlogn" = "dotted", "n" = "dotdash"),
                          labels = c("Actual FOF Time", expression(n^2), expression(n %.% log(n)), expression(n))) +
    labs(
        title = "Extended FOF Runtime vs Galaxy Count (Log-Log Scale)",
        x = "Number of Galaxies (log scale)",
        y = "Time (seconds, log scale)",
        color = "Legend",
        linetype = "Legend"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom")

# Print both plots
print(p1)
print(p2)

# Print summary statistics
cat("\nSummary of galaxy counts and timings:\n")
print(df[, c("galaxies", "time")])

# Calculate and print scaling exponent
log_df <- df[df$time > 0, ]  # Remove any zero times
if (nrow(log_df) > 1) {
    fit <- lm(log(time) ~ log(galaxies), data = log_df)
    scaling_exponent <- coef(fit)[2]
    cat("\nEstimated scaling exponent:", round(scaling_exponent, 3), "\n")
    cat("(1.0 = linear, 1.5 = n*log(n), 2.0 = quadratic)\n")
}
