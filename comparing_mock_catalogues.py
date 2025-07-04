"""
Comparing the galform mocks and the shark mocks
"""

import numpy as np
import polars as pl
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from astropy.coordinates import SkyCoord
import pyvista as pv
import pylab as plt

from plot_wedge_diagrams_3d import add_xyz
from plotting import start_plot, end_plot


cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

INFILE_SHARK = "~/Desktop/GAMA_paper_plotter/mocks/shark_for_R_0.parquet"
INFILE_GALFORM = "~/Desktop/GAMA_paper_plotter/mocks/galform_gals_for_R.parquet"
INFILE_BROKEN_SHARK = "gama_gals_for_R.parquet"
INFILE_ACTUAL_GALAXIES = "~/Desktop/GAMA_paper_plotter/gama_galaxy_catalogs/g09_galaxies.dat"


#Actual gaama g09 region.
df_gama = pl.read_csv(INFILE_ACTUAL_GALAXIES, separator=' ')
df_gama = df_gama.insert_column(-1, pl.Series("co_dist", cosmo.comoving_distance(df_gama['Z'])))
new_names = {"RA": 'ra', "DEC": "dec", "Z": "redshift"}
df_gama = df_gama.rename(new_names)
df_gama = add_xyz(df_gama)



# Broken abundance matching online
df_shark_broken = pl.read_parquet(INFILE_BROKEN_SHARK)
df_shark_broken = df_shark_broken.filter(df_shark_broken['lightcone_gamafield'] == '0g09')
df_shark_broken = df_shark_broken.filter(df_shark_broken['ra'] < 150)
df_shark_broken = df_shark_broken.insert_column(-1, pl.Series("co_dist", cosmo.comoving_distance(df_shark_broken['zobs'])))
df_shark_broken = add_xyz(df_shark_broken)
shark_group_ids_broken, shark_counts_broken = np.unique(df_shark_broken['GroupID'], return_counts=True)
singletons = set(shark_group_ids_broken[np.where(shark_counts_broken==1)])
df_shark_broken = df_shark_broken.with_columns([pl.when(pl.col("GroupID").is_in(singletons)).then(-1).otherwise(pl.col("GroupID")).alias("GroupID")])
biggest_group_id_broken = shark_group_ids_broken[np.where(shark_counts_broken == np.max(shark_counts_broken))]
df_shark_big_group_broken = df_shark_broken.filter(df_shark_broken['GroupID'] == biggest_group_id_broken)

# Current shark mocks
df_shark = pl.read_parquet(INFILE_SHARK)
df_shark = df_shark.filter(df_shark['ra'] < 150)

df_shark = df_shark.insert_column(-1, pl.Series("co_dist", cosmo.comoving_distance(df_shark['zobs'])))
df_shark = add_xyz(df_shark)
shark_group_ids, shark_counts = np.unique(df_shark['id_group_sky'], return_counts=True)
shark_group_ids, shark_counts = shark_group_ids[1:], shark_counts[1:]
biggest_group_id = shark_group_ids[np.where(shark_counts == np.max(shark_counts))]
df_shark_big_group = df_shark.filter(df_shark['id_group_sky'] == biggest_group_id)


# Count group sizes
shark_group_ids, shark_counts = np.unique(df_shark['id_group_sky'], return_counts=True)

# Remove group ID 0 if needed
shark_group_ids, shark_counts = shark_group_ids[1:], shark_counts[1:]

# Get the top 20 group IDs by size
top_20_indices = np.argsort(shark_counts)[-20:]
top_20_group_ids = shark_group_ids[top_20_indices]

# Filter DataFrame to only include those groups
df_shark_big_groups = df_shark.filter(pl.col('id_group_sky').is_in(top_20_group_ids))


# Old galform mocks
df_galform = pl.read_parquet(INFILE_GALFORM)
df_galform = df_galform.filter(df_galform['RA'] < 150)
df_galform = df_galform.filter(df_galform['Volume'] == 1)
mapping = {'RA': 'ra', 'DEC': 'dec', 'Zspec': 'redshift'}
df_galform = df_galform.rename(mapping)
df_galform.insert_column(-1, pl.Series("co_dist", cosmo.comoving_distance(df_galform['redshift'])))
df_galform = add_xyz(df_galform)
galform_group_ids, galform_counts = np.unique(df_galform['GroupID'].to_numpy(), return_counts=True)
galform_group_ids, galform_counts = galform_group_ids[1:], galform_counts[1:]
biggest_galform_id = galform_group_ids[np.where(galform_counts == np.max(galform_counts))]
df_galform_big_group = df_galform.filter(df_galform['GroupID'] == biggest_galform_id)

galform_group_ids, galform_counts = np.unique(df_galform['GroupID'].to_numpy(), return_counts=True)
galform_group_ids, galform_counts = galform_group_ids[1:], galform_counts[1:]

top_20_indices = np.argsort(galform_counts)[-20:]
top_20_group_ids = galform_group_ids[top_20_indices]

df_galform_big_groups = df_galform.filter(pl.col('GroupID').is_in(top_20_group_ids))


df_shark_isolated = df_shark.filter((df_shark['id_group_sky'] == -1))
df_galform_isolated = df_galform.filter((df_galform['GroupID'] == -1))
df_shark_broken_isolated = df_shark_broken.filter((df_shark_broken['GroupID'] == -1) & (df_shark_broken['zobs'] < 0.1))

# Print the numbers first
isolated_galaxies = (len(np.where(df_galform['GroupID'].to_numpy() == -1)[0]), len(np.where(df_shark['id_group_sky'] == -1)[0]))
total_galaxies = (len(df_galform), len(df_shark))
#print(f"total: {total_galaxies[0]} {total_galaxies[1]} {total_galaxies[0]/total_galaxies[1]}")
#print(f"isolated: {isolated_galaxies[0]} {isolated_galaxies[1]} {isolated_galaxies[0]/isolated_galaxies[1]}")

pvpl = pv.Plotter(shape=(1, 3))
pvpl.subplot(0, 0)
pvpl.add_mesh(pv.PointSet(df_shark.select(["X", "Y", "Z"]).to_numpy()), style='points', color='k', point_size=2)
#pvpl.add_mesh(pv.PointSet(df_shark_isolated.select(["X", "Y", "Z"]).to_numpy()), style='points', color='r', point_size=3)
pvpl.add_mesh(pv.PointSet(df_shark_big_groups.select(["X", "Y", "Z"]).to_numpy()), style='points', color='b', point_size=4)


pvpl.subplot(0, 1)
pvpl.add_mesh(pv.PointSet(df_gama.select(["X", "Y", "Z"]).to_numpy()), style = 'points', color='k', point_size=2)


pvpl.subplot(0, 2)
pvpl.add_mesh(pv.PointSet(df_galform.select(["X", "Y", "Z"]).to_numpy()), style='points', color='k', point_size=2)
#pvpl.add_mesh(pv.PointSet(df_galform_isolated.select(["X", "Y", "Z"]).to_numpy()), style='points', color='r', point_size=3)
pvpl.add_mesh(pv.PointSet(df_galform_big_groups.select(["X", "Y", "Z"]).to_numpy()), style='points', color='b', point_size=4)


pvpl.link_views()
pvpl.show()


# histogram the multiplicites
bins = np.arange(5, 30, 1)
start_plot("Multiplicity", "Counts")
plt.hist(shark_counts, bins=bins, lw=3, color='r', label='Shark', histtype='step')
plt.hist(galform_counts, bins=bins, lw=2, color='b', label='Galform', histtype='step')
plt.legend()
plt.yscale('log')
end_plot("shark_galform_mult_comp.png")


# Making histogram of satellite masses.
bins = np.arange(6, 11.5, 0.1)
df_shark_satelites = df_shark.filter((df_shark['id_group_sky'] != -1) & (df_shark['type'] == 1))

plt.hist(np.log10(df_shark_satelites['mstars_bulge'] + df_shark_satelites['mstars_disk']), bins=bins, lw=2, histtype='step')

# Making the 2pcf for the two mocks and the real data.

def angular_separation(ra1, dec1, ra2, dec2):
    """
    Calculate the angular separation between two points on the sky.
    """
    coord1 = SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree, frame='icrs')
    coord2 = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree, frame='icrs')
    return coord1.separation(coord2).degree

def generate_random_catalog(ra_min, ra_max, dec_min, dec_max, size):
    """
    Generate a random catalog of celestial coordinates.
    """
    ra_random = np.random.uniform(ra_min, ra_max, size)
    dec_random = np.random.uniform(dec_min, dec_max, size)
    return ra_random, dec_random

def two_point_correlation_function(ra, dec, ra_random, dec_random, bins):
    """
    Compute the two-point correlation function using the Landy-Szalay estimator.
    """
    # Convert RA and Dec to SkyCoord objects
    coords = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
    coords_random = SkyCoord(ra=ra_random*u.degree, dec=dec_random*u.degree, frame='icrs')
    
    # Initialize pair counts
    DD = np.zeros(len(bins) - 1)
    RR = np.zeros(len(bins) - 1)
    DR = np.zeros(len(bins) - 1)
    
    # Count data-data pairs
    for i in range(len(ra)):
        sep = coords[i].separation(coords).degree
        DD += np.histogram(sep, bins=bins)[0]
    
    # Count random-random pairs
    for i in range(len(ra_random)):
        sep = coords_random[i].separation(coords_random).degree
        RR += np.histogram(sep, bins=bins)[0]
    
    # Count data-random pairs
    for i in range(len(ra)):
        sep = coords[i].separation(coords_random).degree
        DR += np.histogram(sep, bins=bins)[0]
    
    # Landy-Szalay estimator
    nD = len(ra)
    nR = len(ra_random)
    factor = (nR * (nR - 1)) / (nD * (nD - 1))
    xi = (DD - 2 * DR + RR) / RR * factor
    
    return xi
plt.show()

column_order = df_galform.columns
galform_isolated = df_galform.filter(pl.col("GroupID") == -1).select(column_order)
galform_brightest = df_galform.filter(pl.col("GroupID") != -1).sort("Rpetro", descending=True).group_by("GroupID").first().select(column_order)
galform_centrals = pl.concat([galform_isolated, galform_brightest])

redshift_binwidth = 0.05
angular_bins = np.arange(0, 0.02, 0.001)
plotting_angular_bins = (angular_bins[1:] + angular_bins[:-1])/2
for redshift in np.arange(0.1, 0.4, 0.05):
    shark_df_local = df_shark.filter((df_shark['zobs'] > redshift - redshift_binwidth) & (df_shark['zobs'] < redshift + redshift_binwidth)) #& (df_shark['type'] == 0))
    galform_df_local = df_galform.filter((df_galform['redshift'] < redshift + redshift_binwidth) & (df_galform['redshift'] > redshift - redshift_binwidth))
    gama_df_local = df_gama.filter((df_gama['redshift'] < redshift + redshift_binwidth) & (df_gama['redshift'] > redshift - redshift_binwidth))
    np.random.seed(777)
    ra_random_s, dec_random_s = generate_random_catalog(129, 141, -2, 3, len(shark_df_local))
    ra_random_g, dec_random_g = generate_random_catalog(129, 141, -2, 3, len(galform_df_local))
    ra_random_gama, dec_random_gama = generate_random_catalog(129, 141, -2, 3, len(gama_df_local))

    tp_shark = two_point_correlation_function(shark_df_local['ra'].to_numpy(), shark_df_local['dec'].to_numpy(), ra_random_s, dec_random_s, angular_bins)
    tp_galform = two_point_correlation_function(galform_df_local['ra'].to_numpy(), galform_df_local['dec'].to_numpy(), ra_random_g, dec_random_g, angular_bins)
    tp_gama = two_point_correlation_function(gama_df_local['ra'].to_numpy(), gama_df_local['dec'].to_numpy(), ra_random_gama, dec_random_gama, angular_bins)

    start_plot("Angular Separation [deg]", "TPCF")
    plt.plot(plotting_angular_bins, tp_shark, label='Shark', lw=3)
    plt.plot(plotting_angular_bins, tp_galform, label='Galform', lw=2)
    plt.plot(plotting_angular_bins, tp_gama, label='GAMA', lw=1.5)
    plt.title(f"z = {round(redshift, 2)}")
    plt.ylim(0, 17)
    plt.legend()
    end_plot(f"tpcf_{round(redshift, 2)}.png")
