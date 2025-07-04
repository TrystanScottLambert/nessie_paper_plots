"""Three D visualization of the wedge plots"""

import numpy as np
import pyvista as pv
import polars as pr
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import pylab as plt

from plotting import start_plot, end_plot

def read_in_old_gama_galaxies(infile: str, cosmology: FlatLambdaCDM) -> pr.DataFrame:
    """
    Adds columns that we need for three D plotting
    """
    df = pr.read_csv(infile)
    distances = cosmology.comoving_distance(df['Z'].to_numpy())
    df.insert_column(-1, pr.Series("co_dist", distances))

    new_column_names = {'RA': 'ra', "Dec": 'dec', "GroupID": "group_id", "Z": "redshift"}
    df = df.rename(new_column_names)

    df = df.with_columns(group_id=pr.col("group_id").replace(0, -1))
    df = df.filter((df['ra'] > 100) & (df['ra'] < 150))
    return df



# Add Cartesian coordinates
def add_xyz(data_frame: pr.DataFrame) -> pr.DataFrame:
    """Calculate the cartesian xyz values and add them to the dataframe"""
    c = SkyCoord(
        ra=data_frame["ra"].to_numpy() * u.deg,
        dec=data_frame["dec"].to_numpy() * u.deg,
        distance=data_frame["co_dist"].to_numpy() * u.Mpc,
    )
    return data_frame.with_columns(
        [
            pr.Series("X", c.cartesian.x.value),
            pr.Series("Y", c.cartesian.y.value),
            pr.Series("Z", c.cartesian.z.value),
        ]
    )


def create_glyphs(data_frame: pr.DataFrame):
    """Building glyphs that can be used in the pyvista plotter"""

    points = np.array(data_frame.select(["X", "Y", "Z"]).to_numpy())
    #radii = data_frame["r100"].to_numpy()
    radii = data_frame["multiplicity"].to_numpy()

    # Normalize radii for visibility
    radii_scaled = radii / np.max(radii) * 20  # Adjust the 5 to scale glyph size

    point_cloud = pv.PolyData(points)
    point_cloud["radii"] = radii_scaled
    glyphs = point_cloud.glyph(scale="radii", geom=pv.Sphere(theta_resolution=12, phi_resolution=12))
    return glyphs



def get_isolated_galaxies(galaxy_data_frame: pr.DataFrame) -> pr.DataFrame:
    """
    Return a data frane of all the isolated galaxies
    """
    df = galaxy_data_frame.clone()
    df = df.filter(df['group_id'] == -1)
    return df

def get_group_galaxies(df: pr.DataFrame, min_group_size: int) -> pr.DataFrame:
    group_counts = (
        df.select([
            pr.col("group_id").value_counts().alias("counts")
        ])
        .unnest("counts")  # this will give 'group_id' and 'counts'
        .filter(pr.col("count") >= min_group_size)
        .select("group_id")
    )
    df = df.filter(pr.col("group_id").is_in(group_counts["group_id"]))
    df = df.filter(df['group_id'] != -1)
    return df


if __name__ == "__main__":
    cosmo = FlatLambdaCDM(Om0=0.3, H0=70)
    galaxies = pr.read_csv("group_catalogs/group_galaxies_g09.dat")
    new_names = ["ra", "dec", "redshift"]
    mapping = {galaxies.columns[i+1]: new_names[i] for i in range(3)}
    galaxies = galaxies.rename(mapping)
    galaxies.insert_column(-1, pr.Series("co_dist", cosmo.comoving_distance(galaxies['redshift'].to_numpy()).value))
    galaxies = add_xyz(galaxies)
    galaxies = get_isolated_galaxies(galaxies)

    galaxies_gama = read_in_old_gama_galaxies('version_10_gama_catalogs/gama_group_galaxies.csv', cosmo)
    galaxies_gama = add_xyz(galaxies_gama)
    galaxies_gama = get_isolated_galaxies(galaxies_gama)
    galaxy_points = galaxies_gama.select(["X", "Y", "Z"]).to_numpy()

    g09_groups_tuned = pr.read_csv("group_catalogs/group_catalog_tuned_g09.dat")
    g09_groups_tuned = add_xyz(g09_groups_tuned)

    g09_groups_old = pr.read_csv("group_catalogs/group_catalog_old_g09.dat")
    g09_groups_old = add_xyz(g09_groups_old)

    g09_tuned_glyphs = create_glyphs(g09_groups_tuned)
    g09_old_glyphs = create_glyphs(g09_groups_old)

    pl = pv.Plotter(shape=(1, 2))
    pl.subplot(0, 0)
    pl.add_mesh(g09_tuned_glyphs, color="blue", opacity=0.3)
    pl.add_mesh(pv.PointSet(galaxy_points), style='points', color='k', point_size=2)

    pl.subplot(0, 1)
    pl.add_mesh(g09_old_glyphs, color="red", opacity=0.3)
    pl.add_mesh(pv.PointSet(galaxy_points), style='points', color='k', point_size=2)

    pl.link_views()
    pl.show()

    # Plotting the multiplicity distributions
    mult_bin = np.arange(3, 50, 1)
    start_plot('Multiplicity', "Frequency")
    plt.hist(g09_groups_tuned['multiplicity'], bins=mult_bin, lw=3, color='blue', label='new: (0.05, 36)', histtype='step')
    plt.hist(g09_groups_old['multiplicity'], bins=mult_bin, lw=2, color='red', label='old: (0.06, 18)', histtype='step')
    plt.legend()
    plt.yscale('log')

    # Customizing ticks
    plt.tick_params(axis='both', which='major', labelsize=14, direction='in', length=6, width=2)
    plt.tick_params(axis='both', which='minor', labelsize=12, direction='in', length=4, width=1.5)
    end_plot("effects_of_links_on_multiplicity.png")

    pl = pv.Plotter()
