"""
Comparing GAMA catalogues Nessie vs old values.
"""

import pandas as pd
import pylab as plt
import numpy as np
import plotting


def read_in_gama_cats(file_name: str) -> pd.DataFrame:
    """
    Reading in the group catalogs that were made using Nessie
    """
    return pd.read_csv(file_name)


def histogram_plot(
    nessie: pd.DataFrame, other: pd.DataFrame, column: str, **kwargs
) -> None:
    """Compares the two data frames with the given column name"""
    plt.hist(nessie[column], lw=3, color="r", histtype="step", **kwargs, label="Nessie")
    plt.hist(other[column], lw=2, color="k", histtype="step", **kwargs, label="G3Cv10")
    plt.legend()


def scatter_plot(nessie: pd.DataFrame, other: pd.DataFrame) -> None:
    """
    Scatter plot between nessie and another survey
    """

    plt.scatter(
        nessie["ra"], nessie["dec"], c=nessie["redshift"], s=10 * nessie["multiplicity"]
    )
    plt.scatter(other["ra"], other["dec"], marker="+", c="r", alpha=0.5)
    plt.show()


def convert_polar_to_cartesian(ra_array: np.ndarray, redshift_array: np.ndarray, convert_rad:bool = True):
    """
    poo poo pee pee
    """
    faq = np.pi/180
    if convert_rad is False:
        faq = 1

    x = np.cos(ra_array * faq) * redshift_array
    y = np.sin(ra_array * faq) * redshift_array
    return x, y

def wedge_scatter(
    data_frame: pd.DataFrame,
    central_ra,
    z_limits: tuple[float, float],
    dec_limits: tuple[float, float],
    ra_limits: tuple[float, float],
    **kwargs
) -> None:
    """
    Creates a wedge plot for the dataframe centralizing with 0 upwards
    """

    df = data_frame[
        (data_frame["dec"] >= dec_limits[0])
        & (data_frame["dec"] < dec_limits[1])
        & (data_frame["redshift"] > z_limits[0])
        & (data_frame["redshift"] < z_limits[1])
        & (data_frame['ra'] > ra_limits[0])
        & (data_frame['ra'] < ra_limits[1])
    ]
    df["ra"] = df["ra"] - central_ra + 90
    x, y = convert_polar_to_cartesian(df['ra'], df['redshift'])

    n=1000
    top_line = max(df['redshift']) + 0.001
    bottom_line = max(min(df['redshift']) - 0.001, 0)
    right_line = max(df['ra']) + 0.1
    left_line = min(df['ra']) - 0.1

    plt.scatter(x, y, s=10 * df["multiplicity"], **kwargs)
    
    top_line_redshifts = np.ones(n) * top_line
    top_line_ras = np.linspace(left_line, right_line, n)

    bottom_line_redshifts = np.ones(n) * bottom_line
    bottom_line_ras = np.linspace(left_line, right_line, n)

    left_line_ras = np.ones(n) * left_line
    left_line_redshifts = np.linspace(bottom_line, top_line, n)

    right_line_ras = np.ones(n) * right_line
    right_line_redshifts = np.linspace(bottom_line, top_line, n)

    top_x, top_y = convert_polar_to_cartesian(top_line_ras, top_line_redshifts)
    bottom_x, bottom_y = convert_polar_to_cartesian(bottom_line_ras, bottom_line_redshifts)
    left_x, left_y = convert_polar_to_cartesian(left_line_ras, left_line_redshifts)
    right_x, right_y = convert_polar_to_cartesian(right_line_ras, right_line_redshifts)

    plt.plot(top_x, top_y, color='k', lw=1.5)
    plt.plot(bottom_x, bottom_y, color='k', lw=1.5)
    plt.plot(left_x, left_y, color='k', lw=1.5)
    plt.plot(right_x, right_y, color='k', lw=1.5)


if __name__ == "__main__":
    nessie_df = pd.read_csv("nessie_on_gama.csv")
    gama_df = pd.read_csv("g_v10.csv")

    # Multiplicity
    plotting.start_plot("Multiplicity", "log(counts)")
    histogram_plot(nessie_df, gama_df, "multiplicity", bins=np.arange(5, 100))
    plt.yscale("log")
    plotting.end_plot("gama_comp_multiplicity.png")
    plt.show()

    # N(z)
    plotting.start_plot("Redshift", "Counts")
    histogram_plot(
        nessie_df[nessie_df["multiplicity"] >= 5],
        gama_df[gama_df["multiplicity"] >= 5],
        "redshift",
        bins=np.arange(0, 0.5, 0.01),
    )
    plotting.end_plot("gama_comp_n_of_z.png")
    plt.show()

    scatter_plot(nessie_df, gama_df)

    zlims = (0.1, 0.3)
    dec_lim = (0, 90)
    ra_lim = (200, 300)
    cen_ra = np.mean(gama_df[(gama_df['ra'] > ra_lim[0]) & (gama_df['ra'] < ra_lim[1])]['ra'])

    fig = plt.figure(figsize=(4.8, 2*3.54), dpi=600)
    wedge_scatter(nessie_df[nessie_df["multiplicity"] >= 3],cen_ra, zlims, dec_lim, ra_lim, facecolor='none', edgecolor='r', label='Nessie')
    wedge_scatter(gama_df[gama_df["multiplicity"] >=3], cen_ra, zlims, dec_lim, ra_lim, color = 'k', alpha=0.5, label= "G3Cv10")
    #plt.legend(loc=4, frameon=False, markerscale=0.5, fontsize=12)
    plt.text(-0.03, 0.1, r"$\delta \geq 0$", fontsize=15)
    plt.axis("off")
    plotting.end_plot("middle_wedge_positive_gama_15.png")
    plt.show()

    dec_lim = (-90, 0)
    fig = plt.figure(figsize=(4.8, 2*3.54), dpi=600)
    wedge_scatter(nessie_df[nessie_df["multiplicity"] >= 3],cen_ra, zlims, dec_lim, ra_lim, facecolor='none', edgecolor='r', label='Nessie')
    wedge_scatter(gama_df[gama_df["multiplicity"] >=3], cen_ra, zlims, dec_lim, ra_lim, color = 'k', alpha=0.5, label= "G3Cv10")
    plt.axis("off")
    plt.legend(loc=4, frameon=False, markerscale=0.5, fontsize=12)
    plt.text(-0.03, 0.1, r"$\delta < 0$", fontsize=15)
    plotting.end_plot("middle_wedge_negative_gama_15.png")
    plt.show()
