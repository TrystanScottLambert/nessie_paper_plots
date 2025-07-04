"""
Comparing group properties found by Nessie versus the published ones.
"""
import pandas as pd
import pylab as plt
import numpy as np
import plotting


def histogram_plot(nessie: pd.DataFrame, other: pd.DataFrame, column: str, **kwargs) -> None:
    """Compares the two data frames with the given column name"""
    plt.hist(nessie[column], lw=3, color='r', histtype='step', **kwargs, label='Nessie')
    plt.hist(other[column], lw=2, color='k', histtype='step', **kwargs, label='Tempel+2017')

    plt.legend()

def scatter_plot(nessie: pd.DataFrame, other:pd.DataFrame) -> None:
    """
    Scatter plot between nessie and another survey
    """

    plt.scatter(nessie['ra'], nessie['dec'], c = nessie['redshift'], s = 10*nessie['multiplicity'])
    plt.scatter(other['ra'], other['dec'], marker='+', c='r', alpha=0.5)
    plt.show()


if __name__ == '__main__':
    nessie_df = pd.read_csv('nessie_on_sdss.csv')
    sdss_df = pd.read_csv('sdss_groups_calculated.csv')

    #Multiplicity
    plotting.start_plot("Multiplicity", "log(counts)")
    histogram_plot(nessie_df, sdss_df, "multiplicity", bins=np.arange(5, 100))
    plt.yscale('log')
    plotting.end_plot("sdss_comp_multiplicity.png")
    plt.show()

    #N(z)
    plotting.start_plot("Redshift", "Counts")
    histogram_plot(nessie_df[nessie_df['multiplicity'] >=5], sdss_df[sdss_df['multiplicity'] >=5], "redshift", bins = np.arange(0, 0.5, 0.01))
    plt.xlim(0, 0.3)
    plotting.end_plot("sdss_comp_n_of_z.png")
    plt.show()

    scatter_plot(nessie_df, sdss_df)
