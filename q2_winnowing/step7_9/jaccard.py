

# <><><> SETUP IMPORTS <><><>
import os
import pandas as pd
import numpy as np
from rpy2.robjects.packages import importr
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
from rpy2.robjects import pandas2ri


# Setup R packages used
pandas2ri.activate()

utils = rpackages.importr('utils')
utils.chooseCRANmirror(ind=1)

packnames = ['irr']
names_to_install = []

for x in packnames:
    if (rpackages.isinstalled(x)==False):
        names_to_install.append(x)

if len(names_to_install) > 0:
    utils.install_packages(StrVector(names_to_install))

rirr = importr('irr')


# <><><> DEFINE FUNCTIONS <><><>

def _print_summaries( kappa_df, dump ):
    """
    Output table and column stats to dump
    :param kappa_df: Dataframe which is refrenced for prints
    :param dump: file where output is written to
    :return: nothing is returned function is simply for generating verbose output
    """

    kappa_df['select_iter'] = kappa_df['select_iter'].astype('int')
    # <><><> OUTPUT TABLE <><><>
    dump.write( str( kappa_df ) + "\n\n")

    # <><><> PRINTS STATS OF COLUMNS <><><>
    dump.write(str( kappa_df.groupby(['conditioning', 'centrality', 'correl', 'threshold'])['agreement'].mean() ) + "\n")
    dump.write(str( kappa_df.loc[( kappa_df["conditioning"] == "Add1") & ( kappa_df["correl"] == "MIC")].groupby(
        ['centrality', 'correl', 'threshold'])['agreement'].mean() ) + "\n")

    return  # finished print



def jaccard_coefficient(x,y):
    """
    Jaccard index used to display the similarity between sample sets.
    :param A: set x
    :param y: set y
    :return: similarity of A and B or A intersect B
    """
    return len(set(x) & set(y)) / len(set(x) | set(y)) # Set gets rid of duplicates


def main( leave_one_out_df, name, detailed=False, verbose=False ):
    """
    perform jaccard index on data to find most influential taxa.
    :param leave_one_out_df: result of leave one out feature ordering method. Ordered OTU
    :param name: name that is attached to output data for identification
    :param detailed: Output helper tables
    :param verbose: Output helper prints
    :return:
    """

    out_dir = f"{os.path.dirname(os.path.realpath(__file__))}/output"
    # allows for cleaner execution and use of relative paths

    try:
        leave_one_out_otu_df = leave_one_out_df.iloc[:, leave_one_out_df.columns.get_loc( 1 ):] # Get first OTU to last
    except: # if the column is not an int try using a string
        leave_one_out_otu_df = leave_one_out_df.iloc[:, leave_one_out_df.columns.get_loc( "1" ):] # Get first OTU to last
    leave_one_out_otu_df = leave_one_out_otu_df.T
    # get blank cells and replace them with NA
    leave_one_out_otu_df = leave_one_out_otu_df.astype('str')
    leave_one_out_otu_df = leave_one_out_otu_df.replace("nan", np.nan)

    # create data frame in order to plot as well as pass data easier
    kappa_df = pd.DataFrame(
        columns=['conditioning', 'centrality', 'correl', 'threshold', 'select_iter', 'kappa', 'agreement'])

    j = leave_one_out_otu_df.shape[1] -1 # index of highest iteration selection
    for i in range(0, j): # iterate through columns

        _kConditioning = leave_one_out_df["conditioning"].iloc[i]
        _kCentrality = leave_one_out_df["centrality"].iloc[i]
        _kCorrel = leave_one_out_df["correlation"].iloc[i]
        _kThreshold = leave_one_out_df["keep threshold"].iloc[i]
        _kSelect_iter = leave_one_out_df["iteration select"].iloc[i]
        _kKappa = None  # DEFINED LATER
        _kAgreement = None  # DEFINED LATER

        if (len( leave_one_out_otu_df.iloc[:, [i, j]].dropna()) >= 1):
            _kKappa = (rirr.kappa2((leave_one_out_otu_df.iloc[:, [i, j]]).dropna())[4][0])
            _kAgreement = jaccard_coefficient(leave_one_out_otu_df.iloc[:, i].dropna(), leave_one_out_otu_df.iloc[:, j].dropna())

        elif (len( leave_one_out_otu_df.iloc[:, [i, j]].dropna()) <= 0):
            _kKappa = np.nan
            _kAgreement = np.nan

        kappa_df.loc[i + 1] = [_kConditioning, _kCentrality, _kCorrel, _kThreshold, _kSelect_iter, _kKappa, _kAgreement]
    # It should be evident which result the kappa and agreement are being compared to
    kappa_df.loc[j+1] = [ leave_one_out_df["conditioning"].iloc[j],leave_one_out_df["centrality"].iloc[j],
                       leave_one_out_df["correlation"].iloc[j],leave_one_out_df["keep threshold"].iloc[j],
                       leave_one_out_df["iteration select"].iloc[j], 1.0, 1.0 ]

    kappa_df.reset_index( drop=True, inplace=True ) # index should start at 0 for consistency

    if( detailed ):
        out_file = f"{out_dir}/{name}_Jaccard_result.csv"
        # Create new files for output
        out_file = open(out_file, "w+", encoding="utf-8")
        kappa_df.to_csv( out_file )
        out_file.close()

    if( verbose ):
        # write all summaries to dump file
        dump = open( f"{out_dir}/step7_9_dump.txt", "w", encoding="utf-8" )
        _print_summaries( kappa_df, dump )
        dump.close()

    return kappa_df





