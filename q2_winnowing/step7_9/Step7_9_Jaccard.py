

# <><><> SETUP IMPORTS <><><>
import os
import pandas as pd
import numpy as np
from rpy2.robjects.packages import importr
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
from rpy2.robjects import r
from scipy.spatial.distance import pdist,squareform
from rpy2.robjects import pandas2ri
from rpy2.robjects import Formula
import rpy2.robjects as robjects
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

# Setup R packages used
pandas2ri.activate()

utils = rpackages.importr('utils')
utils.chooseCRANmirror(ind=1)

packnames = ('vegan', 'scales','data.table','zoo','dplyr','irr')
names_to_install = []
#names_to_install = [x for packnames if not rpackages.isinstalled(x)]

for x in packnames:
    if (rpackages.isinstalled(x)==False):
        names_to_install.append(x)

if len(names_to_install) > 0:
    utils.install_packages(StrVector(names_to_install))

rvegan = importr('vegan')
rscales = importr('scales')
rzoo = importr('zoo')
rdatatable = importr('data.table')
rstats = importr('stats')
rdplyr = importr('dplyr',on_conflict="warn")
rirr = importr('irr')
rsummary = r['summary']


# <><><> DEFINE FUNCTIONS <><><>

def _print_summaries( df_kappa, dump ):

    df_kappa['select_iter'] = df_kappa['select_iter'].astype('int')
    # <><><> OUTPUT TABLE <><><>
    dump.write( str( df_kappa ) + "\n\n")

    # <><><> PRINTS STATS OF COLUMNS <><><>
    dump.write(str( df_kappa.groupby(['conditioning', 'centrality', 'correl', 'threshold'])['agreement'].mean() ) + "\n")
    dump.write(str( df_kappa.loc[( df_kappa["conditioning"] == "Add1") & ( df_kappa["correl"] == "MIC")].groupby(
        ['centrality', 'correl', 'threshold'])['agreement'].mean() ) + "\n")

    return  # finished print



def jaccard_coefficient(x,y):
    return len(set(x) & set(y)) / len(set(x) | set(y)) # Set gets rid of duplicates


def main( df_leaveOneOut, name, detailed=False, verbose=False ):

    outDir = f"{os.path.dirname(os.path.realpath(__file__))}/output"
    # allows for cleaner execution and use of relative paths

    try:
        df_leaveOneOut_OTUs = df_leaveOneOut.iloc[:, df_leaveOneOut.columns.get_loc( 1 ):] # Get first OTU
    except: # if the column is not an int try using a string
        df_leaveOneOut_OTUs = df_leaveOneOut.iloc[:, df_leaveOneOut.columns.get_loc( "1" ):] # Get first OTU
    df_leaveOneOut_OTUs = df_leaveOneOut_OTUs.T
    # get blank cells and replace them with NA
    df_leaveOneOut_OTUs = df_leaveOneOut_OTUs.astype('str')
    df_leaveOneOut_OTUs = df_leaveOneOut_OTUs.replace("nan", np.nan)

    # create data frame in order to plot as well as pass data easier
    kappa_df = pd.DataFrame(
        columns=['conditioning', 'centrality', 'correl', 'threshold', 'select_iter', 'kappa', 'agreement'])

    j = df_leaveOneOut_OTUs.shape[1] -1 # index of highest iteration selection
    for i in range(0, j ): # iterate through columns

        _kConditioning = df_leaveOneOut["conditioning"].iloc[i]
        _kCentrality = df_leaveOneOut["centrality"].iloc[i]
        _kCorrel = df_leaveOneOut["correlation"].iloc[i]
        _kThreshold = df_leaveOneOut["keep threshold"].iloc[i]
        _kSelect_iter = df_leaveOneOut["iteration select"].iloc[i]
        _kKappa = None  # DEFINED LATER
        _kAgreement = None  # DEFINED LATER

        if (len( df_leaveOneOut_OTUs.iloc[:, [i, j]].dropna()) >= 1):
            _kKappa = (rirr.kappa2((df_leaveOneOut_OTUs.iloc[:, [i, j]]).dropna())[4][0])
            _kAgreement = jaccard_coefficient(df_leaveOneOut_OTUs.iloc[:, i].dropna(), df_leaveOneOut_OTUs.iloc[:, j].dropna())

        elif (len( df_leaveOneOut_OTUs.iloc[:, [i, j]].dropna()) <= 0):
            _kKappa = np.nan
            _kAgreement = np.nan

        kappa_df.loc[i + 1] = [_kConditioning, _kCentrality, _kCorrel, _kThreshold, _kSelect_iter, _kKappa, _kAgreement]


    if( detailed ):
        outFile = f"{outDir}/{name}_Jaccard_result.csv"
        # Create new files for output
        outFile = open(outFile, "w+", encoding="utf-8")
        kappa_df.to_csv( outFile )
        outFile.close()

    if( verbose ):
        # write all summaries to dump file
        dump = open( f"{outDir}/step7_9_dump.txt", "w", encoding="utf-8" )
        _print_summaries( kappa_df, dump )
        dump.close()

    return kappa_df




# <><><> TEST <><><>
# bfa = pd.read_csv("./test_data/arc_bac_fun_abundances.csv") # Hmmmmmmm
# bfa.rename(columns={'hor.plot': 'OTUid'}, inplace=True)  # Hmmmmmmm
#
# loo = pd.read_csv("./test_data/dataFromStep1.csv")
# main( loo, "test_results", True, True )
#
# loo = pd.read_csv("./test_data/firstValues.csv")
# main( loo, "test_firstValues", True, True )






















