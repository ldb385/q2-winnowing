

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

    # <><><> PRINTS STATS OF COLUMNS <><><>
    dump.write(str( df_kappa.groupby(['conditioning', 'centrality', 'correl', 'threshold'])['agreement'].mean() ) + "\n")
    dump.write(str( df_kappa.loc[( df_kappa["conditioning"] == "Add1") & ( df_kappa["correl"] == "MIC")].groupby(
        ['centrality', 'correl', 'threshold'])['agreement'].mean() ) + "\n")
    raov = rstats.aov(Formula('agreement ~ conditioning + centrality + correl + as.factor(threshold)'), df_kappa)
    dump.write(str( rsummary(raov) ) + "\n")

    # <><><> PRINTS STATS OF COLUMNS <><><>
    dump.write( "Conditioning" + "\n")
    dump.write(str( rstats.TukeyHSD(raov)[0] ) + "\n")
    dump.write( "centrality" + "\n")
    dump.write(str( rstats.TukeyHSD(raov)[1] ) + "\n")
    dump.write( "correl" + "\n")
    dump.write(str( rstats.TukeyHSD(raov)[2] ) + "\n")
    dump.write( "as.factor(threshold)" + "\n")
    dump.write(str( rstats.TukeyHSD(raov)[1] ) + "\n")

    # <><><> CREATE TEMP DFs FOR EASY OUTPUT SUMMARY OF DIFFERENT CONDITIONINGS <><><>
    rformula = Formula('agreement ~ centrality + correl + as.factor(threshold)')
    dump.write(str( rsummary(rstats.aov(rformula, df_kappa.loc[(df_kappa["conditioning"] == "Add1")])) ) + "\n")
    dump.write(str( rsummary(rstats.aov(rformula, df_kappa.loc[(df_kappa["conditioning"] == "Hellinger")])) ) + "\n")

    # <><><> SEARCH THROUGH DF FOR DIFFERENT MIC VALUES <><><>

    # MIC > 0.2
    rformula_1 = Formula('agreement ~ centrality + select_iter')
    dump.write(str( rsummary(rstats.aov(rformula_1, df_kappa.loc[
        ((df_kappa["conditioning"] == "Add1") & (df_kappa["correl"] == "MIC") & (df_kappa["threshold"] == 0.2))])) ) + "\n")
    dump.write(str( rsummary(rstats.aov(rformula_1, df_kappa.loc[
        ((df_kappa["conditioning"] == "Hellinger") & (df_kappa["correl"] == "MIC") & (df_kappa["threshold"] == 0.2))])) ) + "\n")

    rformula = Formula('agreement ~ centrality*select_iter')
    dump.write(str( rsummary(rstats.aov(rformula, df_kappa.loc[
        ((df_kappa["conditioning"] == "Add1") & (df_kappa["correl"] == "Spearman") & (df_kappa["threshold"] == 0.2))])) ) + "\n")
    dump.write(str( rsummary(rstats.aov(rformula, df_kappa.loc[(
                (df_kappa["conditioning"] == "Hellinger") & (df_kappa["correl"] == "Spearman") & (
                    df_kappa["threshold"] == 0.2))])) ) + "\n")

    # MIC > 0.3
    dump.write(str( rsummary(rstats.aov(rformula, df_kappa.loc[
        ((df_kappa["conditioning"] == "Add1") & (df_kappa["correl"] == "MIC") & (df_kappa["threshold"] == 0.3))])) ) + "\n")
    dump.write(str( rsummary(rstats.aov(rformula_1, df_kappa.loc[
        ((df_kappa["conditioning"] == "Hellinger") & (df_kappa["correl"] == "MIC") & (df_kappa["threshold"] == 0.3))])) ) + "\n")
    dump.write(str( rsummary(rstats.aov(rformula, df_kappa.loc[
        ((df_kappa["conditioning"] == "Add1") & (df_kappa["correl"] == "Spearman") & (df_kappa["threshold"] == 0.3))])) ) + "\n")
    dump.write(str( rsummary(rstats.aov(rformula, df_kappa.loc[(
                (df_kappa["conditioning"] == "Hellinger") & (df_kappa["correl"] == "Spearman") & (
                    df_kappa["threshold"] == 0.3))])) ) + "\n")

    # MIC > 0.4
    dump.write(str( rsummary(rstats.aov(rformula, df_kappa.loc[
        ((df_kappa["conditioning"] == "Add1") & (df_kappa["correl"] == "MIC") & (df_kappa["threshold"] == 0.4))])) ) + "\n")
    dump.write(str( rsummary(rstats.aov(rformula, df_kappa.loc[
        ((df_kappa["conditioning"] == "Hellinger") & (df_kappa["correl"] == "MIC") & (df_kappa["threshold"] == 0.4))])) ) + "\n")
    dump.write(str( rsummary(rstats.aov(rformula, df_kappa.loc[
        ((df_kappa["conditioning"] == "Add1") & (df_kappa["correl"] == "Spearman") & (df_kappa["threshold"] == 0.4))])) ) + "\n")
    dump.write(str( rsummary(rstats.aov(rformula, df_kappa.loc[(
                (df_kappa["conditioning"] == "Hellinger") & (df_kappa["correl"] == "Spearman") & (
                    df_kappa["threshold"] == 0.4))])) ) + "\n")

    return # Nothing returned this just write to a dump file



def jaccard_coefficient(x,y):
    return len(set(x) & set(y)) / len(set(x) | set(y)) # Set gets rid of duplicates


def main( df_leaveOneOut, name, detailed=False, verbose=False ):

    outDir = f"{os.path.dirname(os.path.realpath(__file__))}/output"
    # allows for cleaner execution and use of relative paths

    print( df_leaveOneOut, "\n\n\n")
    try:
        df_leaveOneOut_OTUs = df_leaveOneOut.iloc[:, df_leaveOneOut.columns.get_loc( 1 ):] # Get first OTU
    except: # if the column is not an int try using a string
        df_leaveOneOut_OTUs = df_leaveOneOut.iloc[:, df_leaveOneOut.columns.get_loc( "1" ):] # Get first OTU
    df_leaveOneOut_OTUs = df_leaveOneOut_OTUs.T
    # get blank cells and replace them with NA
    df_leaveOneOut_OTUs = df_leaveOneOut_OTUs.astype('str')
    df_leaveOneOut_OTUs = df_leaveOneOut_OTUs.replace("nan", np.nan)
    print( df_leaveOneOut_OTUs, "\n\n\n")

    # create data frame in order to plot as well as pass data easier
    kappa_df = pd.DataFrame(
        columns=['conditioning', 'centrality', 'correl', 'threshold', 'select_iter', 'kappa', 'agreement'])

    _kConditioning = df_leaveOneOut["conditioning"].iloc[0]
    _kCentrality = df_leaveOneOut["centrality"].iloc[0]
    _kCorrel = df_leaveOneOut["correlation"].iloc[0]
    _kThreshold = df_leaveOneOut["keep threshold"].iloc[0]
    _kSelect_iter = df_leaveOneOut["iteration select"].iloc[0]
    _kKappa = None # DEFINED LATER
    _kAgreement = None # DEFINED LATER


    if (len( df_leaveOneOut_OTUs.iloc[:, 0].dropna()) >= 1): # all rows, columns i and j
        print( rirr.kappa2( df_leaveOneOut_OTUs.iloc[:, 0] ) )
        _kKappa = ( rirr.kappa2( df_leaveOneOut_OTUs.iloc[:, 0] )[4][0] )
        _kAgreement = jaccard_coefficient(df_leaveOneOut_OTUs.iloc[:, 0].dropna(), df_leaveOneOut_OTUs.iloc[:, :].dropna())

    elif (len( df_leaveOneOut_OTUs.iloc[:, :].dropna()) <= 0):
        _kKappa = np.nan
        _kAgreement = np.nan

    # input kappa and agreement values
    kappa_df.loc[1] = [_kConditioning, _kCentrality, _kCorrel, _kThreshold, _kSelect_iter, _kKappa, _kAgreement]

    if( detailed ):
        outFile = f"{outDir}/{name}_Jaccard_result.csv"
        # Create new files for output
        outFile = open(outFile, "w+", encoding="utf-8")
        kappa_df.to_csv( outFile )
        outFile.close()

    if( verbose ):
        # write all summaries to dump file
        dump = open( f"{outDir}/step7_9_dump.txt", "w", encoding="utf-8" )
        # _print_summaries( kappa_df, dump )
        dump.close()

    return kappa_df




# <><><> TEST <><><>
# bfa = pd.read_csv("./test_data/arc_bac_fun_abundances.csv") # Hmmmmmmm
# bfa.rename(columns={'hor.plot': 'OTUid'}, inplace=True)  # Hmmmmmmm
#
loo = pd.read_csv("./test_data/dataFromStep1.csv")

main( loo, "test_results", True, True )






















