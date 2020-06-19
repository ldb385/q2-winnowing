

# <><><> SETUP IMPORTS <><><>

import pandas as pd
import numpy as np
from rpy2.robjects.packages import importr
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
from rpy2.robjects import r
from scipy.spatial.distance import pdist,squareform
from rpy2.robjects import pandas2ri
from rpy2.robjects import IntVector, Formula
import rpy2.robjects as robjects
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
pandas2ri.activate()

# Setup R packages used

utils = rpackages.importr('utils')
utils.chooseCRANmirror(ind=1)

packnames = ('vegan', 'scales','data.table','zoo','DescTools','dplyr','irr')
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
rdesctools = importr('DescTools')
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


def _check_agreement( df_kappa, row, agreement_index ):

    if (df_kappa.iloc[row, agreement_index] == "" or df_kappa.iloc[row, agreement_index] == " "):
        df_kappa.iloc[row, agreement_index] = np.nan

    return # Nothing df is altered directly


def _check_kappa( df_kappa, row, kappa_index ):

    if (df_kappa.iloc[row, kappa_index] == 1 or df_kappa.iloc[row, kappa_index] == "" or df_kappa.iloc[row, kappa_index] == " "):
        df_kappa.iloc[row, kappa_index] = np.nan

    return  # Nothing df is altered directly


def _check_conditioning( df_kappa, row, conditioning_index, conditioning_flag ):

    if (conditioning_flag == True):
        df_kappa.iloc[ row, conditioning_index ] = "Hellinger"
    elif (conditioning_flag == False):
        df_kappa.iloc[ row, conditioning_index ] = "Add1"

    if (( row + 1) % 32 == 0):
        conditioning_flag = not conditioning_flag

    return conditioning_flag # new flag value


def _check_centrality( df_kappa, row, centrality_index, centrality_flag ):

    if (centrality_flag <= 4):
        df_kappa.iloc[ row, centrality_index ] = "bw"  # Betweeness
        centrality_flag += 1

    elif (centrality_flag > 4 and centrality_flag <= 8):
        df_kappa.iloc[ row, centrality_index ] = "cl"  # Closeness
        centrality_flag += 1

    elif (centrality_flag > 8 and centrality_flag <= 12):
        df_kappa.iloc[ row, centrality_index ] = "dg"  # Degree
        centrality_flag += 1

    elif (centrality_flag > 12 and centrality_flag <= 16):
        df_kappa.iloc[ row, centrality_index ] = "ev"  # Eigenvector
        centrality_flag += 1

        if (centrality_flag == 17):
            centrality_flag = 1

    return centrality_flag # new flag value


def _check_coorelation( df_kappa, row, correlation_index, correlation_flag):

    if ( correlation_flag == True):
        df_kappa.iloc[row, correlation_index] = "Spearman"

    elif (correlation_flag == False):
        df_kappa.iloc[row, correlation_index] = "MIC"

    if (( row + 1) % 16 == 0):
        correlation_flag = not correlation_flag

    return correlation_flag # new flag value


def _check_theshold( df_kappa, row, threshold_index, threshold_flag ):

    if (threshold_flag == 1):
        df_kappa.iloc[ row, threshold_index ] = 0.2

    elif (threshold_flag == 2):
        df_kappa.iloc[ row, threshold_index ] = 0.3

    elif (threshold_flag == 3):
        df_kappa.iloc[ row, threshold_index ] = 0.4

    if (( row + 1) % 64 == 0):
        threshold_flag += 1

    return threshold_flag # new flag value


def _check_iterations( df_kappa, row, iteration_index, iteration_flag ):

    if ( iteration_flag == 1):
        df_kappa.iloc[ row, iteration_index ] = 1
        iteration_flag += 1

    elif (iteration_flag == 2):
        df_kappa.iloc[ row, iteration_index ] = 4
        iteration_flag += 1

    elif (iteration_flag == 3):
        df_kappa.iloc[ row, iteration_index ] = 16
        iteration_flag += 1

    elif (iteration_flag == 4):
        df_kappa.iloc[ row, iteration_index ] = 64
        iteration_flag = 1

    return iteration_flag # new flag value


def main( df_leaveOneOut, name, detailed=False, verbose=False ):

    # TODO: TEST IF THIS PART IS NECCESSARY, MIGHT NOT BE
    df_leaveOneOut = df_leaveOneOut.iloc[:, 15:]
    df_leaveOneOut = df_leaveOneOut.T
    df_leaveOneOut = df_leaveOneOut.astype('str')
    df_leaveOneOut = df_leaveOneOut.replace("nan", np.nan)
    # ^^^^^ THIS PART ^^^^^

    j = 4
    # create data frame in order to plot as well as pass data easier
    kappa_df = pd.DataFrame(
        columns=['conditioning', 'centrality', 'correl', 'threshold', 'select_iter', 'kappa', 'agreement'])

    print( df_leaveOneOut.shape )

    for i in range(0, df_leaveOneOut.shape[1] - 1): # iterate through columns
        if (i == j):
            j += 5
            i += 1

        _kConditioning = "0"
        _kCentrality = "0"
        _kCorrel = "0"
        _kThreshold = 0.0
        _kSelect_iter = 0
        _kKappa = None # DEFINED LATER
        _kAgreement = None # DEFINED LATER

        if (len( df_leaveOneOut.iloc[:, [i, j]].dropna()) >= 1):
            _kKappa = (rirr.kappa2((loo.iloc[:, [i, j]]).dropna())[4][0])
            _kAgreement = jaccard_coefficient(loo.iloc[:, i].dropna(), loo.iloc[:, j].dropna())

        elif (len( df_leaveOneOut.iloc[:, [i, j]].dropna()) <= 0):
            _kKappa = np.nan
            _kAgreement = np.nan

        kappa_df.loc[i + 1] = [_kConditioning, _kCentrality, _kCorrel, _kThreshold, _kSelect_iter, _kKappa, _kAgreement]


    # Flags will be used in order to make sure proper results are output
    cond_flag = True # Conditioning
    cent_flag = 1 # Centrality
    corre_flag = True # Correlation
    thres_flag = 1 # Threshold
    iter_flag = 1 # Iterations Selected

    for i in range(0, len( kappa_df )): # Fill Kappa df with data

        _check_agreement( kappa_df, i, 6)
        _check_kappa( kappa_df, i, 5 )

        cond_flag = _check_conditioning( kappa_df, i, 0, cond_flag )
        cent_flag = _check_centrality( kappa_df, i, 1, cent_flag )
        corre_flag = _check_coorelation( kappa_df, i, 2, corre_flag )
        thres_flag = _check_theshold( kappa_df, i, 3, thres_flag )
        iter_flag = _check_iterations( kappa_df, i, 4, iter_flag )

    if( detailed ):
        outFile = f"output/{name}_Jaccard_result.csv"
        # Create new files for output
        outFile = open(outFile, "w+", encoding="utf-8")
        kappa_df.to_csv( outFile )
        outFile.close()

    if( verbose ):
        # write all summaries to dump file
        dump = open("output/step7_9_dump.txt", "w", encoding="utf-8")
        _print_summaries( kappa_df, dump )
        dump.close()

    return kappa_df




# <><><> TEST <><><>
bfa = pd.read_csv("./test_data/arc_bac_fun_abundances.csv") # Hmmmmmmm
bfa.rename(columns={'hor.plot': 'OTUid'}, inplace=True)  # Hmmmmmmm

loo = pd.read_csv("./test_data/results-brome_total_bad_removed.csv")

main( loo, "test_results", True, True )






















