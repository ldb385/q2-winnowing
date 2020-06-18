
# <><><> SETUP IMPORTS <><><>

import os

import pandas as pd
import numpy as np
from rpy2.robjects.packages import importr
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector, FloatVector
from rpy2.robjects import r
from rpy2.robjects import default_converter, conversion
from scipy.spatial.distance import pdist,squareform
from rpy2.robjects import pandas2ri
from rpy2.robjects import IntVector, Formula
import rpy2.robjects as robjects
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from rpy2.robjects.conversion import localconverter

# Install packages
pandas2ri.activate()

utils = rpackages.importr('utils')
utils.chooseCRANmirror(ind=1)

packnames = ('vegan', 'scales','data.table','zoo')
names_to_install = []
#names_to_install = [x for packnames if not rpackages.isinstalled(x)]

for x in packnames:
    if (rpackages.isinstalled(x)==False):
        names_to_install.append(x)

if len(names_to_install) > 0:
    utils.install_packages(StrVector(names_to_install))

# Setup R packages
rvegan = importr('vegan')
rscales = importr('scales')
rzoo = importr('zoo')
rdatatable = importr('data.table')
rstats = importr('stats')
rprint = r['print']
rpar = r['par']
rplot = r['plot']


# <><><> DEFINE FUNCTIONS <><><>

def _generate_figures( dataFrame_permanova, outdir ):

    rf_model_scale = pandas2ri.py2rpy(dataFrame_permanova["F.model.scale"])
    sliding_sd = rzoo.rollapply(rf_model_scale, width=5, FUN=rstats.sd, fill='NA')

    # <><> PLOTTING DEGREE <><>
    plt.figure()
    plt.axis([0, 100, 0.0, 1.0])
    plt.xlabel('AUC%')
    plt.ylabel('F-Score(scaled)')
    plt.plot(
        dataFrame_permanova.iloc[[x for x in range(39, 45)] + [y for y in range(45, 101, 5)], [2]].values.tolist(),
        dataFrame_permanova.iloc[[x for x in range(39, 45)] + [y for y in range(45, 101, 5)], [9]].values.tolist(),
        'ro', color="black")
    plt.plot(dataFrame_permanova.iloc[[x for x in range(0, 36, 5)] + [37, 38, 39], [2]].values.tolist(),
             dataFrame_permanova.iloc[[x for x in range(0, 36, 5)] + [36, 37, 38], [9]].values.tolist(), 'ro',
             color="red")
    spl = UnivariateSpline(dataFrame_permanova["auc"].values.tolist(),
                           dataFrame_permanova["F.model.scale"].values.tolist())
    xs = np.linspace(0, 100, 10)
    plt.plot(xs, spl(xs), 'g', lw=1)
    # plt.plot(brome_permanova["auc"].values.tolist(),brome_permanova["F.model.scale"],color="green")
    plt.plot([x for x in range(0, 101)], sliding_sd, color="red")
    plt.plot(dataFrame_permanova.iloc[36, 2], dataFrame_permanova.iloc[36, 9], 'ro', color="blue")
    plt.text(dataFrame_permanova.iloc[36, 2] + 40, dataFrame_permanova.iloc[36, 9],
             str(dataFrame_permanova.iloc[36, 8]) + " OTUs", color="blue")
    plt.legend(loc=1, labels="E")
    plt.savefig( os.path.join( outdir,"F-Score_AUC_degree.png"))

    # <><> PLOTTING CLOSENESS <><>
    plt.figure()
    plt.axis([0, 100, 0.0, 1.1])
    plt.xlabel('AUC%')
    plt.ylabel('F-Score(scaled)')
    plt.plot(dataFrame_permanova.iloc[[8, 9, 10] + [x for x in range(11, 101, 5)], [2]].values.tolist(),
             dataFrame_permanova.iloc[[8, 9, 10] + [x for x in range(11, 101, 5)], [9]].values.tolist(), 'ro',
             color="black")
    plt.plot(dataFrame_permanova.iloc[[x for x in range(1, 7)], [2]].values.tolist(),
             dataFrame_permanova.iloc[[x for x in range(1, 7)], [9]].values.tolist(), 'ro', color="red")
    spl = UnivariateSpline(dataFrame_permanova["auc"].values.tolist(),
                           dataFrame_permanova["F.model.scale"].values.tolist())
    xs = np.linspace(0, 101, 10)
    plt.plot(xs, spl(xs), 'g', lw=1)
    plt.plot([x for x in range(0, 101)], sliding_sd, color="red")
    plt.plot(dataFrame_permanova.iloc[4, 2], dataFrame_permanova.iloc[4, 9], 'ro', color="blue")
    plt.text(dataFrame_permanova.iloc[4, 2] + 40, dataFrame_permanova.iloc[4, 9],
             str(dataFrame_permanova.iloc[4, 8]) + " OTUs", color="blue")
    plt.legend(loc=1, labels="F")
    plt.savefig( os.path.join( outdir,"F-Score_AUC_closeness.png"))

    # <><> PLOTTING BETWEENESS <><>
    plt.figure()
    plt.axis([0, 100, 0.0, 1.0])
    plt.xlabel('AUC%')
    plt.ylabel('F-Score(scaled)')
    plt.plot(
        dataFrame_permanova.iloc[[x for x in range(36, 40)] + [y for y in range(40, 101, 5)], [2]].values.tolist(),
        dataFrame_permanova.iloc[[x for x in range(36, 40)] + [y for y in range(40, 101, 5)], [9]].values.tolist(),
        'ro', color="black")
    plt.plot(dataFrame_permanova.iloc[[x for x in range(0, 31, 5)] + [32, 33, 34, 35], [2]].values.tolist(),
             dataFrame_permanova.iloc[[x for x in range(0, 31, 5)] + [32, 33, 34, 35], [9]].values.tolist(), 'ro',
             color="red")
    spl = UnivariateSpline(dataFrame_permanova["auc"].values.tolist(),
                           dataFrame_permanova["F.model.scale"].values.tolist())
    xs = np.linspace(0, 100, 10)
    plt.plot(xs, spl(xs), 'g', lw=1)
    plt.plot([x for x in range(0, 101)], sliding_sd, color="red")
    plt.plot(dataFrame_permanova.iloc[35, 2], dataFrame_permanova.iloc[35, 9], 'ro', color="blue")
    plt.text(dataFrame_permanova.iloc[35, 2] + 40, dataFrame_permanova.iloc[35, 9] + 0.05,
             str(dataFrame_permanova.iloc[35, 8]) + " OTUs", color="blue")
    plt.legend(loc=1, labels="G")
    plt.savefig( os.path.join( outdir,"F-Score_AUC_betweeness.png") )


    # <><> PLOTTING EIGENVECTOR <><>
    plt.figure()
    plt.axis([0, 100, 0.0, 1.0])
    plt.xlabel('AUC%')
    plt.ylabel('F-Score(scaled)')
    plt.plot(
        dataFrame_permanova.iloc[[x for x in range(26, 31)] + [y for y in range(30, 101, 5)], [2]].values.tolist(),
        dataFrame_permanova.iloc[[x for x in range(26, 31)] + [y for y in range(30, 101, 5)], [9]].values.tolist(),
        'ro', color="black")
    plt.plot(dataFrame_permanova.iloc[[x for x in range(0, 21, 5)] + [22, 23], [2]].values.tolist(),
             dataFrame_permanova.iloc[[x for x in range(0, 21, 5)] + [22, 23], [9]].values.tolist(), 'ro',
             color="red")
    spl = UnivariateSpline(dataFrame_permanova["auc"].values.tolist(),
                           dataFrame_permanova["F.model.scale"].values.tolist())
    xs = np.linspace(0, 100, 5)
    plt.plot(xs, spl(xs), 'g', lw=1)
    plt.plot([x for x in range(0, 101)], sliding_sd, color="red")
    plt.plot(dataFrame_permanova.iloc[21, 2], dataFrame_permanova.iloc[21, 9], 'ro', color="blue")
    plt.text(dataFrame_permanova.iloc[21, 2] + 40, dataFrame_permanova.iloc[21, 9],
             str(dataFrame_permanova.iloc[21, 8]) + " OTUs", color="blue")
    plt.legend(loc=1, labels="H")
    plt.savefig( os.path.join( outdir,"F-Score_AUC_eigenvector.png") )

    return #Nothing is returned this states end of function


def perform_permanova_dataFrame( sample_file, data_frame_1, data_frame_2, output_file, detailed=False, verbose=False, dump=None):

    if( verbose ):
        dump.write( f"Processing Input File: {str(data_frame_1)} \n")
        dump.write( f"Processing Input file: {str(data_frame_2)} \n")

    df_sample = pd.read_csv(sample_file) # AB horizon, BRA

    df_dg_auc = data_frame_1 # AUC curves for subsetting
    df_dg_auc100 = data_frame_2 # Abundances used for Hellinger

    # This is what our dataframe output should be formatted as, will need to build as each index is iterated
    df_dg_premanova = pd.DataFrame(columns=['test', 'order', 'auc','SumsOfSqs','MeanSqs','F.model','R2','Pval','N.taxa','F.model.scale'])

    # Setup variables before loop to avoid extra computation
    typeIndex = int( df_sample.columns.get_loc("type") ) # Since Rpy doesn't allow for string index must get column int index
    with localconverter( default_converter + pandas2ri.converter ):
        rdf_sample = conversion.py2rpy( df_sample )
    rdf_data_sample = rdf_sample[ typeIndex ]

    rformula = Formula('x ~ y') # Note, y is independent while x is dependant

    for i in range( 0, len( df_dg_auc) ): # all auc values are used to compute

        # This is STEP 1
        df_data_dg = df_dg_auc100.iloc[:, 0:int(df_dg_auc.iloc[i,2]) ] # get all columns with rows 0 to value at [i,2]
        with localconverter(default_converter + pandas2ri.converter):
            rdf_data_dg = conversion.py2rpy(df_data_dg)

        # This is STEP 2
        # Convert to Hellinger distance matrix
        rdf_data_dg_hel = rvegan.vegdist(rvegan.decostand( rdf_data_dg, "hellinger"), "euclidean")


        # This is STEP 3
        # Setup formula for the adonis calculation
        rdf_data_dg_hel = rdf_data_dg_hel.reshape(( len(rdf_data_sample),-1)) # -1 infers columns based on rows
        renv = rformula.environment
        renv['x'] = rdf_data_dg_hel # x should have equal amount of rows as y has values
        renv['y'] = rdf_data_sample

        radonis = rvegan.adonis( rformula, permutations=999)

        # This is STEP 4
        # Must build the permanova table stated earlier
        _pTest = f"auc{str(i+1)}"
        _pOrder = f"{i+1}"
        _pAUC = i+1
        _pSumOfSquares = radonis[0].iloc[0, 1]
        _pMeanSquares = radonis[0].iloc[0, 2]
        _pFModel = radonis[0].iloc[0, 3]
        _pR2 = radonis[0].iloc[0, 4]
        _pPVal = radonis[0].iloc[0, 5]
        _pNTaxa = len(df_data_dg.columns)
        _pFModelScale = -1.0 # this just a place holder

        df_dg_premanova.loc[i] = [_pTest] + [_pOrder] + [_pAUC] + [_pSumOfSquares] + [_pMeanSquares] + \
                                 [_pFModel] + [_pR2] + [_pPVal] + [_pNTaxa] + [_pFModelScale]

    # This is STEP 4.5
    # Convert F.model to a scaled version of F.model
    df_dg_premanova["F.model.scale"] = (df_dg_premanova["F.model"] - min(df_dg_premanova["F.model"]))/(max(df_dg_premanova["F.model"]) - min(df_dg_premanova["F.model"]))
    df_dg_premanova.index = range(1,len(df_dg_premanova)+1)

    if( detailed ):
        df_dg_premanova.to_csv( output_file )

    # df_dg_premanova.loc[0] = ["auc0"]+[0]+[0]+[0]+[0]+[0]+[0]+[0]+[0]+[0]
    # df_dg_premanova.sort_index(axis=0,inplace=True)

    if( verbose and detailed ):
        dump.write( f"Output is Written in File: {str(output_file)} \n\n" )

    return df_dg_premanova



# <><><> DEFINE EXECUTION FUNCTION <><><>
def main_dataFrame( dataFrame1, dataFrame2, sampleFile,  name, detailed=False, verbose=False):

    if( detailed ):
        outFile = f"output/{name}_PERMANOVA_result.csv"
        # Create new files for output
        outFile = open(outFile, "w+", encoding="utf-8")
        outDir = "./output/"

        if( verbose ):

            dump = open("output/step6_dump.txt", "w", encoding="utf-8")

            # Call PERMANOVA calculation
            df_permanova = perform_permanova_dataFrame( sampleFile, dataFrame1, dataFrame2, outFile, detailed, verbose, dump )
            dump.write( f"Plots generated to {outDir}.\n" )
            dump.close()

        else:
            # Call PERMANOVA calculation
            df_permanova = perform_permanova_dataFrame( sampleFile, dataFrame1, dataFrame2, outFile, detailed )

        # Since this is detailed must generate plots
        _generate_figures( df_permanova, outDir )

    elif( verbose ):

        dump = open("output/step6_dump.txt", "w", encoding="utf-8")

        # Call PERMANOVA calculation
        df_permanova = perform_permanova_dataFrame( sampleFile, dataFrame1, dataFrame2, verbose=verbose, dump=dump )
        dump.close()

    else:
        # No excess files necessary just generate dataframe to pass on
        df_permanova = perform_permanova_dataFrame( sampleFile, dataFrame1, dataFrame2 )

    return df_permanova


# <><> TEST <><>
# Testing dataframe function
test_sample = "./test_data/Brome_BFA_AB_sample_info.csv"
test_abundance = pd.read_csv("./test_data/ADD1_AUC100_MIC0.2_Brome_bacfunarc_dw_otu_table-graph_centrality-degree-selectallbyall-abundances.csv")
test_AUC = pd.read_csv("./test_data/brome.dg.auc.csv")
main_dataFrame( test_AUC, test_abundance, test_sample, "Test_Step6", True, True)

# test_sample = "./test_data/test_sample.csv"
# test_abundance = pd.read_csv("./test_data/bromeA_all-abundances-0.csv")
# test_AUC = pd.read_csv("./test_data/testframe_1_auc_result.csv")
# main_dataFrame( test_abundance, test_AUC, test_sample, "Test_Step6", True, True)











