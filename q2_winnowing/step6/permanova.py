
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

def _generate_figures( permanova_df, cent_type, outdir, name ):
    """
    Generate figures depending on centrality type was used
    :param permanova_df: PERMANOVA table generated
    :param cent_type: centrality type used ( None, degree, closeness, eigenvector, betweenness )
    :param outdir: directory where figures are written to
    :param name: name attached to figures
    :return:
    """

    if( cent_type is None ):
        # cent type was not given don't attempt to graph
        return # Nothing simply termination of function

    permanova_df_copy = permanova_df.copy( deep=True ) # deep copy to avoid altering original dataframe

    permanova_df_copy.loc[0] = ["auc0"]+[0]+[0]+[0]+[0]+[0]+[0]+[0]+[0]+[0]
    permanova_df_copy.sort_index(axis=0,inplace=True)

    rf_model_scale = pandas2ri.py2rpy(permanova_df_copy["F.model.scale"])
    sliding_sd = rzoo.rollapply(rf_model_scale, width=5, FUN=rstats.sd, fill='NA')

    if( cent_type == "degree" ):
        # <><> PLOTTING DEGREE <><>
        plt.figure()
        plt.axis([0, 100, 0.0, 1.0])
        plt.xlabel('AUC%')
        plt.ylabel('F-Score(scaled)')
        plt.plot(
            permanova_df_copy.iloc[[x for x in range(39, 45)] + [y for y in range(45, 101, 5)], [2]].values.tolist(),
            permanova_df_copy.iloc[[x for x in range(39, 45)] + [y for y in range(45, 101, 5)], [9]].values.tolist(),
            'ro', color="black")
        plt.plot(permanova_df_copy.iloc[[x for x in range(0, 36, 5)] + [37, 38, 39], [2]].values.tolist(),
                 permanova_df_copy.iloc[[x for x in range(0, 36, 5)] + [36, 37, 38], [9]].values.tolist(), 'ro',
                 color="red")
        spl = UnivariateSpline(permanova_df_copy["auc"].values.tolist(),
                               permanova_df_copy["F.model.scale"].values.tolist())
        xs = np.linspace(0, 100, 10)
        plt.plot(xs, spl(xs), 'g', lw=1)
        # plt.plot(brome_permanova["auc"].values.tolist(),brome_permanova["F.model.scale"],color="green")
        plt.plot([x for x in range(0, 101)], sliding_sd, color="red")
        plt.plot(permanova_df_copy.iloc[36, 2], permanova_df_copy.iloc[36, 9], 'ro', color="blue")
        plt.text(permanova_df_copy.iloc[36, 2] + 40, permanova_df_copy.iloc[36, 9],
                 str(permanova_df_copy.iloc[36, 8]) + " OTUs", color="blue")
        plt.legend(loc=1, labels="E")
        plt.savefig( os.path.join( outdir, f"{name}_F-Score_AUC_degree.png"))

    elif( cent_type == "closeness" ):
        # <><> PLOTTING CLOSENESS <><>
        plt.figure()
        plt.axis([0, 100, 0.0, 1.1])
        plt.xlabel('AUC%')
        plt.ylabel('F-Score(scaled)')
        plt.plot(permanova_df_copy.iloc[[8, 9, 10] + [x for x in range(11, 101, 5)], [2]].values.tolist(),
                 permanova_df_copy.iloc[[8, 9, 10] + [x for x in range(11, 101, 5)], [9]].values.tolist(), 'ro',
                 color="black")
        plt.plot(permanova_df_copy.iloc[[x for x in range(1, 7)], [2]].values.tolist(),
                 permanova_df_copy.iloc[[x for x in range(1, 7)], [9]].values.tolist(), 'ro', color="red")
        spl = UnivariateSpline(permanova_df_copy["auc"].values.tolist(),
                               permanova_df_copy["F.model.scale"].values.tolist())
        xs = np.linspace(0, 101, 10)
        plt.plot(xs, spl(xs), 'g', lw=1)
        plt.plot([x for x in range(0, 101)], sliding_sd, color="red")
        plt.plot(permanova_df_copy.iloc[4, 2], permanova_df_copy.iloc[4, 9], 'ro', color="blue")
        plt.text(permanova_df_copy.iloc[4, 2] + 40, permanova_df_copy.iloc[4, 9],
                 str(permanova_df_copy.iloc[4, 8]) + " OTUs", color="blue")
        plt.legend(loc=1, labels="F")
        plt.savefig( os.path.join( outdir, f"{name}_F-Score_AUC_closeness.png"))

    elif( cent_type == "betweenness" ):
        # <><> PLOTTING BETWEENESS <><>
        plt.figure()
        plt.axis([0, 100, 0.0, 1.0])
        plt.xlabel('AUC%')
        plt.ylabel('F-Score(scaled)')
        plt.plot(
            permanova_df_copy.iloc[[x for x in range(36, 40)] + [y for y in range(40, 101, 5)], [2]].values.tolist(),
            permanova_df_copy.iloc[[x for x in range(36, 40)] + [y for y in range(40, 101, 5)], [9]].values.tolist(),
            'ro', color="black")
        plt.plot(permanova_df_copy.iloc[[x for x in range(0, 31, 5)] + [32, 33, 34, 35], [2]].values.tolist(),
                 permanova_df_copy.iloc[[x for x in range(0, 31, 5)] + [32, 33, 34, 35], [9]].values.tolist(), 'ro',
                 color="red")
        spl = UnivariateSpline(permanova_df_copy["auc"].values.tolist(),
                               permanova_df_copy["F.model.scale"].values.tolist())
        xs = np.linspace(0, 100, 10)
        plt.plot(xs, spl(xs), 'g', lw=1)
        plt.plot([x for x in range(0, 101)], sliding_sd, color="red")
        plt.plot(permanova_df_copy.iloc[35, 2], permanova_df_copy.iloc[35, 9], 'ro', color="blue")
        plt.text(permanova_df_copy.iloc[35, 2] + 40, permanova_df_copy.iloc[35, 9] + 0.05,
                 str(permanova_df_copy.iloc[35, 8]) + " OTUs", color="blue")
        plt.legend(loc=1, labels="G")
        plt.savefig( os.path.join( outdir, f"{name}_F-Score_AUC_betweeness.png") )

    elif( cent_type == "eigenvector"):
        # <><> PLOTTING EIGENVECTOR <><>
        plt.figure()
        plt.axis([0, 100, 0.0, 1.0])
        plt.xlabel('AUC%')
        plt.ylabel('F-Score(scaled)')
        plt.plot(
            permanova_df_copy.iloc[[x for x in range(26, 31)] + [y for y in range(30, 101, 5)], [2]].values.tolist(),
            permanova_df_copy.iloc[[x for x in range(26, 31)] + [y for y in range(30, 101, 5)], [9]].values.tolist(),
            'ro', color="black")
        plt.plot(permanova_df_copy.iloc[[x for x in range(0, 21, 5)] + [22, 23], [2]].values.tolist(),
                 permanova_df_copy.iloc[[x for x in range(0, 21, 5)] + [22, 23], [9]].values.tolist(), 'ro',
                 color="red")
        spl = UnivariateSpline(permanova_df_copy["auc"].values.tolist(),
                               permanova_df_copy["F.model.scale"].values.tolist())
        xs = np.linspace(0, 100, 5)
        plt.plot(xs, spl(xs), 'g', lw=1)
        plt.plot([x for x in range(0, 101)], sliding_sd, color="red")
        plt.plot(permanova_df_copy.iloc[21, 2], permanova_df_copy.iloc[21, 9], 'ro', color="blue")
        plt.text(permanova_df_copy.iloc[21, 2] + 40, permanova_df_copy.iloc[21, 9],
                 str(permanova_df_copy.iloc[21, 8]) + " OTUs", color="blue")
        plt.legend(loc=1, labels="H")
        plt.savefig( os.path.join( outdir, f"{name}_F-Score_AUC_eigenvector.png") )


    return #Nothing is returned this states end of function



def _convert_to_dist_hel_matrix( array, length ):
    """
    Convert array to a distance hellinger matrix
    :param array: the array with all values
    :param length: the length x width of what the new matrix will be
    :return: output distance hellinger matrix
    """

    # initialize an empty matrix to fill
    matrix_new = np.zeros( (length, length), dtype=float )

    print( len( array ) , len( matrix_new ) )

    col_index = 0
    row_index = 1

    # needs to fill by columns instead of by rows so need to do this manually
    for value in array:
        matrix_new[row_index][col_index] = value

        if( row_index < length -1 ):
            row_index += 1
        else:
            col_index += 1
            row_index = col_index + 1

    # it needs the upper triangle to mirror the lower
    matrix_new[ np.triu_indices(length, 1 )] = array

    return matrix_new



def perform_permanova( auc_df, auc100_df, sample_df, out_file, detailed=False, verbose=False, dump=None):
    """
    Perform ANOVA calculation with permutations
    :param sample_df: Sample that provide data on whether taxom were invaded or natural
    :param auc_df: AUC ordering
    :param auc100_df: abundances generated
    :param out_file: this is the result folder where output is written to if detailed
    :param detailed: Output helper tables
    :param verbose: Output helper prints
    :param dump: file were verbose is being written to
    :return: data frame containing calculations performed
    """

    if( verbose ):
        dump.write( f"Processing Input File: {str(auc_df)} \n") # AUC curves for subsetting
        dump.write( f"Processing Input file: {str(auc100_df)} \n") # Abundances used for Hellinger

    # Make sure indices begin at 0 so loop works
    auc_df.reset_index( drop=True, inplace=True )
    auc100_df.reset_index(  drop=True, inplace=True )

    # This is what our dataframe output should be formatted as, will need to build as each index is iterated
    premanova_df = pd.DataFrame(columns=['test', 'order', 'auc','SumsOfSqs','MeanSqs','F.model','R2','Pval','N.taxa','F.model.scale'])

    # Setup variables before loop to avoid extra computation
    if ("type" in sample_df.columns): # check for capitalization in index
        type_index = int( sample_df.columns.get_loc("type"))  # Since Rpy doesn't allow for string index must get column int index
    else:
        type_index = int( sample_df.columns.get_loc("Type"))  # Since Rpy doesn't allow for string index must get column int index
    sample_rdf = pandas2ri.py2rpy(sample_df)[type_index]

    rformula = Formula('x ~ y') # Note, y is independent while x is dependant

    for i in range( 0, len( auc_df) ): # all auc values are used to compute

        # This is STEP 1
        range_of_decay = int(auc_df.iloc[i, len(auc_df.columns)-1])
        data_df = auc100_df.iloc[:, 0:range_of_decay ] # get all rows with columns 0 to value at [i,2]
        data_rdf = pandas2ri.py2rpy(data_df)

        # This is STEP 2
        # Convert to Hellinger distance matrix
        data_hel_rdf = rvegan.vegdist( rvegan.decostand( data_rdf, "hellinger"), "euclidean")
        # this should be reformatted to lower triangular matrix where x rows == y values
        print( data_hel_rdf, len( sample_df ) )
        data_hel_tri_rdf = _convert_to_dist_hel_matrix( data_hel_rdf, len( sample_rdf ) )

        # This is STEP 3
        # Setup formula for the adonis calculation.
        renv = rformula.environment
        renv['x'] = data_hel_tri_rdf # LHS
        renv['y'] = sample_rdf # RHS

        radonis = rvegan.adonis( rformula, permutations=999 )

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
        _pNTaxa = len(data_df.columns)
        _pFModelScale = -1.0 # this just a place holder

        premanova_df.loc[i] = [_pTest, _pOrder, _pAUC, _pSumOfSquares, _pMeanSquares, _pFModel, _pR2, _pPVal, _pNTaxa, _pFModelScale]

    # This is STEP 4.5
    # Convert F.model to a scaled version of F.model
    premanova_df["F.model.scale"] = (premanova_df["F.model"] - min(premanova_df["F.model"]))/(max(premanova_df["F.model"]) - min(premanova_df["F.model"]))
    premanova_df.index = range(1,len(premanova_df)+1)

    if( detailed ):
        premanova_df.to_csv( out_file )

    if( verbose and detailed ):
        dump.write( f"Output is Written in File: {str(out_file)} \n\n" )

    return premanova_df



# <><><> DEFINE EXECUTION FUNCTION <><><>
def main( auc_ordering, abundances, sample_df, centrality_type, name, detailed=False, verbose=False):
    """
    take in auc, abundance, and sample information and output values with permanova table
    :param auc_ordering: auc ordering generated
    :param abundances: abundances generated
    :param sample_df: samples with cooresponding types of invaded/natural
    :param centrality_type: used for generating detailed output ( betweenness, degree, eigenvector, closseness )
    :param name: name attached to all detailed output
    :param detailed: Output helper tables
    :param verbose: Output helper prints
    :return: table of all values generated during permanova calculation
    """

    out_dir = f"{os.path.dirname(os.path.realpath(__file__))}/output"
    # allows for cleaner execution and use of relative paths

    if( detailed ):
        out_file = f"{out_dir}/{name}_PERMANOVA_result.csv"
        # Create new files for output

        if( verbose ):

            dump = open(f"{out_dir}/step6_dump.txt", "w", encoding="utf-8")

            # Call PERMANOVA calculation
            permanova_df = perform_permanova( auc_ordering, abundances, sample_df, out_file, detailed, verbose, dump )
            dump.write( f"Plots generated to {out_dir}.\n" )
            dump.close()

        else:
            # Call PERMANOVA calculation
            permanova_df = perform_permanova( auc_ordering, abundances, sample_df, out_file, detailed )

        # Since this is detailed must generate plots
        _generate_figures( permanova_df, centrality_type, out_dir, name )

    elif( verbose ):

        dump = open(f"{out_dir}/step6_dump.txt", "w", encoding="utf-8")

        # Call PERMANOVA calculation
        permanova_df = perform_permanova( auc_ordering, abundances, sample_df, None, verbose=verbose, dump=dump )
        dump.close()

    else:
        # No excess files necessary just generate dataframe to pass on
        permanova_df = perform_permanova( auc_ordering, abundances, sample_df, None )

    permanova_df.reset_index( drop=True, inplace=True ) # reset indicis as a precautionary to make sure all df's start at index 0

    return permanova_df




