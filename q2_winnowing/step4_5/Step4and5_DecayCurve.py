import pandas as pd
import numpy as np
import itertools
from rpy2.robjects.packages import importr
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
from rpy2.robjects import r
from scipy.spatial.distance import pdist, squareform
from rpy2.robjects import pandas2ri
from rpy2.robjects import IntVector, Formula
import rpy2.robjects as robjects
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

# import biom

pandas2ri.activate()

utils = rpackages.importr('utils')
utils.chooseCRANmirror(ind=1)

packnames = ('PKNCA', 'colorRamps')
names_to_install = []
# names_to_install = [x for packnames if not rpackages.isinstalled(x)]

for x in packnames:
    if (rpackages.isinstalled(x) == False):
        names_to_install.append(x)

if len(names_to_install) > 0:
    utils.install_packages(StrVector(names_to_install))

rpkcna = importr('PKNCA')
raucx = r['pk.calc.auc']
rc = r['c']


def calc_auc_percentile(input_file, output_file, parameter_file, verbose=False, dump=None ):
    brome_dg = pd.read_csv(input_file)
    brome_dg = brome_dg.sort_values("metric", axis=0, ascending=False)
    brome_dg.index = range(1, len(brome_dg) + 1)
    brome_auc = raucx(brome_dg["metric"], brome_dg.index, interval=rc(1, len(brome_dg["metric"])))
    result_df = pd.DataFrame(columns=['auc', 'otu.num'])
    parameter_df = pd.DataFrame(columns=['x', 'y'])

    for factor in np.arange(0.01, 1.00, 0.01):
        area = 0.0
        end_range = 2
        while (area <= round(factor, 2) * brome_auc):
            area = raucx(brome_dg["metric"], brome_dg.index, interval=rc(1, end_range))
            end_range += 1

        if( verbose ):
            dump.write( f"The point at which we reach {str(round(factor * 100, 2))}% of the AUC is = {str(end_range)}\n" )

        result_df.loc[int(round(factor * 100, 2))] = ["auc" + str(int(round(factor * 100, 2)))] + [end_range]

    result_df.loc[100] = ["auc100"] + [len(brome_dg["metric"])]
    parameter_df['x'] = brome_dg.index - 1
    parameter_df['y'] = brome_dg["metric"]
    parameter_df.loc[len(brome_dg)] = [len(brome_dg)] + [parameter_df.iloc[len(brome_dg) - 1, 1]]
    result_df.to_csv(output_file)
    parameter_df.iloc[1:, :].to_csv(parameter_file)


# def main_silent( inFiles: [biom.Table], parameterFiles: [biom.Table] ):
#
#     if (len(inFiles) != len(parameterFiles)):
#         raise Exception("Each input file must have a corresponding parameter file.")
#
#     outFiles = [f"{x.split('.')[:-2].split('/')[-1]}_auc_result.csv" for x in inFiles]
#
#     # Create new files for output
#     for x in outFiles:
#         x = open(x, "w+", encoding="utf-8")
#         x.close()
#
#
#     for inp, out, param in zip(inFiles, outFiles, parameterFiles):
#         calc_auc_percentile(inp, out, param)
#
#     return


def main_verbose(inFiles, parameterFiles):

    if (len(inFiles) != len(parameterFiles)):
        raise Exception(f"Each input file must have a corresponding parameter file.\n"
                        f"you have provided {len(inFiles)} input and {len(parameterFiles)} parameters")

    outFiles = []

    # Extract names
    for x in parameterFiles:
        path = "".join( x.split('.')[:-1] )
        name = {path.split('/')[-1]}
        outFiles.append( f"output/{list(name)[0]}_auc_result.csv")

    # Create new files for output
    for x in outFiles:
        x = open(x, "w+", encoding="utf-8")
        x.close()

    # Since this is verbose we must also write to a dump
    dump = open("output/step4_5_dump.txt", "w", encoding="utf-8")

    for inp, out, param in zip(inFiles, outFiles, parameterFiles):
        dump.write(f"\n\nProcessing Input file: {inp}\n")
        calc_auc_percentile(inp, out, param, True, dump)
        dump.write(f"Output is written in file: {out}\n")
        dump.write(f"Parameters are written in file: {param}\n")

    dump.close()

    return


# <><><> TEST <><><>
main_verbose(["test_data/ADD1_AUC100_MIC0.2_Brome_bacfunarc_dw_otu_table-graph_centrality-degree-selectallbyall.csv"],
             ["test_data/brome.add1.degree.fig2-new.csv"])

# main_verbose(["test_data/ADD1_AUC100_MIC0.2_Brome_bacfunarc_dw_otu_table-graph_centrality-degree-selectallbyall.csv",
#               "test_data/ADD1_AUC100_MIN3_Brome_bacfunarc_dw_otu_table-graph_centrality-closeness-selectallbyall.csv",
#               "test_data/ADD1_AUC100_MIN3Brome_bacfunarc_dw_otu_table-graph_centrality-eigenvector-selectallbyall.csv",
#               "test_data/ADD1_AUC100_MIN3_MIC0.2__Brome_bacfunarc_dw_otu_table-graph_centrality-betweenness-selectallbyall.csv"],
#              ["test_data/brome.add1.degree.fig2-new.csv", "test_data/brome.add1.closeness-new.csv",
#               "test_data/brome.add1.eigenvector-new.csv", "test_data/brome.add1.betweenness-new.csv"])
