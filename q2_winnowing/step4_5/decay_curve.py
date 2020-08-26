

# <><><> SETUP IMPORTS <><><>

import os
import pandas as pd
import numpy as np
from rpy2.robjects.packages import importr
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
from rpy2.robjects import r
from rpy2.robjects import pandas2ri

pandas2ri.activate()

utils = rpackages.importr('utils')
utils.chooseCRANmirror(ind=1)

packnames = ['PKNCA', 'colorRamps']
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


# <><><> DEFINE FUNCTIONS <><><>


def calc_auc_percentile( input_df ):
    """
    calculate each AUC percentile of decay curve
    :param input_df: Important features and OTUs that have been passed on to be ordered
    :return: ordered AUC percentiles of which OTU is most influential for percentile
    """

    input_df = input_df.sort_values("metric", axis=0, ascending=False)
    input_df.index = range(1, len(input_df) + 1)
    input_auc = raucx(input_df["metric"], input_df.index, interval=rc(1, len(input_df["metric"])))
    result_df = pd.DataFrame(columns=['auc', 'otu.num'])
    parameter_df = pd.DataFrame(columns=['x', 'y'])

    for factor in np.arange(0.01, 1.00, 0.01):
        area = 0.0
        end_range = 2
        # 1. calculate the area of each trapezoid
        while (area <= round(factor, 2) * input_auc):
            area = raucx(input_df["metric"], input_df.index, interval=rc(1, end_range))
            end_range += 1

        print( f"The point at which we reach {str(round(factor * 100, 2))}% of the AUC is = {str(end_range)}\n" )

        #2. sum trapezoid areas to get AUC
        result_df.loc[int(round(factor * 100, 2))] = ["auc" + str(int(round(factor * 100, 2)))] + [end_range]

    result_df.loc[100] = ["auc100"] + [len(input_df["metric"])]
    parameter_df['x'] = input_df.index - 1
    parameter_df['y'] = input_df["metric"]
    parameter_df.loc[len(input_df)] = [len(input_df)] + [parameter_df.iloc[len(input_df) - 1, 1]]

    return result_df, parameter_df.iloc[1:, :]



# <><><> DEFINE EXECUTION FUNCTION <><><>


def main( input_df, name, detailed=False ):
    """
    Each OTU is now ordered by centrality and the AUC of each is calculated.
    :param input_df: Important features and OTUs that have been passed on to be ordered
    :param name: name attached to all detailed output
    :param detailed: Output helper tables
    :return: ordered AUC percentiles of which OTU is most influential for percentile
    """

    out_dir = f"{os.path.dirname(os.path.realpath(__file__))}/output"
    # allows for cleaner execution and use of relative paths


    if( detailed ):
        out_file = f"{out_dir}/{name}_auc_result.csv"
        parameter_file = f"{out_dir}/{name}_auc_parameter.csv"

        # Create new files for output
        out_file = open( out_file, "w+", encoding="utf-8")
        parameter_file = open( parameter_file, "w+", encoding="utf-8" )

        print(f"Processing Input dataFrame: {name}")
        result, param = calc_auc_percentile(input_df )
        print(f"Output is written in file: {out_file}")
        print(f"Parameters are written in file: {parameter_file}")

        # Write to CSV since this is detailed
        result.to_csv(out_file)
        param.to_csv(parameter_file)

        out_file.close()
        parameter_file.close()


    else:
        print(f"Processing Input dataFrame: {name}")
        result, param = calc_auc_percentile( input_df )

    # Return results dataframe along with the parameters dataframe
    return result, param



