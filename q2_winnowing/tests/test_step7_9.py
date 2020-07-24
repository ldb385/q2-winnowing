
from unittest import TestCase, main as unittest_main
import pandas as pd
import numpy as np
import os

from q2_winnowing.step7_9.jaccard import main as step7_9_main
from q2_winnowing.step7_9.jaccard import jaccard_coefficient

class Step7_9Tests( TestCase ):
    # <><><> Testing class for Step 7 to 9 <><><>

    package = 'q2_winnowing.tests'

    # <><> read input values and expected output values for testing <><>
    testing_data = []
    testing_data.append((
        f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step7_9/test_in_metric_results_graph_spearman_0.5_bw.csv",
        f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step7_9/test_out_jaccard_graph_spearman_0.5_bw.csv"
    ))
    testing_data.append((
        f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step7_9/test_in_metric_results_graph_spearman_0.5_ei.csv",
        f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step7_9/test_out_jaccard_graph_spearman_0.5_ei.csv"
    ))
    testing_data.append((
        f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step7_9/test_in_metric_results_graph_MIC_0.5_cl.csv",
        f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step7_9/test_out_jaccard_graph_MIC_0.5_cl.csv",
    ))
    testing_data.append((
        f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step7_9/test_in_metric_results_graph_MIC_0.5_dg.csv",
        f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step7_9/test_out_jaccard_graph_MIC_0.5_dg.csv"
    ))


    def test_jaccard_coefficient(self):
        a = (0,1,2,5,6)
        b = (0,2,3,4,5,7,9)
        out = 0.3333333 # this is repeating
        jaccard_coefficient_result = jaccard_coefficient(a,b)
        np.testing.assert_almost_equal(
            jaccard_coefficient_result,
            out,
            decimal=7
        )


    def test_step7_9_main(self):

        for metric_results_path, jaccard_path in self.testing_data:

            metric_results_in = pd.read_csv( metric_results_path )
            jaccard_out = pd.read_csv( jaccard_path, index_col=0 )

            jaccard_result = step7_9_main( metric_results_in, "", False, False )

            for row in range(0, len(jaccard_result)):
                # <><> verify tables are in proper format <><>
                self.assertEqual(
                    jaccard_result.loc[ jaccard_result.index[[row]], "conditioning" ].values.item(),
                    jaccard_out.loc[ jaccard_out.index[[row]], "conditioning" ].values.item()
                )
                self.assertEqual(
                    jaccard_result.loc[jaccard_result.index[[row]], "centrality"].values.item(),
                    jaccard_out.loc[jaccard_out.index[[row]], "centrality"].values.item()
                )
                self.assertEqual(
                    jaccard_result.loc[jaccard_result.index[[row]], "correl"].values.item(),
                    jaccard_out.loc[jaccard_out.index[[row]], "correl"].values.item()
                )
                self.assertEqual(
                    jaccard_result.loc[jaccard_result.index[[row]], "threshold"].values.item(),
                    jaccard_out.loc[jaccard_out.index[[row]], "threshold"].values.item()
                )
                self.assertEqual(
                    jaccard_result.loc[jaccard_result.index[[row]], "select_iter"].values.item(),
                    jaccard_out.loc[jaccard_out.index[[row]], "select_iter"].values.item()
                )
                # <><> verify output is similar <><>
                np.testing.assert_almost_equal(
                    jaccard_result.loc[jaccard_result.index[[row]], "kappa"].values.item(),
                    jaccard_out.loc[jaccard_out.index[[row]], "kappa"].values.item(),
                    decimal=15
                )
                np.testing.assert_almost_equal(
                    jaccard_result.loc[jaccard_result.index[[row]], "agreement"].values.item(),
                    jaccard_out.loc[jaccard_out.index[[row]], "agreement"].values.item(),
                    decimal=15
                )


if __name__ == '__main__':
    unittest_main()