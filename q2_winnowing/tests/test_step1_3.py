
from unittest import TestCase, main as unittest_main
import os
import pandas as pd
import numpy as np

from q2_winnowing.step1_3.pipeline import main as step1_3_main

class Step1_3Tests( TestCase ):
    # <><><> Testing class for Step 1 to 3 <><><>

    package = 'q2_winnowing.tests'

    # <><> read input values and expected output values for testing <><>
    testing_data = []
    testing_data.append((
        f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step1_3/test_in_data.csv",
        f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step1_3/test_out_features.csv",
        f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step1_3/test_out_important.csv",
        f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step1_3/test_out_abundances.csv",
    ))

    def test_step1_3_main_graph_bw(self):

        for data_path, features_path, important_path, adundances_path in self.testing_data:
            data_in = pd.read_csv( data_path )
            data_in.name = "TESTING_FRAME"
            features_out = pd.read_csv( features_path )
            important_out = pd.read_csv( important_path, index_col=0 )
            abundances_out = pd.read_csv( adundances_path )

            features_result, important_result, abundances_result = \
                step1_3_main( False, data_in, None, metric_name="graph_centrality", c_type="add_one", min_count=3,
                              total_select=25, iteration_select=128, pca_components=4, smooth_type="sliding_window",
                              window_size=3, centrality_type="betweenness", keep_threshold=0.5, correlation="spearman",
                              weighted=True, corr_prop="both", evaluation_type="kl_divergence", min_connected=0,
                              detailed=False, verbose_p=False )

            run_out = features_out.pop("run time")
            run_result = features_result.pop("run time")
            self.assertLess( abs( run_result.tail(1).iloc[0] - run_out.tail(1).iloc[0] ), 100 )
                # test that time for step is reasonably close to expected ( within range of 100 sec )

            # Remove names column
            features_out.pop("dataframe1")
            features_result.pop("dataframe1")

            np.testing.assert_array_equal(
                features_result.to_numpy().astype(object),
                features_out.to_numpy().astype(object),
            )
            np.testing.assert_array_almost_equal(
                important_result.to_numpy().astype(object),
                important_out.to_numpy().astype(object),
                15
            )
            np.testing.assert_array_equal(
                abundances_result.to_numpy().astype(object),
                abundances_out.to_numpy().astype(object),
            )






if __name__ == '__main__':
    unittest_main()
