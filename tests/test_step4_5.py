
from unittest import TestCase, main as unittest_main
import os
import pandas as pd
import numpy as np

from q2_winnowing.step4_5.decay_curve import main as step4_5_main

class Step4_5Tests( TestCase ):
    # <><><> Testing class for Step 4 to 5 <><><>

    package = 'q2_winnowing.tests'

    # <><> read input values and expected output values for testing <><>
    testing_data = {
        "degree": (
        f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step4_5/test_in_MIC_0.2_dw_graph_centrality_dg.csv",
        f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step4_5/test_out_dg_param.csv",
        f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step4_5/test_out_dg_auc.csv"
        ),
        "closeness": (
        f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step4_5/test_in_MIN_3_dw_graph_centrality_cl.csv",
        f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step4_5/test_out_cl_param.csv",
        f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step4_5/test_out_cl_auc.csv"
        ),
        "eigenvector": (
        f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step4_5/test_in_MIN_3_dw_graph_centrality_ei.csv",
        f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step4_5/test_out_ei_param.csv",
        f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step4_5/test_out_ei_auc.csv"
        ),
        "betweeness": (
        f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step4_5/test_in_MIN_3_MIC_0.2_dw_graph_centrality_bw.csv",
        f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step4_5/test_out_bw_param.csv",
        f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step4_5/test_out_bw_auc.csv"
        )
    }


    def test_step4_5_main_degree(self):

        important_features_path, parameter_path, auc_path = self.testing_data["degree"]
        important_features_in = pd.read_csv( important_features_path )
        auc_result, auc_param = step4_5_main( important_features_in, "", False )
        parameter_out = pd.read_csv( parameter_path, index_col=0 )
        auc_out = pd.read_csv( auc_path, index_col=0 )
        # Read in data as needed to avoid unnecessary memory consumption

        np.testing.assert_array_equal(
            auc_result.to_numpy().astype(object),
            auc_out.to_numpy().astype(object),
        )
        np.testing.assert_array_almost_equal(
            auc_param.to_numpy().astype(type(0.0)),
            parameter_out.to_numpy().astype(type(0.0)),
            decimal=15
        )


    def test_step4_5_main_closeness(self):

        important_features_path, parameter_path, auc_path = self.testing_data["closeness"]
        important_features_in = pd.read_csv( important_features_path )
        auc_result, auc_param = step4_5_main( important_features_in, "", False )
        parameter_out = pd.read_csv( parameter_path, index_col=0 )
        auc_out = pd.read_csv( auc_path, index_col=0 )
        # Read in data as needed to avoid unnecessary memory consumption

        np.testing.assert_array_equal(
            auc_result.to_numpy().astype(object),
            auc_out.to_numpy().astype(object),
        )
        np.testing.assert_array_almost_equal(
            auc_param.to_numpy().astype(type(0.0)),
            parameter_out.to_numpy().astype(type(0.0)),
            decimal=15
        )


    def test_step4_5_main_eigenvector(self):

        important_features_path, parameter_path, auc_path = self.testing_data["eigenvector"]
        important_features_in = pd.read_csv( important_features_path )
        auc_result, auc_param = step4_5_main( important_features_in, "", False )
        parameter_out = pd.read_csv( parameter_path, index_col=0 )
        auc_out = pd.read_csv( auc_path, index_col=0 )
        # Read in data as needed to avoid unnecessary memory consumption

        np.testing.assert_array_equal(
            auc_result.to_numpy().astype(object),
            auc_out.to_numpy().astype(object),
        )
        np.testing.assert_array_almost_equal(
            auc_param.to_numpy().astype(type(0.0)),
            parameter_out.to_numpy().astype(type(0.0)),
            decimal=15
        )


    def test_step4_5_main_betweeness(self):

        important_features_path, parameter_path, auc_path = self.testing_data["betweeness"]
        important_features_in = pd.read_csv( important_features_path )
        auc_result, auc_param = step4_5_main( important_features_in, "", False )
        parameter_out = pd.read_csv( parameter_path, index_col=0 )
        auc_out = pd.read_csv( auc_path, index_col=0 )
        # Read in data as needed to avoid unnecessary memory consumption

        np.testing.assert_array_equal(
            auc_result.to_numpy().astype(object),
            auc_out.to_numpy().astype(object),
        )
        np.testing.assert_array_almost_equal(
            auc_param.to_numpy().astype(type(0.0)),
            parameter_out.to_numpy().astype(type(0.0)),
            decimal=15
        )



if __name__ == '__main__':
    unittest_main()