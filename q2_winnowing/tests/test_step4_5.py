
from unittest import TestCase, main as unittest_main
import os
import pandas as pd

from q2_winnowing.step4_5.Step4and5_DecayCurve import main as step4_5_main

class Step4_5Tests( TestCase ):

    # <><> read input values and expected output values for testing <><>
    testing_data = []
    testing_data.append((
        pd.read_csv(f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step4_5/test_in_MIC_0.2_dw_graph_centrality_dg.csv"),
        pd.read_csv(f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step4_5/test_out_dg_param.csv"),
        pd.read_csv(f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step4_5/test_out_dg_auc.csv")
    ))
    testing_data.append((
        pd.read_csv(f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step4_5/test_in_MIN_3_dw_graph_centrality_cl.csv"),
        pd.read_csv(f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step4_5/test_out_cl_param.csv"),
        pd.read_csv(f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step4_5/test_out_cl_auc.csv")
    ))
    testing_data.append((
        pd.read_csv(f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step4_5/test_in_MIN_3_dw_graph_centrality_ei.csv"),
        pd.read_csv(f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step4_5/test_out_ei_param.csv"),
        pd.read_csv(f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step4_5/test_out_ei_auc.csv")
    ))
    testing_data.append((
        pd.read_csv(f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step4_5/test_in_MIN_3_MIC_0.2_dw_graph_centrality_bw.csv"),
        pd.read_csv(f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step4_5/test_out_bw_param.csv"),
        pd.read_csv(f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step4_5/test_out_bw_auc.csv")
    ))


    def test_step4_5_main(self):

        for important_features_in, parameter_out, auc_out in self.testing_data:
            auc_result, auc_param = step4_5_main( important_features_in, "", False, False )

            self.assertEqual(
                auc_result.to_numpy(),
                auc_out.values.to_numpy()
            )
            self.assertEqual(
                auc_param.values.to_numpy(),
                parameter_out.values.to_numpy()
            )



if __name__ == '__main__':
    unittest_main()