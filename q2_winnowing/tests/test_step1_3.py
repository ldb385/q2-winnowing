
from unittest import TestCase, main as unittest_main
import os
import pandas as pd
import numpy as np

from q2_winnowing.step1_3.Step1_3_Pipeline import main as step1_3_main

class Step1_3Tests( TestCase ):

    package = 'q2_winnowing.tests'

    # <><> read input values and expected output values for testing <><>
    testing_data = []
    testing_data.append((
        f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step1_3/test_in_data.csv",
        f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step1_3/test_out_features.csv",
        f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step1_3/test_out_important.csv",
        f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step1_3/test_out_abundances.csv",
    ))

    def test_step1_3_main(self):

        for data_path, features_path, important_path, adundances_path in self.testing_data:
            data_in = pd.read_csv( data_path )
            features_out = pd.read_csv( features_path, index_col=0 )
            important_out = pd.read_csv( important_path, index_col=0 )
            adundances_out = pd.read_csv( adundances_path, index_col=0 )






if __name__ == '__main__':
    unittest_main()

