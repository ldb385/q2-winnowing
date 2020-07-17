
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
        pd.read_csv(f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step4_5/test_in_data.csv"),
        pd.read_csv(f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step4_5/test_out_features.csv",
                    index_col=0),
        pd.read_csv(f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step4_5/test_out_important.csv",
                    index_col=0),
        pd.read_csv(f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/step4_5/test_out_abundances.csv",
                    index_col=0),
    ))

