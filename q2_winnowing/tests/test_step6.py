
from unittest import TestCase, main as unittest_main
import pandas as pd

from q2_winnowing.step6.Step6_Permanova import main as step6_main

class Step6Tests( TestCase ):

    # <><> read input values for testing <><>
    auc_df_in = pd.read_csv("./sample_data/step6/test_in_AUCs.csv")
    abundances_df_in = pd.read_csv("./sample_data/step6/test_in_abundances.csv")
    samples_df_in = pd.read_csv("./sample_data/step6/test_in_samples.csv")
    # <><> read output values for comparison <><>
    permanova_df_out = pd.read_csv("./sample_data/step6/test_out_PERMANOVA.csv", index_col=0 )

    def test_step6_main(self):

        permanova_output = step6_main( self.auc_df_in, self.abundances_df_in, self.samples_df_in, "", False, False )

        for row in range(0, len(permanova_output) ):
            self.assertEqual(
                permanova_output.loc[ permanova_output.index[[row]], "test" ].values.item(),
                self.permanova_df_out.loc[ self.permanova_df_out.index[[row]], "test" ].values.item()
            )
            self.assertEqual(
                int( permanova_output.loc[ permanova_output.index[[row]], "order" ].values.item() ),
                int( self.permanova_df_out.loc[ self.permanova_df_out.index[[row]], "order" ].values.item() )
            )
            self.assertEqual(
                int( permanova_output.loc[ permanova_output.index[[row]], "auc" ].values.item() ),
                int( self.permanova_df_out.loc[ self.permanova_df_out.index[[row]], "auc" ].values.item() )
            )
            self.assertEqual(
                str( permanova_output.loc[ permanova_output.index[[row]], "SumsOfSqs" ].values.item() )[:15],
                str( self.permanova_df_out.loc[ self.permanova_df_out.index[[row]], "SumsOfSqs" ].values.item() )[:15]
            )
            self.assertEqual(
                str( permanova_output.loc[ permanova_output.index[[row]], "MeanSqs" ].values.item() )[:15],
                str( self.permanova_df_out.loc[ self.permanova_df_out.index[[row]], "MeanSqs" ].values.item() )[:15]
            )
            self.assertEqual(
                str( permanova_output.loc[ permanova_output.index[[row]], "F.model" ].values.item() )[:15],
                str( self.permanova_df_out.loc[ self.permanova_df_out.index[[row]], "F.model" ].values.item() )[:15]
            )
            self.assertEqual(
                str( permanova_output.loc[ permanova_output.index[[row]], "R2" ].values.item() )[:15],
                str( self.permanova_df_out.loc[ self.permanova_df_out.index[[row]], "R2" ].values.item() )[:15]
            )
            self.assertEqual(
                int( permanova_output.loc[ permanova_output.index[[row]], "N.taxa" ].values.item() ),
                int( self.permanova_df_out.loc[ self.permanova_df_out.index[[row]], "N.taxa" ].values.item() )
            )
            self.assertEqual(
                str( permanova_output.loc[ permanova_output.index[[row]], "F.model.scale" ].values.item() )[:15],
                str( self.permanova_df_out.loc[ self.permanova_df_out.index[[row]], "F.model.scale" ].values.item() )[:15]
            )


if __name__ == '__main__':
    unittest_main()