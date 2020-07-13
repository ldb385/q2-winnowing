
import pkg_resources
import unittest
import uuid

import pandas as pd
import qiime2.core.archive as archive
from qiime2.sdk import Artifact

from q2_winnowing import Winnowed

exp_featureOrdering = pd.DataFrame(
    data=[[False,"-name-_1_1_","graph_centrality","betweenness",25,3,"sliding_window","add_one",0.5,"spearman",True,"both",2760.5702665330027,0.09688013136288998,0.47058823529411764,"Otu5330","Otu4626","Otu6188","Otu6434","Otu4238","Otu5737","Otu3245","Otu5484","Otu6390",'Otu6747',"Otu6264","Otu6602","Otu6306",'Otu6284','Otu5981',"Otu6671","Otu4727","Otu6326","Otu6199","Otu4164","Otu5267","Otu6344","Otu5313","Otu0728","Otu5584"],
          [False,"-name-_1_4_","graph_centrality","betweenness",25,3,"sliding_window","add_one",0.5,"spearman",True,"both",421.0375310439995,0.13651315789473684,0.5151515151515151,"Otu5330","Otu4626","Otu6434","Otu6188","Otu4238","Otu5737","Otu3245","Otu5484","Otu6390",'Otu6326',"Otu6747","Otu6602","Otu6264",'Otu6284','Otu4727',"Otu5981","Otu6671","Otu6306","Otu4164","Otu4706","Otu5270","Otu6199","Otu6344","Otu5267","Otu0728"],
          [False,"-name-_1_16_","graph_centrality","betweenness",25,3,"sliding_window","add_one",0.5,"spearman",True,"both",128.871447255995,0.66996699669967,0.6129032258064516,"Otu5330","Otu4626","Otu6434","Otu6188","Otu5737","Otu3245","Otu6326","Otu4238","Otu5484","Otu6747","Otu5538","Otu6264","Otu4706","Otu4164","Otu6284","Otu6210","Otu6390","Otu6602","Otu5981","Otu6306","Otu6199","Otu6671","Otu5270","Otu6448","Otu6544"],
          [False,"-name-_1_64_","graph_centrality","betweenness",25,3,"sliding_window","add_one",0.5,"spearman",True,"both",68.40203193300113,1.0,1.0,"Otu5330","Otu4626","Otu6434","Otu6188","Otu5737","Otu3245","Otu6326","Otu4238","Otu5484","Otu6747","Otu5538","Otu6264","Otu4706","Otu4164","Otu6284","Otu6210","Otu3798","Otu6602","Otu5520","Otu6697","Otu6736","Otu6390","Otu1133","Otu6306","Otu6279"],
          [False,"-name-_1_128_","graph_centrality","betweenness",25,3,"sliding_window","add_one",0.5,"spearman",True,"both",68.50248471400118,1.0,1.0,"Otu5330","Otu4626","Otu6434",'Otu6188',"Otu5737","Otu3245","Otu6326","Otu4238","Otu5484","Otu6747","Otu5538","Otu6264","Otu4706","Otu4164","Otu6284","Otu6210","Otu3798","Otu6602","Otu5520","Otu6697","Otu6736","Otu6390","Otu1133","Otu6306","Otu6279"]],
    columns=["ab_comp","dataframe1","metric","centrality","total select","min count","smooth type","conditioning","keep threshold","correlation","weighted","correlation property","run time","kappa","agreement",1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]
)
exp_permanova = pd.DataFrame(
    data=[["auc1",1,1,1.517518014,1.517518014,5.186579863,0.084766625,0.004,3,0.572329632],
          ["auc2",2,2,1.517518014,1.517518014,5.186579863,0.084766625,0.011,3,0.572329632],
          ["auc3",3,3,1.517518014,1.517518014,5.186579863,0.084766625,0.004,3,0.572329632],
          ["auc4",4,4,1.517518014,1.517518014,5.186579863,0.084766625,0.007,3,0.572329632],
          ["auc5",5,5,1.517518014,1.517518014,5.186579863,0.084766625,0.009,3,0.572329632]],
    columns=["test","order",'auc',"SumsOfSqs",'MeanSqs',"F.model","R2","Pval","N.taxa","F.model.scale"]
) # Suffices to test with first five lines of file
exp_auc = pd.DataFrame(
    data=[["auc1",3],
          ["auc2",3],
          ["auc3",3],
          ["auc4",3],
          ["auc5",3]],
    columns=["auc","otu.num"]
) # Suffices to test with first five lines of file


class WinnowedFormatTests(unittest.TestCase):
    # test class layout referenced from:
    #   https://github.com/qiime2-graveyard/q2-dummy-types/blob/master/q2_dummy_types/tests/test_mapping.py

    def test_data_import(self):
        fp = pkg_resources.resource_filename(
            'q2_winnowing.tests', 'sample_data/test_in_dir')

        # `Artifact.import_data` copies `test_in_featureData.tsv` into the artifact after
        # performing validation on the file.
        artifact = Artifact.import_data(Winnowed, fp)

        self.assertEqual(artifact.type, Winnowed)
        self.assertIsInstance(artifact.uuid, uuid.UUID)

        # <><> END OF FUNCTION <><>


    def test_reader_transformer(self):
        fp = pkg_resources.resource_filename(
            'q2_winnowing.tests', 'sample_data/test_in_dir')

        artifact = Artifact.import_data(Winnowed, fp)
        featureOrdering_df, permanova_df, auc_df = artifact.view( type((pd.DataFrame(),pd.DataFrame(),pd.DataFrame())) )
        # `Artifact.view` invokes the transformer that handles the
        # `WinnowedFormat` -> `dataframe` transformation.
        # print( featureOrdering_df, exp_featureOrdering )
        pd.testing.assert_frame_equal(
            featureOrdering_df.astype(type("")),
            exp_featureOrdering.astype( type("")),
            check_dtype=False
        ) # Avoid checking values since reading df stores as objects while, hard coding in does not
        # ex) bool(False) == Object(False) in pandas is False although the values function the same.
        pd.testing.assert_frame_equal(
            permanova_df.iloc[:5].astype(type("")),
            exp_permanova.astype(type("")),
            check_dtype=False,
            check_exact=False
        )
        pd.testing.assert_frame_equal(
            auc_df.iloc[:5].astype(type("")),
            exp_auc.astype(type("")),
            check_dtype=False
        )

        # <><> END OF FUNCTION <><>



    def test_writer_transformer(self):
        # `Artifact._from_view` invokes transformer that handles `dataframe` ->
        # `WinnowedFormat`, because the `WinnowedDirectoryFormat` has
        # been registered as the directory format for the semantic type.
        # Directory formatting isn't created since single file formatting is invoked
        artifact = Artifact._from_view(Winnowed, exp_featureOrdering,
                                       pd.DataFrame, archive.ImportProvenanceCapture())
        # Test that the directory and file format can be read again.
        self.assertEqual(artifact.view(pd.DataFrame), exp_featureOrdering)

        # <><> END OF FUNCTION <><>



if __name__ == "__main__":
    unittest.main()