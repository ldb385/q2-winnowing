
import pkg_resources
import unittest
import uuid

import pandas as pd
import qiime2.core.archive as archive
from qiime2.sdk import Artifact

from q2_winnowing import Winnowed

input_to_dataframe = pd.DataFrame(
    data=[[0,False,"-name-_1_1_","graph_centrality","betweenness",25,3,"sliding_window","add_one",0.5,"spearman",True,"both",2760.5702665330027,0.09688013136288998,0.47058823529411764,"Otu5330","Otu4626","Otu6188","Otu6434","Otu4238","Otu5737","Otu3245","Otu5484","Otu6390",'Otu6747',"Otu6264","Otu6602","Otu6306",'Otu6284','Otu5981',"Otu6671","Otu4727","Otu6326","Otu6199","Otu4164","Otu5267","Otu6344","Otu5313","Otu0728","Otu5584"],
          [1,False,"-name-_1_4_","graph_centrality","betweenness",25,3,"sliding_window","add_one",0.5,"spearman",True,"both",421.0375310439995,0.13651315789473684,0.5151515151515151,"Otu5330","Otu4626","Otu6434","Otu6188","Otu4238","Otu5737","Otu3245","Otu5484","Otu6390",'Otu6326',"Otu6747","Otu6602","Otu6264",'Otu6284','Otu4727',"Otu5981","Otu6671","Otu6306","Otu4164","Otu4706","Otu5270","Otu6199","Otu6344","Otu5267","Otu0728"],
          [2,False,"-name-_1_16_","graph_centrality","betweenness",25,3,"sliding_window","add_one",0.5,"spearman",True,"both",128.871447255995,0.66996699669967,0.6129032258064516,"Otu5330","Otu4626","Otu6434","Otu6188","Otu5737","Otu3245","Otu6326","Otu4238","Otu5484","Otu6747","Otu5538","Otu6264","Otu4706","Otu4164","Otu6284","Otu6210","Otu6390","Otu6602","Otu5981","Otu6306","Otu6199","Otu6671","Otu5270","Otu6448","Otu6544"],
          [3,False,"-name-_1_64_","graph_centrality","betweenness",25,3,"sliding_window","add_one",0.5,"spearman",True,"both",68.40203193300113,1.0,1.0,"Otu5330","Otu4626","Otu6434","Otu6188","Otu5737","Otu3245","Otu6326","Otu4238","Otu5484","Otu6747","Otu5538","Otu6264","Otu4706","Otu4164","Otu6284","Otu6210","Otu3798","Otu6602","Otu5520","Otu6697","Otu6736","Otu6390","Otu1133","Otu6306","Otu6279"],
          [4,False,"-name-_1_128_","graph_centrality","betweenness",25,3,"sliding_window","add_one",0.5,"spearman",True,"both",68.50248471400118,1.0,1.0,"Otu5330","Otu4626","Otu6434",'Otu6188',"Otu5737","Otu3245","Otu6326","Otu4238","Otu5484","Otu6747","Otu5538","Otu6264","Otu4706","Otu4164","Otu6284","Otu6210","Otu3798","Otu6602","Otu5520","Otu6697","Otu6736","Otu6390","Otu1133","Otu6306","Otu6279"]],
    columns=["","ab_comp","dataframe1","metric","centrality","total select","min count","smooth type","conditioning","keep threshold","correlation","weighted","correlation property","run time","kappa","agreement",1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]
)


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


    def test_reader_transformer(self):
        fp = pkg_resources.resource_filename(
            'q2_winnowing.tests', 'sample_data/test_in_dir')

        artifact = Artifact.import_data(Winnowed, fp)
        # `Artifact.view` invokes the transformer that handles the
        # `WinnowedFormat` -> `dataframe` transformation.
        pd.testing.assert_frame_equal( artifact.view(pd.DataFrame), input_to_dataframe )


    def test_writer_transformer(self):
        # `Artifact._from_view` invokes transformer that handles `dataframe` ->
        # `WinnowedFormat`, because the `WinnowedDirectoryFormat` has
        # been registered as the directory format for the semantic type.
        # Directory formatting isn't created since single file formatting is invoked
        artifact = Artifact._from_view( Winnowed, input_to_dataframe,
                                       pd.DataFrame, archive.ImportProvenanceCapture())
        # Test that the directory and file format can be read again.
        self.assertEqual(artifact.view(pd.DataFrame), input_to_dataframe)



if __name__ == "__main__":
    unittest.main()