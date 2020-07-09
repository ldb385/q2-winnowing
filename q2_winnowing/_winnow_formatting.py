

import qiime2.plugin
import qiime2.plugin.model as model

from .plugin_setup import plugin

# Since type of output produced by q2_winnowing is only supported as metadata a new
# semantic type will be defined in order to account for data passed by q2_winnowing plugin

# Define semantic type
Winnowed = qiime2.plugin.SemanticType("Winnowed")
# Register type
plugin.register_semantic_types(Winnowed)

# Multiple file formats will need to be used to account for different output
# Final output will be stored in the form of a directory format

# Define a file format representing the bulk of the output
# This will be the feature ordering
class WinnowedFeatureOrderingFormat( model.TextFileFormat ):

    # this will just make sure that the file is in a expected format
    def sniff(self):
        expected_headers = ['ab_comp','dataframe1','metric','centrality','total select','min count',
                            'smooth type','conditioning','keep threshold','correlation','weighted',
                            'correlation property','run time','kappa','agreement']
        with self.open() as file:
            for line, idx in zip( file, range(2) ):
                delimited_line = list( line.rstrip("\n").split("\t") )
                if( idx == 0 ):
                    # first line, verify headers
                    for header in expected_headers:
                        if( not delimited_line.__contains__( header ) ):
                            return False
                # make sure following rows aren't empty
                if( len( delimited_line ) <= 0 ): # must have atleast 2 data filled rows
                    return False
        return True # Passed test







# Define transformers for reading in and writing out format information