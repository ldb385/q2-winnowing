
import pandas as pd

import qiime2.plugin
import qiime2.plugin.model as model

from .plugin_setup import plugin

_EXPECTED_HEADERS_ = ['ab_comp', 'dataframe1', 'metric', 'centrality', 'total select', 'min count',
                    'smooth type', 'conditioning', 'keep threshold', 'correlation', 'weighted',
                    'correlation property', 'run time', 'kappa', 'agreement']
    # Placed here since is used in multiple places

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
        with self.open() as file:
            for line, idx in zip( file, range(2) ):
                delimited_line = list( line.rstrip("\n").split("\t") )
                if( idx == 0 ):
                    # first line, verify headers
                    for header in _EXPECTED_HEADERS_:
                        if( not delimited_line.__contains__( header ) ):
                            return False
                # make sure following rows aren't empty
                if( len( delimited_line ) <= 0 ): # must have atleast 2 data filled rows
                    return False
        return True # Passed test

    # <><> END OF CLASS <><>


# Define a directory format to be used
# this is necessary since multpile files will be stored in the directory
class WinnowedDirectoryFormat( model.DirectoryFormat ):
    # this is an example of a fixed layout since it will always include Feature ordering w/ Jaccard results
    # as well as complementary metadata files AUC values and Permanova files
    featureOrdering = model.File( r"feature_ordered.tsv", format=WinnowedFeatureOrderingFormat ) # Feature ordering w/ Jaccard
    auc = model.File( r"auc_ordered.tsv", format=model.TextFileFormat ) # AUC values with ordering
    permanova = model.File( r"permanova_ordered.tsv", format=model.TextFileFormat ) # PERMANOVA values with ordering

    # <><> END OF CLASS <><>


# Register the formats defined
plugin.register_formats( WinnowedFeatureOrderingFormat, WinnowedDirectoryFormat )

# Register directory format with semantic type
plugin.register_semantic_type_to_format(
    Winnowed,
    artifact_format=WinnowedDirectoryFormat
)


# <><><> Define transformers for reading in and writing out format information <><><>
# Define transformer to convert file format into dataframe
    # Standard to use a non-meaningful name for plugin transformer
@plugin.register_transformer
def _1( WinnowedFile: WinnowedFeatureOrderingFormat ) -> pd.DataFrame:
    # No Doc String since annotation describes functionality
    with WinnowedFile.open() as wf:
        expected_col = _EXPECTED_HEADERS_ + range(1, abs( len(_EXPECTED_HEADERS_) - max([ len(lwf) for lwf in wf ]) ))
        orderedFeatures_df = pd.DataFrame( columns=expected_col )
        for line in wf:
            row = line.rstrip("\n").split("\t")
            if( len(row) > len( expected_col ) ):
                raise ValueError( f"header does not have enough columns for row:\n\t{row}" )
            orderedFeatures_df.append( row, ignore_index=True )

        # Make sure dataframe has proper index values
        orderedFeatures_df.reset_index( drop=True, inplace=True )
        return orderedFeatures_df


# Define transformer to convert dataframe to file format
    # Standard to use a non-meaningful name for plugin transformer
@plugin.register_transformer
def _2( data: pd.DataFrame ) -> WinnowedFeatureOrderingFormat:
    # No Doc String since annotation describes functionality
    wf = WinnowedFeatureOrderingFormat()
    data.to_csv( wf, sep="\t", header=0 ) # pandas conveniantly converts to tsv

    return wf



















