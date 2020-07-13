
import pandas as pd
import shutil

import qiime2.plugin
import qiime2.plugin.model as model
from qiime2.sdk import Artifact

from q2_winnowing.plugin_setup import plugin

# <><><> DEFINE STATIC HEADERS FOR SNIFFING <><><>
_EXPECTED_FEATURE_HEADERS_ = ['ab_comp', 'dataframe1', 'metric', 'centrality', 'total select', 'min count',
                    'smooth type', 'conditioning', 'keep threshold', 'correlation', 'weighted',
                    'correlation property', 'run time', 'kappa', 'agreement']
_EXPECTED_AUC_HEADERS_ = ["auc", "otu.num"]
_EXPECTED_PERMANOVA_HEADERS_ = ["test", "order", "auc", "SumsOfSqs", "MeanSqs", "F.model", "R2", "Pval",
                                "N.taxa", "F.model.scale"]
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
                    for header in _EXPECTED_FEATURE_HEADERS_:
                        if( not delimited_line.__contains__( header ) ):
                            return False
                # make sure following rows aren't empty
                if( len( delimited_line ) <= 0 ): # must have atleast 2 data filled rows
                    return False
        return True # Passed test

    # <><> END OF CLASS <><>


# Define a file format representing the AUC ordered
# most accurate ordering returned
class WinnowedAucOrderingFormat( model.TextFileFormat ):

    # this will just make sure that the file is in a expected format
    def sniff(self):
        with self.open() as file:
            if( len(list(file)) != 101 ):
                return False # 1% of AUC to 100%
            file.readline() # ignore header
            for line, _ in zip( file, range(5) ):
                delimited_line = line.rstrip('\n').split('\t')
                if len(delimited_line) != 3:
                    return False # file must be 3 columns

            return True # Test passed

    # <><> END OF CLASS <><>


# Define a file format representing the PERMANOVA ordered
# most accurate ordering returned
class WinnowedPermanovaOrderingFormat(model.TextFileFormat):

    # this will just make sure that the file is in a expected format
    def sniff(self):
        with self.open() as file:
            if( len(list(file)) != 101 ):
                return False  # AUC output was used {1% of AUC to 100%}
            for line, idx in zip(file, range(3)):
                delimited_line = list(line.rstrip("\n").split("\t"))
                if (idx == 0):
                    # first line, verify headers
                    for header in _EXPECTED_PERMANOVA_HEADERS_:
                        if (not delimited_line.__contains__(header)):
                            return False
                # make sure following rows aren't empty
                if (len(delimited_line) <= 0):  # row must contain data
                    return False

            return True  # Test passed

    # <><> END OF CLASS <><>


# Define a directory format to be used
# this is necessary since multpile files will be stored in the directory
class WinnowedDirectoryFormat( model.DirectoryFormat ):
    # this is an example of a fixed layout since it will always include Feature ordering w/ Jaccard results
    # as well as complementary metadata files AUC values and Permanova files
    featureOrdering = model.File( r"feature_ordered.tsv", format=WinnowedFeatureOrderingFormat ) # Feature ordering w/ Jaccard
    auc = model.File( r"auc_ordered.tsv", format=WinnowedAucOrderingFormat ) # AUC values with ordering
    permanova = model.File( r"permanova_ordered.tsv", format=WinnowedPermanovaOrderingFormat ) # PERMANOVA values with ordering

    # <><> END OF CLASS <><>


# Register the formats defined
plugin.register_formats( WinnowedFeatureOrderingFormat, WinnowedAucOrderingFormat,
                         WinnowedPermanovaOrderingFormat, WinnowedDirectoryFormat )

# Register directory format with semantic type
plugin.register_semantic_type_to_format(
    Winnowed,
    artifact_format=WinnowedDirectoryFormat
)


# <><><> Define transformers for reading in and writing out format information <><><>
# Define transformer to convert featureOrdering file format into dataframe
    # Standard to use a non-meaningful name for plugin transformer
@plugin.register_transformer
def _1( WinnowedFile: WinnowedFeatureOrderingFormat ) -> pd.DataFrame:
    # No Doc String since annotation describes functionality
    with WinnowedFile.open() as wf:
        expected_col = _EXPECTED_FEATURE_HEADERS_ + \
                       list( range(1, abs( len(_EXPECTED_FEATURE_HEADERS_) - max([ len(lwf.rstrip("\n").split("\t")) for lwf in wf ]) )) )
        wf.seek( 0 ) # Reset read cursor to start of data
        wf.readline() # skip headers
        orderedFeatures_df = pd.DataFrame()
        for line in wf:
            row = line.rstrip("\n").split("\t")
            if( str(row[0]).isdigit() ):
                row = row[1:] # index column can be removed
            if( len(row) > len( expected_col ) ):
                raise ValueError( f"FeatureOrdering: header does not have enough columns for row:\n\t{row}" )
            orderedFeatures_df = orderedFeatures_df.append( [row], ignore_index=True )

        # Make sure dataframe has proper index and column values
        orderedFeatures_df.columns = expected_col # Append changes columns so got to reset
        orderedFeatures_df.reset_index( drop=True, inplace=True )
        return orderedFeatures_df


# Define transformer to convert featureOrdering dataframe to file format
    # Standard to use a non-meaningful name for plugin transformer
@plugin.register_transformer
def _2( data: pd.DataFrame ) -> WinnowedFeatureOrderingFormat:
    # No Doc String since annotation describes functionality
    wf = WinnowedFeatureOrderingFormat()
    path = str( wf.path )
    data.to_csv( path, sep="\t" ) # pandas conveniantly converts to tsv

    return wf


# Define transformer to convert auc file format into dataframe
    # Standard to use a non-meaningful name for plugin transformer
@plugin.register_transformer
def _3( AucFile: WinnowedAucOrderingFormat ) -> pd.DataFrame:
    # No Doc String since annotation describes functionality
    with AucFile.open() as wf:
        expected_col = _EXPECTED_AUC_HEADERS_
        auc_df = pd.DataFrame()
        wf.readline() # skip headers
        for line in wf:
            row = line.rstrip("\n").split("\t")
            if( str(row[0]).isdigit() ):
                row = row[1:] # index column can be removed
            if( len(row) > len( expected_col ) ):
                raise ValueError( f"AUC: header does not have enough columns for row:\n\t{row}" )
            auc_df = auc_df.append( [row], ignore_index=True )

        # Make sure dataframe has proper index values
        auc_df.columns = expected_col # Append changes columns so got to reset
        auc_df.reset_index( drop=True, inplace=True )
        return auc_df


# Define transformer to convert auc dataframe to file format
    # Standard to use a non-meaningful name for plugin transformer
@plugin.register_transformer
def _4( data: pd.DataFrame ) -> WinnowedAucOrderingFormat:
    # No Doc String since annotation describes functionality
    wf = WinnowedAucOrderingFormat()
    path = str( wf.path )
    data.to_csv( path, sep="\t" ) # pandas conveniantly converts to tsv

    return wf


# Define transformer to convert permanova file format into dataframe
    # Standard to use a non-meaningful name for plugin transformer
@plugin.register_transformer
def _5( PermanovaFile: WinnowedPermanovaOrderingFormat ) -> pd.DataFrame:
    # No Doc String since annotation describes functionality
    with PermanovaFile.open() as wf:
        expected_col = _EXPECTED_PERMANOVA_HEADERS_
        permanova_df = pd.DataFrame()
        wf.readline() # skip headers
        for line in wf:
            row = line.rstrip("\n").split("\t")
            if( str(row[0]).isdigit() ):
                row = row[1:] # index column can be removed
            if( len(row) > len( expected_col ) ):
                raise ValueError( f"PERMANOVA: header does not have enough columns for row:\n\t{row}" )
            permanova_df = permanova_df.append( [row], ignore_index=True )

        # Make sure dataframe has proper index values
        permanova_df.columns = expected_col # Append changes columns so got to reset
        permanova_df.reset_index( drop=True, inplace=True )
        return permanova_df


# Define transformer to convert permanova dataframe to file format
    # Standard to use a non-meaningful name for plugin transformer
@plugin.register_transformer
def _6( data: pd.DataFrame ) -> WinnowedPermanovaOrderingFormat:
    # No Doc String since annotation describes functionality
    wf = WinnowedPermanovaOrderingFormat()
    path = str( wf.path )
    data.to_csv( path, sep="\t" ) # pandas conveniantly converts to tsv

    return wf


# <><><> Define directory transformers for reading in and writing out format information <><><>
# Define transformer to convert permanova file format into dataframe
    # Standard to use a non-meaningful name for plugin transformer
@plugin.register_transformer
def _7( dirfmt: WinnowedDirectoryFormat ) -> tuple:
    # No Doc String since annotation describes functionality

        # Ordered feature read of Dir
    artifact_ordered_features = dirfmt.featureOrdering.view( WinnowedFeatureOrderingFormat )
    # artifact_ordered_features = Artifact.import_data(Winnowed, f_orderedFeature )
        # PERMANOVA read of Dir
    artifact_permanova = dirfmt.permanova.view( WinnowedPermanovaOrderingFormat )
    # artifact_permanova = Artifact.import_data(Winnowed, f_permanova )
        # AUC feature read of Dir
    artifact_auc = dirfmt.auc.view( WinnowedAucOrderingFormat )
    # artifact_auc = Artifact.import_data(Winnowed, f_auc )

    dirTuple = (
        artifact_ordered_features.view(pd.DataFrame),
        artifact_auc.view(pd.DataFrame),
        artifact_permanova.view(pd.DataFrame)
    )

    return dirTuple


# # Define transformer to convert permanova dataframe to file format
#     # Standard to use a non-meaningful name for plugin transformer
@plugin.register_transformer
def _8( data: tuple ) -> WinnowedDirectoryFormat:
    # No Doc String since annotation describes functionality
    result = WinnowedDirectoryFormat()
    path = result.path
    of, oa, op = data

    features_fp = str(path / "feature_ordered.tsv")
    permanova_fp = str(path / "permanova_ordered.tsv")
    auc_fp = str(path / "auc_ordered.tsv")

    shutil.copyfile( _2( of ).path, features_fp ) # feature ordering
    shutil.copyfile( _4( oa ).path, auc_fp ) # AUC
    shutil.copyfile( _6( op ).path, permanova_fp ) # PERMANOVA

    return result











