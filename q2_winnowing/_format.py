
import qiime2.plugin.model as model

# <><><> DEFINE STATIC HEADERS <><><>
_EXPECTED_FEATURE_HEADERS_ = ['ab_comp', 'dataframe1', 'metric', 'centrality', 'iteration select', 'total select', 'min count',
                    'smooth type', 'conditioning', 'keep threshold', 'correlation', 'weighted',
                    'correlation property', 'run time', 'kappa', 'agreement']
_EXPECTED_AUC_HEADERS_ = ["auc", "otu.num"]
_EXPECTED_PERMANOVA_HEADERS_ = ["test", "order", "auc", "SumsOfSqs", "MeanSqs", "F.model", "R2", "Pval",
                                "N.taxa", "F.model.scale"]
    # Placed here since is used in multiple places

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










