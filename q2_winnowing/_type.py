

from qiime2.plugin import SemanticType

# Since type of output produced by q2_winnowing is only supported as metadata a new
# semantic type will be defined in order to account for data passed by q2_winnowing plugin

# Define semantic type
Winnowed = SemanticType("Winnowed")

