__version__ = "2020.0.0"

# Make the types defined in this plugin importable from the top-level package
# so they can be easily imported by other plugins relying on these types.
from ._format import Winnowed

__all__ = ["Winnowed"]