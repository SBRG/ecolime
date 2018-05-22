from distutils.version import StrictVersion

from ecolime import (chaperones, dna_replication, generics, flat_files,
                     modifications, ribosome, transcription, translation,
                     translocation, trna_charging, corrections)
import cobrame

__version__ = "0.0.9"


if StrictVersion(__version__) != StrictVersion(cobrame.__version__):
    raise UserWarning('COBRAme version (%s) and ECOLIme version (%s) are not'
                      ' equivalent' % (cobrame.__version__, __version__))
