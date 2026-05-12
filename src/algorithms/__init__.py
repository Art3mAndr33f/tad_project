"""
algorithms package
==================
Все алгоритмы возвращают pd.DataFrame с колонками: chrom, start, end
(start и end в базовых парах, координатная система hg19).
"""

from .run_armatus        import run_armatus
from .run_topdom         import run_topdom
from .run_scktld         import run_scktld
from .run_coitad         import run_coitad
from .run_dihmm          import run_dihmm
from .run_ontad          import run_ontad
from .run_modularity_tad import run_modularity_tad

ALGORITHM_REGISTRY = {
    "armatus":        run_armatus,
    "topdom":         run_topdom,
    "scktld":         run_scktld,
    "coitad":         run_coitad,
    "dihmm":          run_dihmm,
    "ontad":          run_ontad,
    "modularity_tad": run_modularity_tad,
}

__all__ = ["run_armatus", "run_topdom", "run_scktld", "run_coitad",
           "run_dihm", "run_ontad", "run_modularity_tad", 
           "ALGORITHM_REGISTRY"]