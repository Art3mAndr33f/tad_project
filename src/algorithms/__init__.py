"""
algorithms package
==================
Все алгоритмы возвращают pd.DataFrame с колонками: chrom, start, end
(start и end в базовых парах, координатная система hg19).
"""

from .run_armatus import run_armatus
from .run_topdom  import run_topdom
from .run_scktld  import run_scktld
from .run_coitad  import run_coitad

ALGORITHM_REGISTRY = {
    "armatus": run_armatus,
    "topdom":  run_topdom,
    "scktld":  run_scktld,
    "coitad":  run_coitad,
}

__all__ = ["run_armatus", "run_topdom", "run_scktld", "run_coitad",
           "ALGORITHM_REGISTRY"]