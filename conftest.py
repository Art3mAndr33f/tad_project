from __future__ import annotations

"""Root conftest.py — добавляет корень проекта в sys.path.

Нужно чтобы pytest мог находить пакет src/ без установки через pip install -e .
Не трогать: этот файл нужен всем воркстримам.
"""

import sys
from pathlib import Path

# Добавить корень проекта (где лежит src/) в начало sys.path
_PROJECT_ROOT = str(Path(__file__).resolve().parent)
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)
