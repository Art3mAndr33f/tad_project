#!/usr/bin/env python3
"""
update_rules.py
===============
Утилита для управления rules.md.

Команды:
  check                 Проверить структуру и версию rules.md
  bump [major|minor]    Обновить версию
  add-algorithm NAME    Сгенерировать чеклист + шаблон для нового алгоритма
  add-resolution RES    Сгенерировать чеклист для нового разрешения
  changelog             Добавить новую запись в Changelog (интерактивно)
  show-version          Показать текущую версию
  diff-check            Сравнить секции rules.md с реальным состоянием кода

Примеры:
  python scripts/update_rules.py check
  python scripts/update_rules.py bump minor
  python scripts/update_rules.py add-algorithm HiCseg
  python scripts/update_rules.py add-resolution 5000
  python scripts/update_rules.py changelog
"""

from __future__ import annotations

import argparse
import re
import sys
from datetime import date
from pathlib import Path

RULES_PATH = Path("rules.md")
ALGORITHMS_DIR = Path("src/algorithms")
CONFIG_PATH = Path("config/config.yaml")

# ─── Regex-паттерны ───────────────────────────────────────────────────────────

VERSION_PATTERN = re.compile(
    r"\*Версия rules\.md:\s*(\d+)\.(\d+)"
)
SECTION_PATTERN = re.compile(
    r"^## (\d+)\. (.+)$", re.MULTILINE
)
ALGO_REGISTRY_PATTERN = re.compile(
    r'(ALGORITHM_REGISTRY\s*=\s*\{[^}]+\})', re.DOTALL
)
CHANGELOG_ANCHOR = "## 19. Changelog"
VERSION_LINE_ANCHOR = "*Версия rules.md:"


# ─── Чтение / запись rules.md ────────────────────────────────────────────────

def read_rules() -> str:
    if not RULES_PATH.exists():
        print(f"[ERROR] {RULES_PATH} не найден. Запускай из корня проекта.")
        sys.exit(1)
    return RULES_PATH.read_text(encoding="utf-8")


def write_rules(content: str) -> None:
    RULES_PATH.write_text(content, encoding="utf-8")
    print(f"[OK] {RULES_PATH} обновлён.")


# ─── Версия ──────────────────────────────────────────────────────────────────

def get_version(content: str) -> tuple[int, int]:
    """Извлечь текущую версию из rules.md."""
    m = VERSION_PATTERN.search(content)
    if not m:
        print("[WARN] Строка версии не найдена, считаем 1.0")
        return 1, 0
    return int(m.group(1)), int(m.group(2))


def set_version(content: str, major: int, minor: int) -> str:
    """Обновить строку версии в rules.md."""
    new_version_str = f"*Версия rules.md: {major}.{minor}"
    result = VERSION_PATTERN.sub(new_version_str, content)
    if result == content:
        # Паттерн не сработал — заменяем вручную
        result = content.replace(
            VERSION_LINE_ANCHOR,
            f"{VERSION_LINE_ANCHOR} {major}.{minor}",
            1,
        )
    return result


def cmd_show_version(args: argparse.Namespace) -> None:
    content = read_rules()
    major, minor = get_version(content)
    print(f"Текущая версия rules.md: {major}.{minor}")


def cmd_bump(args: argparse.Namespace) -> None:
    kind = args.kind  # "major" или "minor"
    content = read_rules()
    major, minor = get_version(content)

    if kind == "major":
        major += 1
        minor = 0
    else:
        minor += 1

    content = set_version(content, major, minor)
    write_rules(content)
    print(f"Версия обновлена: {major}.{minor}")


# ─── Проверка структуры ──────────────────────────────────────────────────────

REQUIRED_SECTIONS = [
    "Что это за проект",
    "Структура проекта",
    "Данные",
    "Технический стек",
    "Алгоритмы",
    "Консенсус границ",
    "Метрики и статистика",
    "Конфигурация",
    "Соглашения по именованию",
    "Правила и запреты",
    "data_prep.py",
    "CLI пайплайна",
    "Тестирование",
    "Размеры хромосом",
    "Визуализация",
    "Типичные задачи",
    "Частые подводные камни",
    "Как обновлять этот файл",
    "Changelog",
]


def cmd_check(args: argparse.Namespace) -> None:
    """Проверить структуру rules.md."""
    content = read_rules()
    errors: list[str] = []
    warnings: list[str] = []

    # Проверить версию
    major, minor = get_version(content)
    print(f"[INFO] Версия: {major}.{minor}")

    # Проверить обязательные секции
    found_sections = {m.group(2) for m in SECTION_PATTERN.finditer(content)}
    for required in REQUIRED_SECTIONS:
        match_found = any(required.lower() in s.lower() for s in found_sections)
        if not match_found:
            errors.append(f"Отсутствует секция: «{required}»")

    # Проверить наличие CONSENSUS_COLORS
    if "#FFD700" not in content:
        warnings.append("Не найден цвет #FFD700 (CONSENSUS_COLORS для support=2)")
    if "#FF8C00" not in content:
        warnings.append("Не найден цвет #FF8C00 (CONSENSUS_COLORS для support=3)")
    if "#00C800" not in content:
        warnings.append("Не найден цвет #00C800 (CONSENSUS_COLORS для support=4)")

    # Проверить Changelog
    if CHANGELOG_ANCHOR not in content:
        errors.append("Отсутствует секция '## 19. Changelog'")

    # Проверить строку версии в конце
    if VERSION_LINE_ANCHOR not in content:
        errors.append("Отсутствует строка версии в конце файла")

    # Проверить соответствие алгоритмов в rules.md и коде
    _check_algorithms_sync(content, warnings)

    # Итог
    print()
    if warnings:
        print(f"⚠️  Предупреждения ({len(warnings)}):")
        for w in warnings:
            print(f"   - {w}")
    if errors:
        print(f"\n❌ Ошибки ({len(errors)}):")
        for e in errors:
            print(f"   - {e}")
        sys.exit(1)
    else:
        print("✅ rules.md корректен.")


def _check_algorithms_sync(content: str, warnings: list[str]) -> None:
    """Проверить что алгоритмы в rules.md соответствуют файлам в src/algorithms/."""
    if not ALGORITHMS_DIR.exists():
        warnings.append(f"Директория {ALGORITHMS_DIR} не найдена")
        return

    # Алгоритмы в коде
    code_algos = {
        f.stem.replace("run_", "")
        for f in ALGORITHMS_DIR.glob("run_*.py")
    }

    # Алгоритмы в rules.md (ищем шаблон "run_<algo>.py")
    rules_algos = set(re.findall(r'run_(\w+)\.py', content))
    rules_algos.discard("algorithm")   # шаблонное слово

    for algo in code_algos:
        if algo not in rules_algos:
            warnings.append(
                f"Алгоритм '{algo}' есть в коде (run_{algo}.py), "
                f"но не упомянут в rules.md"
            )

    for algo in rules_algos:
        if algo not in code_algos and algo not in ("algo", "new"):
            warnings.append(
                f"Алгоритм '{algo}' упомянут в rules.md, "
                f"но run_{algo}.py не найден в {ALGORITHMS_DIR}"
            )


# ─── Добавить алгоритм ───────────────────────────────────────────────────────

ALGORITHM_SECTION_TEMPLATE = """
### 5.{n} {AlgoName} (`src/algorithms/run_{algo}.py`)

- **Входные данные**: TODO — [dense numpy matrix / RAWobserved / другое]
- **Нормализация**: TODO — [RAW / KR / VC]
- **Параметры** (из конфига `algorithms.{algo}`):
  - `param1`: TODO — описание, диапазон
- **Выбор параметров**: TODO — описание стратегии подбора
- **Ограничения памяти**: TODO — [нет / указать хромосомы и разрешения]
- **Внешние зависимости**: TODO — [subprocess / pip / нет]
- **Fallback при ошибке**: возвращает `pd.DataFrame(columns=["chrom","start","end"])`
- **Особые требования**: TODO
"""

ALGORITHM_CONFIG_TEMPLATE = """  {algo}:
    # TODO: добавить параметры алгоритма {AlgoName}
    param1: default_value
    # Ограничения по хромосомам (если нужны):
    # chrom_limits:
    #   10000:  ["chr21", "chr22"]
    #   25000:  ["chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]
    #   50000:  []
    #   100000: []
"""

ALGORITHM_CHECKLIST_TEMPLATE = """
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  Чеклист: добавить алгоритм «{AlgoName}»
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

ШАГ 1 — Обновить rules.md:
  [ ] Секция 5.{n}: подраздел «{AlgoName}» добавлен ниже ↓
  [ ] Секция 5.6 (ALGORITHM_REGISTRY): добавить строку
        "{algo}": run_{algo},
  [ ] Секция 8 (config.yaml): добавить блок algorithms.{algo}
  [ ] Секция 15 (ALGO_COLORS): добавить цвет для «{algo}»
  [ ] Секция 19 (Changelog): добавить запись
  [ ] Обновить версию: python scripts/update_rules.py bump minor

ШАГ 2 — Создать код:
  [ ] src/algorithms/run_{algo}.py
        def run_{algo}(chrom, resolution, data_path, cfg, **kwargs) -> pd.DataFrame
  [ ] src/algorithms/__init__.py
        from .run_{algo} import run_{algo}
        ALGORITHM_REGISTRY["{algo}"] = run_{algo}
  [ ] config/config.yaml — блок algorithms.{algo}
  [ ] src/visualization.py — ALGO_COLORS["{algo}"] = "#XXXXXX"
  [ ] requirements.txt / environment.yml — новые зависимости (если есть)
  [ ] tests/test_{algo}.py — unit-тесты

ШАГ 3 — Проверить:
  [ ] python scripts/update_rules.py check
  [ ] pytest tests/ -v
  [ ] python pipeline/run_pipeline.py --resolution 100000 --chroms chr22 \\
        --algorithms {algo}

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Шаблон для секции 5.{n} (добавлен в rules.md):
{section}

Шаблон для config.yaml (добавь вручную):
{config}
"""


def _find_last_algorithm_section_number(content: str) -> int:
    """Найти номер последнего подраздела алгоритма в секции 5."""
    matches = re.findall(r"### 5\.(\d+)", content)
    return max((int(m) for m in matches), default=5)


def cmd_add_algorithm(args: argparse.Namespace) -> None:
    algo_name = args.name                        # например 'HiCseg'
    algo_key  = algo_name.lower().replace("-", "_")  # например 'hicseg'

    content = read_rules()
    n = _find_last_algorithm_section_number(content) + 1

    section = ALGORITHM_SECTION_TEMPLATE.format(
        n=n, AlgoName=algo_name, algo=algo_key
    )
    config_block = ALGORITHM_CONFIG_TEMPLATE.format(
        algo=algo_key, AlgoName=algo_name
    )

    # Вставить секцию перед 5.6 (ALGORITHM_REGISTRY) или в конец секции 5
    insert_anchor = "### 5.6 ALGORITHM_REGISTRY"
    if insert_anchor in content:
        content = content.replace(insert_anchor, section + "\n" + insert_anchor)
        write_rules(content)
    else:
        print(
            "[WARN] Якорь '### 5.6 ALGORITHM_REGISTRY' не найден. "
            "Добавь шаблон вручную."
        )

    print(
        ALGORITHM_CHECKLIST_TEMPLATE.format(
            n=n,
            AlgoName=algo_name,
            algo=algo_key,
            section=section.strip(),
            config=config_block.strip(),
        )
    )


# ─── Добавить разрешение ─────────────────────────────────────────────────────

RESOLUTION_CHECKLIST_TEMPLATE = """
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  Чеклист: добавить разрешение {res} bp ({res_kb} kb)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

ШАГ 1 — Обновить rules.md:
  [ ] Секция 3.3: добавить {res_kb} kb в таблицу «Разрешения»
  [ ] Секция 3.3: обновить таблицу scktld_limits для {res} bp
        — какие хромосомы допустимы? (учитывай RAM: ~n² * 4 bytes)
        — n = chrom_size / {res}
        — chr1 @ {res_kb}kb ≈ {chr1_bins} бинов → {chr1_mem_gb:.1f} GB
        — chr22 @ {res_kb}kb ≈ {chr22_bins} бинов → {chr22_mem_gb:.1f} GB
  [ ] Секция 8: добавить {res} в resolutions: [...]
  [ ] Секция 8: добавить {res}: [...] в scktld_limits
  [ ] Секция 19 (Changelog): добавить запись
  [ ] Обновить версию: python scripts/update_rules.py bump minor

ШАГ 2 — Обновить конфиг:
  [ ] config/config.yaml → resolutions: [..., {res}]
  [ ] config/config.yaml → chromosomes.scktld_limits.{res}: [...]

ШАГ 3 — Подготовить данные:
  [ ] bash download_data.sh  (если нужны новые RAWobserved-файлы)
  [ ] python pipeline/run_pipeline.py --prep-data --resolution {res}

ШАГ 4 — Проверить:
  [ ] python scripts/update_rules.py check
  [ ] python pipeline/run_pipeline.py --resolution {res} --chroms chr22 \\
        --algorithms armatus topdom

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Оценка памяти для scKTLD (dense float32, n×n матрица):
  chr1  @ {res_kb} kb: {chr1_bins:>6} бинов → {chr1_mem_gb:.2f} GB
  chr17 @ {res_kb} kb: {chr17_bins:>6} бинов → {chr17_mem_gb:.2f} GB
  chr21 @ {res_kb} kb: {chr21_bins:>6} бинов → {chr21_mem_gb:.2f} GB
  chr22 @ {res_kb} kb: {chr22_bins:>6} бинов → {chr22_mem_gb:.2f} GB

  Рекомендация: при >16 GB оставлять только chr21, chr22
                при 4–16 GB — chr17–chr22
                при <4 GB — все хромосомы
"""

CHROM_SIZES = {
    "chr1":  249_250_621,
    "chr17":  81_195_210,
    "chr21":  48_129_895,
    "chr22":  51_304_566,
}


def _mem_gb(n_bins: int) -> float:
    """Оценка памяти для dense float32 матрицы (GB)."""
    return (n_bins ** 2 * 4) / (1024 ** 3)


def cmd_add_resolution(args: argparse.Namespace) -> None:
    res = int(args.resolution)
    res_kb = res // 1000

    def bins(chrom: str) -> int:
        return int(CHROM_SIZES[chrom] / res) + 1

    print(
        RESOLUTION_CHECKLIST_TEMPLATE.format(
            res=res,
            res_kb=res_kb,
            chr1_bins=bins("chr1"),
            chr1_mem_gb=_mem_gb(bins("chr1")),
            chr17_bins=bins("chr17"),
            chr17_mem_gb=_mem_gb(bins("chr17")),
            chr21_bins=bins("chr21"),
            chr21_mem_gb=_mem_gb(bins("chr21")),
            chr22_bins=bins("chr22"),
            chr22_mem_gb=_mem_gb(bins("chr22")),
        )
    )


# ─── Добавить запись в Changelog ─────────────────────────────────────────────

CHANGELOG_ENTRY_TEMPLATE = """
### v{major}.{minor} — {today}

**Добавлено**
{added}

**Изменено**
{changed}

**Удалено**
{removed}

**Затронутые секции**: {sections}
**Затронутые файлы**: {files}

"""


def _prompt(label: str, default: str = "—") -> str:
    val = input(f"  {label} [{default}]: ").strip()
    return val if val else default


def _prompt_list(label: str) -> str:
    """Собрать многострочный список через ввод строк (пустая строка = конец)."""
    print(f"  {label} (пустая строка = завершить):")
    items = []
    while True:
        line = input("    - ").strip()
        if not line:
            break
        items.append(f"- {line}")
    return "\n".join(items) if items else "- —"


def cmd_changelog(args: argparse.Namespace) -> None:
    content = read_rules()
    major, minor = get_version(content)

    print(f"\nДобавление записи в Changelog (текущая версия: {major}.{minor})")
    print("─" * 50)

    bump_kind = input("Тип изменения [major/minor/none]: ").strip().lower()
    if bump_kind == "major":
        major += 1
        minor = 0
    elif bump_kind == "minor":
        minor += 1
    # else — версия не меняется

    print()
    added   = _prompt_list("Добавлено")
    changed = _prompt_list("Изменено")
    removed = _prompt_list("Удалено")
    sections = _prompt("Затронутые секции", "—")
    files    = _prompt("Затронутые файлы",  "—")

    entry = CHANGELOG_ENTRY_TEMPLATE.format(
        major=major,
        minor=minor,
        today=date.today().isoformat(),
        added=added,
        changed=changed,
        removed=removed,
        sections=sections,
        files=files,
    )

    # Вставить запись сразу после заголовка ## 19. Changelog
    if CHANGELOG_ANCHOR in content:
        content = content.replace(
            CHANGELOG_ANCHOR,
            CHANGELOG_ANCHOR + "\n" + entry,
            1,
        )
    else:
        print(f"[WARN] Якорь '{CHANGELOG_ANCHOR}' не найден, добавляю в конец.")
        content += "\n" + entry

    # Обновить версию
    content = set_version(content, major, minor)
    write_rules(content)
    print(f"\n✅ Запись добавлена. Новая версия: {major}.{minor}")


# ─── diff-check ──────────────────────────────────────────────────────────────

def cmd_diff_check(args: argparse.Namespace) -> None:
    """
    Проверить соответствие rules.md реальному состоянию проекта.
    Сравнивает алгоритмы, пути, разрешения.
    """
    import yaml  # type: ignore

    content = read_rules()
    issues: list[str] = []

    # 1. Алгоритмы в __init__.py vs rules.md
    init_path = ALGORITHMS_DIR / "__init__.py"
    if init_path.exists():
        init_text = init_path.read_text()
        registry_algos = set(re.findall(r'"(\w+)":\s*run_\w+', init_text))
        rules_algos = set(re.findall(r'"(\w+)":\s*run_\w+', content))
        for a in registry_algos - rules_algos:
            issues.append(
                f"Алгоритм '{a}' в ALGORITHM_REGISTRY кода, но не в rules.md"
            )
        for a in rules_algos - registry_algos:
            issues.append(
                f"Алгоритм '{a}' в rules.md, но не в ALGORITHM_REGISTRY кода"
            )

    # 2. Разрешения в config.yaml vs rules.md
    if CONFIG_PATH.exists():
        with open(CONFIG_PATH) as f:
            cfg = yaml.safe_load(f)

        config_res = set(cfg.get("resolutions", []))
        rules_res  = set(
            int(m) for m in re.findall(r'\b(5000|10000|25000|50000|100000)\b', content)
        )
        for r in config_res - rules_res:
            issues.append(
                f"Разрешение {r} bp есть в config.yaml, но не упомянуто в rules.md"
            )

    # 3. Пути из config.yaml vs rules.md
    if CONFIG_PATH.exists():
        paths_section = cfg.get("paths", {})
        for key, val in paths_section.items():
            if isinstance(val, str) and val not in content:
                issues.append(
                    f"Путь '{val}' (config.paths.{key}) не найден в rules.md"
                )

    # Итог
    print()
    if issues:
        print(f"⚠️  Расхождения rules.md с кодом/конфигом ({len(issues)}):")
        for issue in issues:
            print(f"   - {issue}")
        print(
            "\nОбнови rules.md вручную или через:"
            "\n  python scripts/update_rules.py changelog"
        )
    else:
        print("✅ rules.md синхронизирован с кодом и конфигом.")


# ─── CLI ─────────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Утилита управления rules.md",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("check",        help="Проверить структуру rules.md")
    sub.add_parser("show-version", help="Показать текущую версию")
    sub.add_parser("diff-check",   help="Сравнить rules.md с кодом/конфигом")
    sub.add_parser("changelog",    help="Добавить запись в Changelog (интерактивно)")

    p_bump = sub.add_parser("bump", help="Обновить версию")
    p_bump.add_argument(
        "kind",
        choices=["major", "minor"],
        help="Тип изменения версии",
    )

    p_algo = sub.add_parser("add-algorithm", help="Чеклист + шаблон для нового алгоритма")
    p_algo.add_argument("name", help="Название алгоритма (например: HiCseg)")

    p_res = sub.add_parser("add-resolution", help="Чеклист для нового разрешения")
    p_res.add_argument("resolution", type=int, help="Разрешение в bp (например: 5000)")

    return parser


COMMANDS = {
    "check":         cmd_check,
    "show-version":  cmd_show_version,
    "bump":          cmd_bump,
    "add-algorithm": cmd_add_algorithm,
    "add-resolution": cmd_add_resolution,
    "changelog":     cmd_changelog,
    "diff-check":    cmd_diff_check,
}


def main() -> None:
    parser = build_parser()
    args   = parser.parse_args()
    COMMANDS[args.command](args)


if __name__ == "__main__":
    main()