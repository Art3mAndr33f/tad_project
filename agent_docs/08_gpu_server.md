# GPU-сервер brain-lab — правила и настройка

---

## Текущий статус GPU (актуально ~2026-05-14)

```
GPU 0: 73989/81154 MiB занято — ❌ НЕ использовать
GPU 1: 55443/81154 MiB занято — ❌ НЕ использовать
GPU 2: 76549/81154 MiB занято — ❌ НЕ использовать
GPU 3: 41114/81154 MiB занято — ✅ ~40 GB свободно
```

**Проверить актуальный статус:**
```bash
nvidia-smi --query-gpu=index,memory.used,memory.free --format=csv
```

---

## Обязательные переменные окружения

```bash
export CUDA_VISIBLE_DEVICES=3   # GPU 3
export OMP_NUM_THREADS=2        # ограничение CPU threads
export MKL_NUM_THREADS=2        # ограничение MKL threads
```

⚠️ **Не нагружать CPU** — сервер общий. Тяжёлые вычисления только на GPU.

---

## GPU по компонентам

| Компонент | CPU | GPU | Ускорение |
|-----------|-----|-----|-----------|
| scKTLD (n<5000) | scipy.linalg.eigh | torch.linalg.eigh (CUDA) | ~10–20× |
| TopDom, Armatus | numpy | — | — |
| DI+HMM | hmmlearn | — | — |
| OnTAD | subprocess C++ | — | — |
| ModularityTAD | numpy | — | — |

---

## Установка torch (CUDA 12.1)

```bash
pip install torch --index-url https://download.pytorch.org/whl/cu121
python -c "import torch; print(torch.cuda.is_available())"  # → True
```

---

## Типичные команды запуска

```bash
# Активация окружения
conda activate tad_pipeline

# Полный запуск
python pipeline/run_pipeline.py \
    --resolution 25000 \
    --chroms chr17 chr18 chr19 chr20 chr21 chr22 \
    --algorithms armatus topdom scktld coitad dihmm ontad modularity_tad \
    --force

# Запуск в фоне с логом
nohup python pipeline/run_pipeline.py \
    --resolution 25000 \
    --chroms chr17 chr18 chr19 chr20 chr21 chr22 \
    --algorithms armatus topdom scktld coitad dihmm ontad modularity_tad \
    > logs/pipeline_25kb.log 2>&1 &

echo "PID: $!"
```
