# Сборка OnTAD из исходников (C++)

> Нужно только при первоначальной настройке сервера.
> После сборки бинарник лежит в tools/ontad/OnTAD.

---

## Требования

- conda env `tad_pipeline` активирован
- Нет sudo (сервер brain-lab)

---

## Шаг 1: Зависимости через conda

```bash
conda install -c conda-forge libcurl gcc gxx make -y
```

## Шаг 2: Клонировать исходники

```bash
cd tools/
git clone https://github.com/anlin00007/OnTAD.git ontad_src
cd ontad_src/src
```

## Шаг 3: Патч straw.cpp (GCC 15 строже)

```bash
python3 -c "
text = open('straw.cpp').read()
old = '#include "straw.h"'
new = '#include "straw.h"\n#include <cstdint>'
if '<cstdint>' not in text:
    text = text.replace(old, new, 1)
    open('straw.cpp', 'w').write(text)
    print('Patch applied')
else:
    print('Already patched — skip')
"
```

## Шаг 4: Компиляция

```bash
g++ -std=c++11 \
    -I${CONDA_PREFIX}/include \
    -L${CONDA_PREFIX}/lib \
    -Wl,-rpath,${CONDA_PREFIX}/lib \
    main.cpp step1.cpp step2.cpp step3.cpp step4.cpp common.cpp straw.cpp \
    -lm -lcurl -lz -o OnTAD
```

## Шаг 5: Установка

```bash
mkdir -p ~/tad_project/tools/ontad
cp OnTAD ~/tad_project/tools/ontad/OnTAD
chmod +x ~/tad_project/tools/ontad/OnTAD
```

## Проверка

```bash
~/tad_project/tools/ontad/OnTAD --help
# Должен вывести usage без ошибок
```

---

## ❌ Частые ошибки

| Ошибка | Причина | Решение |
|--------|---------|---------|
| `curl/curl.h: No such file` | libcurl не в PATH GCC | conda install libcurl + явный -I${CONDA_PREFIX}/include |
| `uint64_t is not defined` | GCC 15 строже | Применить патч straw.cpp (#include <cstdint>) |
| `make` не находит libcurl | make не передаёт conda пути | Использовать явный g++ (не make) |
| Двойной `#include <cstdint>` | Патч применён дважды | `grep -c cstdint straw.cpp` перед патчем |
| 0 TADs в subprocess | Неправильный парсинг depth | depth = parts[2], не parts[3] |
