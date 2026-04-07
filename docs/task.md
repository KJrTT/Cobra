## Структура проекта
```
cobra-demo/
├── metabolic_analysis.py    # Загрузка модели, основные расчёты, демонстрационные сценарии
└── README.md                # Документация проекта
```

## Описание модулей
### `metabolic_analysis.py` — основной модуль системы
Содержит демонстрацию работы с метаболической моделью *E. coli* core metabolism с использованием библиотеки **cobra** (COBRApy).  
Включает загрузку модели, запуск **Flux Balance Analysis (FBA)**, изменение внешней среды, симуляцию **knockout** генов и реакций, анализ изменчивости потоков (**FVA**) и интерпретацию результатов.

**Основные объекты библиотеки, которые демонстрируются:**
- **`Model`** — метаболическая модель в целом
- **`Reaction`** — биохимическая реакция
- **`Metabolite`** — метаболит
- **`Gene`** — ген
- **`Solution`** — результат оптимизации (FBA)

**Демонстрационные сценарии:**
- Загрузка и инспекция модели
- Запуск базового FBA и просмотр summary
- Изменение состава среды (аэробные / анаэробные условия)
- Симуляция knockout реакций и генов (с использованием контекстного менеджера)
- Сравнение роста до и после knockout
- Запуск Flux Variability Analysis (FVA)
- Анализ входящих/исходящих потоков

## Как использовать
### Установка и запуск
1. Убедитесь, что у вас установлен Python 3.8 или выше.

2. Установите библиотеку COBRApy:
   ```bash
   pip install cobra
   ```

3. Скачайте файл `metabolic_analysis.py`.

4. Запустите демонстрацию:
   ```bash
   python metabolic_analysis.py
   ```

### Работа с приложением
При запуске скрипт последовательно выполняет несколько демонстрационных блоков:
1. Загрузка модели и базовая информация
2. Запуск FBA в стандартных условиях + summary
3. Изменение среды (аэробный и анаэробный рост)
4. Симуляция knockout реакции и гена (с безопасным контекстом)
5. Сравнение скорости роста до и после удаления
6. Запуск Flux Variability Analysis (FVA)
7. Подробный вывод потоков веществ

Каждый блок сопровождается подробным выводом в консоль, показывающим скорость роста (objective value), входящие/исходящие потоки и изменения после модификаций модели.

## Основные объекты COBRApy
### 1. Model — метаболическая модель
Содержит все реакции, метаболиты и гены.  
**Ключевые атрибуты:**
- `model.reactions` — список реакций
- `model.metabolites` — список метаболитов
- `model.genes` — список генов
- `model.objective` — целевая функция (обычно биомасса)

**Пример:**
```python
from cobra.io import load_model

model = load_model("textbook")   # E. coli core model
print(model)
print(f"Реакций: {len(model.reactions)}")
print(f"Метаболитов: {len(model.metabolites)}")
print(f"Генов: {len(model.genes)}")
```

### 2. Reaction — реакция
Представляет биохимическую реакцию с нижней и верхней границами потока (`lower_bound`, `upper_bound`).

**Пример доступа:**
```python
reaction = model.reactions.get_by_id("PFK")  # Фосфофруктокиназа
print(reaction)
print(reaction.bounds)
```

### 3. Gene — ген
Связан с реакциями через Gene-Protein-Reaction (GPR) правила.

**Пример:**
```python
gene = model.genes.get_by_id("b0720")
```

## Демонстрационные сценарии
### 1. Базовый FBA (Flux Balance Analysis)
Демонстрирует запуск оптимизации для максимального роста клетки.
```python
solution = model.optimize()
print("Скорость роста:", solution.objective_value)
print(model.summary())
```

### 2. Изменение внешней среды (medium)
Показывает, как ограничивать поступление веществ (глюкоза, кислород и др.).
```python
model.medium = {"EX_glc__D_e": 10, "EX_o2_e": 1000}  # аэробные условия
# или анаэробные: EX_o2_e = 0
solution = model.optimize()
```

### 3. Симуляция knockout
Демонстрирует удаление реакции или гена и влияние на рост.
```python
# Knockout реакции
with model as mutant:
    mutant.reactions.PFK.knock_out()
    solution = mutant.optimize()
    print("Рост после knockout PFK:", solution.objective_value)

# Knockout гена
with model as mutant:
    mutant.genes.b0720.knock_out()
    solution = mutant.optimize()
```

### 4. Flux Variability Analysis (FVA)
Анализ диапазона возможных потоков для выбранных реакций при оптимальном росте.
```python
from cobra.flux_analysis import flux_variability_analysis

fva_result = flux_variability_analysis(model, model.reactions[:10])
print(fva_result)
```

### 5. Сравнение сценариев
Скрипт сравнивает:
- Нормальный рост (аэробный)
- Анаэробный рост
- Рост после удаления ключевой реакции/гена

## Ключевые концепции COBRApy
### 1. Загрузка моделей
- `load_model("textbook")` — встроенная тестовая модель
- `read_sbml_model("model.xml")` — загрузка из SBML-файла (BiGG, и др.)

### 2. Оптимизация и решение
- `model.optimize()` → возвращает объект `Solution`
- `solution.objective_value` — значение целевой функции (рост)
- `model.summary()` — удобный отчёт по in/out fluxes

### 3. Модификация модели
- Изменение границ: `reaction.lower_bound = -10`
- Knockout: `reaction.knock_out()` или `gene.knock_out()`
- Контекстный менеджер `with model as mutant:` — безопасные временные изменения

### 4. Анализ изменчивости (FVA)
```python
fva = flux_variability_analysis(model, fraction_of_optimum=0.95)
```
Показывает минимальный и максимальный возможный поток для каждой реакции.

### 5. Объекты и доступ
- Поиск по id: `model.reactions.get_by_id("ID")`
- Итерация: `for reaction in model.reactions:`
- Связь ген-реакция: GPR правила автоматически учитываются при knockout генов

### 6. Вывод результатов
```python
print(model.summary())
print(solution.fluxes)          # все потоки
print(solution.status)          # optimal / infeasible и т.д.
```

## Примеры вывода
При запуске скрипта вы увидите:
- Информацию о модели (количество реакций, метаболитов, генов)
- Скорость роста в разных условиях
- Таблицы in/out fluxes из `summary()`
- Изменение роста после knockout
- Диапазоны потоков из FVA

Проект наглядно показывает, как с помощью COBRApy моделировать метаболизм клетки, предсказывать влияние генетических модификаций и условий среды на рост микроорганизма.

Это отличный учебный стенд для знакомства с **Constraint-Based Metabolic Modeling** и методом **Flux Balance Analysis**.

Готов к использованию в образовательных целях!
