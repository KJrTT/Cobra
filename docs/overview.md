# COBRApy

COBRApy (Constraint-Based Reconstruction and Analysis in Python) — это библиотека Python для работы с моделями метаболических сетей на основе ограничений (Constraint-Based Reconstruction and Analysis). Она предоставляет инструменты для загрузки, анализа и симуляции генетически-масштабных метаболических моделей (Genome-Scale Metabolic Models, GSMM) в форматах SBML, JSON и MAT. COBRApy позволяет исследователям и биоинженерам проводить flux balance analysis (FBA), анализ knockout-мутаций, sampling решений и другие вычисления без необходимости использования MATLAB или других проприетарных сред .

---
## Модуль 1. Контент и предпосылки

### Ограничения стандартных средств + Описание проблемы, которые существуют без данного инструмента.

### Проблема номер 1:
Метаболические модели содержат сотни или тысячи реакций, метаболитов и генов. Ручной анализ таких моделей с помощью таблиц Excel или текстовых файлов невозможен из-за масштаба сложности. Каждая реакция имеет стехиометрические коэффициенты, ограничения по потокам, GPR-правила (gene-protein-reaction) и другие атрибуты. Без специализированного инструмента работа с такими данными превращается в управление тысячами взаимосвязанных переменных.

### Проблема номер 2:
Основной метод анализа метаболических сетей — Flux Balance Analysis (FBA) — требует решения задач линейного программирования (LP) на основе стехиометрической матрицы. Реализация этой оптимизации вручную требует глубоких знаний математического программирования и создания сложной инфраструктуры для построения матриц, задания целевой функции и обработки ограничений.

### Проблема номер 3:
В метаболических моделях существует множество взаимосвязей между генами, белками и реакциями (GPR). Knockout гена может инактивировать одну или несколько реакций, что влияет на весь метаболизм. Моделирование таких эффектов требует сложной логики для отслеживания зависимостей, что сложно реализовать без специализированной библиотеки .

### Проблема номер 4:
Стандарты обмена моделями, такие как SBML (Systems Biology Markup Language), имеют сложную XML-структуру. Парсинг SBML-файлов вручную требует понимания этого формата и написания большого объема кода для извлечения реакций, метаболитов и их свойств.

### Проблема номер 5:
Анализ метаболических сетей часто требует интеграции с другими инструментами — визуализации химических структур (RDKit), анализа графов (NetworkX, igraph) и статистической обработки данных. Без единой экосистемы эти задачи оказываются изолированными и труднодоступными для автоматизации .

---

## 2. Решения, которые предоставляет библиотека

### Решение проблемы номер 1:
COBRApy предоставляет объектную модель для работы с метаболическими моделями. Вы работаете с классами `Model` (модель), `Reaction` (реакция), `Metabolite` (метаболит) и `Gene` (ген), а библиотека управляет всеми внутренними связями. Это позволяет интуитивно исследовать и модифицировать сложные модели.

```python
import cobra
from cobra.io import load_model

# Загружаем встроенную модель E. coli core metabolism
model = load_model("textbook")

print(f"Реакций: {len(model.reactions)}")      # 95
print(f"Метаболитов: {len(model.metabolites)}") # 72
print(f"Генов: {len(model.genes)}")             # 137
```
**Принцип работы:** Библиотека загружает модель из файла и создает иерархию Python-объектов, где каждый элемент (`Reaction`, `Metabolite`, `Gene`) содержит свои атрибуты и ссылки на связанные элементы .

### Решение проблемы номер 2:
COBRApy интегрируется с решателями линейного программирования (GLPK, CPLEX, Gurobi и др.) и предоставляет простой API для проведения FBA. Вам не нужно строить матрицы или настраивать оптимизатор вручную — достаточно указать целевую функцию и вызвать метод `optimize()`.

```python
# Устанавливаем целевую функцию (рост биомассы)
model.objective = "Biomass_Ecoli_core"

# Проводим FBA
solution = model.optimize()

print(f"Скорость роста: {solution.objective_value:.2f}")
print(f"Статус: {solution.status}")
for reaction in model.reactions:
    if abs(solution.fluxes[reaction.id]) > 1e-6:
        print(f"{reaction.id}: {solution.fluxes[reaction.id]:.2f}")
```
**Принцип работы:** Библиотека автоматически строит стехиометрическую матрицу, добавляет ограничения на потоки (bounds) и передает задачу оптимизации выбранному LP-решателю, возвращая решение в удобном формате .

### Решение проблемы номер 3:
COBRApy автоматически управляет связями между генами и реакциями через GPR-правила. Knockout гена или реакции выполняется одной командой, а библиотека корректно обновляет ограничения для всех затронутых реакций.

```python
# Knockout гена
gene = model.genes.get_by_id("b4025")
gene.knock_out()

# Или knockout реакции
reaction = model.reactions.get_by_id("PGI")
reaction.knock_out()

# Проводим FBA после knockout
solution = model.optimize()
print(f"Скорость роста после knockout: {solution.objective_value:.2f}")
```
**Принцип работы:** При knockout библиотека обновляет верхние и нижние границы всех реакций, связанных с этим геном через GPR-правила, в соответствии с логикой "И"/"ИЛИ" .

### Решение проблемы номер 4:
COBRApy поддерживает загрузку моделей из различных форматов: SBML, JSON, MAT. Библиотека берет на себя всю сложность парсинга этих форматов.

```python
from cobra.io import load_json_model, load_sbml_model, load_matlab_model

# Загрузка из JSON (рекомендуемый формат)
model = load_json_model("iJO1366.json")

# Загрузка из SBML
model = load_sbml_model("e_coli_core.xml")

# Загрузка из MATLAB .mat файла
model = load_matlab_model("ecoli.mat")
```
**Принцип работы:** Для каждого формата реализован специальный парсер, который извлекает данные и создает соответствующие объекты Python, сохраняя все связи и атрибуты модели .

### Решение проблемы номер 5:
COBRApy легко интегрируется с другими библиотеками научного стека Python: `numpy`, `pandas`, `scipy`, а также с инструментами визуализации (RDKit, NetworkX, igraph). Это позволяет проводить комплексный анализ в единой среде.

```python
import networkx as nx
import igraph
from rdkit import Chem

# Преобразование модели в граф NetworkX для топологического анализа
graph = nx.DiGraph()
for reaction in model.reactions:
    for met in reaction.reactants:
        for prod in reaction.products:
            graph.add_edge(met.id, prod.id)

print(f"Вершин графа: {graph.number_of_nodes()}")
print(f"Ребер графа: {graph.number_of_edges()}")
```
**Принцип работы:** Объектная модель COBRApy предоставляет прямой доступ к атрибутам реакций и метаболитов, что позволяет легко преобразовывать их в структуры, понятные другим библиотекам .

---

## 3. Инженерная/визуальная грамотность

В контексте COBRApy инженерная грамотность — это умение выйти за рамки простой загрузки модели и выполнения FBA, создавая воспроизводимые, масштабируемые и документированные пайплайны анализа.

**Объектно-ориентированный подход против императивных скриптов:** Вместо написания линейных скриптов, смешивающих загрузку данных, анализ и визуализацию, грамотный инженер создает модульную структуру с четким разделением: слой загрузки моделей, слой проведения симуляций, слой анализа результатов и слой визуализации. Это делает код переиспользуемым для разных моделей и экспериментов.

**Воспроизводимость:** Метаболическое моделирование часто включает множество параметров (решатель, толерантность, границы потоков). Грамотный подход требует фиксации всех параметров в конфигурационных файлах или явного указания в коде, чтобы результаты могли быть воспроизведены другими исследователями.

**Интеграция с системами контроля версий:** Модели в форматах JSON или SBML могут храниться в Git, что позволяет отслеживать изменения в структуре модели. COBRApy хорошо работает с текстовыми форматами, что облегчает совместную работу и контроль версий.

---

## 4. Применение в реальных проектах

### Пример 1: Flux Balance Analysis для E. coli core модели

**Задача:** Провести FBA для модели E. coli core metabolism, определить скорость роста и потоки через ключевые метаболические пути.

**Решение:**
```python
import cobra
from cobra.io import load_model
from cobra.flux_analysis import flux_variability_analysis

# Загрузка модели
model = load_model("textbook")

# Настройка среды: глюкоза доступна, кислород ограничен
model.reactions.EX_glc__D_e.lower_bound = -10.0  # глюкоза поступает
model.reactions.EX_o2_e.lower_bound = 0.0        # кислород отсутствует (анаэробные условия)

# Установка цели — максимизация биомассы
model.objective = "Biomass_Ecoli_core"

# Оптимизация
solution = model.optimize()

print("=== Результаты FBA ===")
print(f"Скорость роста: {solution.objective_value:.4f} 1/ч")

# Анализ вариабельности потоков
fva_result = flux_variability_analysis(model, model.reactions[:10])
print("\nВариабельность потоков для первых 10 реакций:")
print(fva_result)
```
**Результат:** Получена оптимальная скорость роста в анаэробных условиях, идентифицированы ключевые реакции с ненулевыми потоками.

### Пример 2: Анализ knockout-мутаций

**Задача:** Исследовать влияние knockout различных генов на рост E. coli, найти летальные мутации.

**Решение:**
```python
import cobra
from cobra.io import load_model
from cobra.flux_analysis import single_gene_deletion

model = load_model("textbook")

# Проводим симуляцию knockout всех генов по отдельности
deletion_results = single_gene_deletion(model, model.genes)

# Анализируем результаты
lethal_genes = []
for result in deletion_results:
    if result.flux < 0.01:  # Считаем рост ниже 0.01 летальным
        lethal_genes.append(result.gene.id)

print(f"Летальные гены: {lethal_genes}")
print(f"Всего проанализировано генов: {len(deletion_results)}")

# Knockout комбинации генов
from cobra.flux_analysis import double_gene_deletion
double_deletions = double_gene_deletion(model, lethal_genes[:2])
```
**Результат:** Идентифицированы гены, критически важные для роста в данных условиях; получена информация о синтетических летальных взаимодействиях.

### Пример 3: Анализ графовой структуры метаболической сети

**Задача:** Преобразовать метаболическую модель в граф и вычислить топологические характеристики сети.

**Решение:**
```python
import cobra
import networkx as nx
import numpy as np
from cobra.io import load_json_model

# Загрузка полной модели E. coli iJO1366
model = load_json_model("iJO1366.json")

# Построение направленного графа (метаболиты -> метаболиты)
graph = nx.DiGraph()

for reaction in model.reactions:
    # Добавляем ребра от каждого реагента к каждому продукту
    for reactant in reaction.reactants:
        for product in reaction.products:
            # Игнорируем внеклеточные метаболиты для упрощения
            if reactant.compartment != 'e' and product.compartment != 'e':
                graph.add_edge(reactant.id, product.id)

print(f"Количество вершин: {graph.number_of_nodes()}")
print(f"Количество ребер: {graph.number_of_edges()}")

# Вычисление характеристик
degrees = [d for n, d in graph.degree()]
print(f"Средняя степень: {np.mean(degrees):.2f}")

# Поиск центральных метаболитов (hub)
centrality = nx.betweenness_centrality(graph, k=100)
top_hubs = sorted(centrality.items(), key=lambda x: x[1], reverse=True)[:5]
print(f"Центральные метаболиты: {[hub[0] for hub in top_hubs]}")
```
**Результат:** Получена графовая модель метаболической сети, идентифицированы ключевые метаболиты (ATP, NADH, глюкоза), играющие центральную роль в метаболизме .

---
# Модуль 2. Основные идеи и механизмы

## 1. Центральные объекты и архитектуры

Основным центральным объектом COBRApy является **Модель (`Model`)**. Она представляет всю метаболическую сеть и служит контейнером для всех компонентов: реакций, метаболитов, генов и групп. Модель содержит информацию о компартментах, целевую функцию (обычно максимизация роста биомассы) и параметры решателя. Создание или загрузка `Model` — это точка входа в библиотеку .

Вторым ключевым объектом является **Реакция (`Reaction`)**. Это математическое представление биохимического превращения. Каждая реакция содержит стехиометрические коэффициенты (через связи с метаболитами), нижние и верхние границы потока (bounds), GPR-правило (gene-protein-reaction) и название. Реакция может быть обратимой или необратимой, что определяется знаками нижней и верхней границ .

Третий фундаментальный объект — **Метаболит (`Metabolite`)**. Это химическое соединение, участвующее в реакциях. Каждый метаболит имеет химическую формулу, заряд, название и принадлежность к определенному компартменту (цитозоль, митохондрии, внеклеточное пространство и т.д.). Метаболиты хранят список реакций, в которых они участвуют .

Четвертый объект — **Ген (`Gene`)**. Он представляет генетическую информацию, кодирующую ферменты. Гены связаны с реакциями через GPR-правила, которые используют логические операторы (AND, OR) для описания того, какие комбинации генов необходимы для функционирования реакции. Например, правило `(b1234 and b5678) or b9012` означает, что реакция активна, если активны оба гена из первого набора ИЛИ активен ген из второго набора .

Архитектурно COBRApy построена по принципу **объектной модели документа с коллекциями типа `DictList`**. Реакции, метаболиты и гены хранятся в специальных списках, которые позволяют обращаться к элементам как по индексу, так и по идентификатору (через метод `get_by_id` или напрямую как атрибут) .

```python
# Доступ по индексу
first_reaction = model.reactions[0]

# Доступ по идентификатору
pgi = model.reactions.get_by_id("PGI")
atp = model.metabolites.atp_c  # прямой доступ как атрибут
```

Важной частью архитектуры является **интеграция с LP-решателями**. COBRApy поддерживает GLPK (встроенный), CPLEX, Gurobi, COIN-OR и другие. Библиотека абстрагирует различия между решателями, предоставляя единый API для построения и решения задач оптимизации.

---

## 2. Ключевые механизмы работы, возможности, принципы использования

### 1) Механизм загрузки моделей
Поддерживаются форматы SBML, JSON и MATLAB (.mat). Рекомендуемый формат — JSON, обеспечивающий быструю загрузку и полную совместимость.

```python
from cobra.io import load_model, load_json_model, load_sbml_model

# Встроенные модели
model = load_model("textbook")      # E. coli core
model = load_model("iJO1366")       # полная E. coli
model = load_model("salmonella")    # Salmonella

# Из файлов
model = load_json_model("path/to/model.json")
model = load_sbml_model("path/to/model.xml")
```

### 2) Механизм работы с реакциями
Реакции можно просматривать, модифицировать, добавлять и удалять. Важное свойство — возможность изменения границ потоков.

```python
# Получение реакции
pgi = model.reactions.get_by_id("PGI")

# Просмотр стехиометрии
print(pgi.reaction)  # "g6p_c <=> f6p_c"

# Изменение границ (делаем необратимой)
pgi.bounds = (0, 1000)

# Добавление метаболита в реакцию
pgi.add_metabolites({model.metabolites.get_by_id("h_c"): -1})
print(pgi.reaction)  # "g6p_c + h_c <=> f6p_c"
```

### 3) Механизм работы с метаболитами
Метаболиты хранят информацию о формуле, заряде и компартменте.

```python
atp = model.metabolites.get_by_id("atp_c")

print(atp.name)          # "ATP"
print(atp.formula)       # "C10H12N5O13P3"
print(atp.charge)        # -4
print(atp.compartment)   # "c" (cytosol)

# Список реакций с участием метаболита
print(len(atp.reactions))  # 13
```

### 4) Механизм работы с генами
Гены связаны с реакциями через GPR-правила. Можно изменять функциональность генов (knockout).

```python
gene = model.genes.get_by_id("b4025")
print(gene.name)           # "pgi"
print(gene.functional)     # True

# Knockout гена
gene.knock_out()
print(gene.functional)     # False

# Восстановление
gene.functional = True
```

### 5) Механизм оптимизации (FBA)
Основной метод анализа — Flux Balance Analysis. Решается задача линейного программирования для нахождения распределения потоков, максимизирующего целевую функцию.

```python
# Установка цели
model.objective = "Biomass_Ecoli_core"

# Оптимизация
solution = model.optimize()

print(solution.objective_value)  # значение целевой функции
print(solution.status)           # "optimal", "infeasible" и т.д.
print(solution.fluxes)           # Series с потоками для всех реакций
```

### 6) Механизм вариабельности потоков (FVA)
Позволяет определить диапазон возможных значений потока для каждой реакции при сохранении оптимального значения целевой функции.

```python
from cobra.flux_analysis import flux_variability_analysis

fva_result = flux_variability_analysis(model, model.reactions)
print(fva_result.head())
```

### 7) Механизм анализа делеций
Позволяет моделировать удаление генов или реакций и оценивать их влияние на рост.

```python
from cobra.flux_analysis import single_gene_deletion, single_reaction_deletion

# Knockout всех генов
gene_deletions = single_gene_deletion(model, model.genes)

# Knockout всех реакций
reaction_deletions = single_reaction_deletion(model, model.reactions)

# Поиск летальных делеций
lethal = [r for r in reaction_deletions if r.flux < 0.01]
```

### 8) Механизм добавления новых реакций
Можно расширять модель, добавляя новые реакции и метаболиты.

```python
from cobra import Reaction, Metabolite

# Создание нового метаболита
x_met = Metabolite("x_c", formula="C6H12O6", name="Compound X", compartment="c")

# Создание новой реакции
new_reaction = Reaction("RX")
new_reaction.name = "Transport of X"
new_reaction.lower_bound = 0
new_reaction.upper_bound = 1000
new_reaction.add_metabolites({x_met: -1})

# Добавление в модель
model.add_reactions([new_reaction])
```

### 9) Механизм анализа производства (Production Envelope)
Позволяет исследовать зависимость между продукцией двух метаболитов.

```python
from cobra.flux_analysis import production_envelope

# Анализ продукции сукцината и ацетата
prod_env = production_envelope(
    model, 
    reactions=["EX_succ_e", "EX_ac_e"], 
    objective="Biomass_Ecoli_core"
)
```

### 10) Механизм выборки (Sampling)
Позволяет генерировать множество возможных распределений потоков внутри допустимого пространства решений.

```python
from cobra.flux_analysis import sample

samples = sample(model, n=1000, method="achr")
print(samples.head())
```

---

### 7. Работа со структурированными данными

COBRApy превращает сложные биологические модели в управляемые структуры данных Python. Модель — это контейнер, содержащий три основные коллекции: реакции, метаболиты и гены.

Для экспорта данных в форматы, удобные для анализа, используются методы преобразования в `pandas.DataFrame`:

```python
import pandas as pd

# Получение таблицы реакций
reactions_df = pd.DataFrame({
    'id': [r.id for r in model.reactions],
    'name': [r.name for r in model.reactions],
    'lower_bound': [r.lower_bound for r in model.reactions],
    'upper_bound': [r.upper_bound for r in model.reactions],
    'reversibility': [r.reversibility for r in model.reactions]
})

# Получение таблицы метаболитов
metabolites_df = pd.DataFrame({
    'id': [m.id for m in model.metabolites],
    'formula': [m.formula for m in model.metabolites],
    'charge': [m.charge for m in model.metabolites],
    'compartment': [m.compartment for m in model.metabolites]
})
```

Для импорта данных из внешних источников (например, экспериментальных данных) можно динамически изменять границы реакций, отражая условия среды:

```python
# Изменение условий среды на основе экспериментальных данных
conditions = {
    "EX_glc__D_e": -10.0,   # лимит по глюкозе
    "EX_o2_e": -20.0,       # лимит по кислороду
    "EX_nh4_e": -100.0,     # аммоний в избытке
}

for reaction_id, bound in conditions.items():
    model.reactions.get_by_id(reaction_id).lower_bound = bound
```

---

### 8. Итеративные элементы

**Итерация по коллекциям:** Все коллекции (`model.reactions`, `model.metabolites`, `model.genes`) поддерживают итерацию, что позволяет проходить по всем элементам модели.

```python
for reaction in model.reactions:
    if reaction.upper_bound > 1000:
        print(f"Высокий поток: {reaction.id}")
```

**Автоматическое обновление связанных объектов:** При изменении реакции (например, добавлении метаболита) автоматически обновляются списки реакций у соответствующих метаболитов.

**Автоматическая проверка масс-баланса:** Метод `check_mass_balance()` позволяет проверить, сбалансирована ли реакция по массе и заряду.

```python
for reaction in model.reactions:
    imbalance = reaction.check_mass_balance()
    if imbalance:
        print(f"{reaction.id}: {imbalance}")
```

---

### 9. Обработка и отладка ошибок

**Основные типы исключений:**

*   `cobra.exceptions.Infeasible`: Возникает при попытке оптимизации, если задача не имеет решения (например, противоречивые ограничения).
*   `ValueError`: Возникает при установке некорректных границ (нижняя граница больше верхней).
*   `KeyError`: Возникает при обращении к несуществующему идентификатору в DictList.

**Отладка несовместных моделей:** Если модель не имеет решения, полезно использовать функцию `check_consistency()` для выявления проблем:

```python
from cobra.util import check_consistency

inconsistencies = check_consistency(model, rtol=1e-6, atol=1e-9)
print(inconsistencies)
```

**Отладка GPR-правил:** GPR-правила могут быть сложными. Для проверки можно использовать метод `genes` реакции:

```python
reaction = model.reactions.get_by_id("PGI")
print(reaction.gene_reaction_rule)  # "b4025"
print([g.id for g in reaction.genes])  # список связанных генов
```

---

## 10. Композиция и организация результата

**Композиция через объектную модель:** Результат анализа может быть представлен как единый объект решения (`Solution`), содержащий значение целевой функции, статус и потоки для всех реакций.

**Сохранение модифицированных моделей:** После анализа модель можно сохранить для дальнейшего использования.

```python
from cobra.io import save_json_model

# Сохранение в JSON
save_json_model(model, "modified_model.json")
```

**Визуализация результатов:** Результаты FBA можно визуализировать с помощью библиотек matplotlib и seaborn.

```python
import matplotlib.pyplot as plt
import seaborn as sns

# Визуализация потоков
fluxes = solution.fluxes
fluxes_nonzero = fluxes[abs(fluxes) > 1e-6].sort_values()

plt.figure(figsize=(10, 6))
fluxes_nonzero.plot(kind='barh')
plt.xlabel("Поток (ммоль/г/ч)")
plt.title("Распределение потоков в модели")
plt.tight_layout()
plt.show()
```

**Интеграция с RDKit для визуализации химии:** Метаболиты и реакции можно визуализировать в виде химических структур с помощью RDKit .

```python
from rdkit import Chem
from rdkit.Chem import Draw

# Создание SMILES для метаболита
met_smiles = "N[C@@H](C)C(=O)O"  # L-Alanine
mol = Chem.MolFromSmiles(met_smiles)
Draw.MolToImage(mol)
```

---
# Модуль 3. Практическое применение

## 11. Формулировка практического задания

**Контекст:** В системной биологии и метаболической инженерии COBRApy является стандартным инструментом для анализа метаболических сетей. Исследователи используют его для предсказания эффектов генетических модификаций, оптимизации продуктивности штаммов и понимания метаболических адаптаций.

**Прикладная задача:** Разработать скрипт "Анализ метаболической сети E. coli", который демонстрирует возможности COBRApy для:
1.  Загрузки модели и исследования её структуры.
2.  Проведения Flux Balance Analysis в различных условиях среды.
3.  Анализа эффектов knockout генов.
4.  Визуализации результатов.

```python
import cobra
from cobra.io import load_model
from cobra.flux_analysis import flux_variability_analysis, single_gene_deletion
import pandas as pd
import matplotlib.pyplot as plt

# 1. Загрузка модели
print("=== Загрузка модели ===")
model = load_model("textbook")
print(f"Модель: {model.id}")
print(f"Реакций: {len(model.reactions)}")
print(f"Метаболитов: {len(model.metabolites)}")
print(f"Генов: {len(model.genes)}")

# 2. FBA в аэробных условиях
print("\n=== FBA в аэробных условиях ===")
model.reactions.EX_glc__D_e.lower_bound = -10.0  # глюкоза
model.reactions.EX_o2_e.lower_bound = -20.0      # кислород
model.objective = "Biomass_Ecoli_core"

solution = model.optimize()
print(f"Скорость роста (аэробные): {solution.objective_value:.4f} 1/ч")

# 3. FBA в анаэробных условиях
print("\n=== FBA в анаэробных условиях ===")
model.reactions.EX_o2_e.lower_bound = 0.0  # без кислорода
solution_anaerobic = model.optimize()
print(f"Скорость роста (анаэробные): {solution_anaerobic.objective_value:.4f} 1/ч")

# 4. Анализ вариабельности потоков для ключевых реакций
print("\n=== Flux Variability Analysis ===")
key_reactions = ["PGI", "PFK", "FBA", "GAPD", "PGK", "ENO", "PYK"]
fva_result = flux_variability_analysis(model, [model.reactions.get_by_id(r) for r in key_reactions])
print(fva_result)

# 5. Анализ knockout генов
print("\n=== Анализ knockout генов ===")
deletion_results = single_gene_deletion(model, list(model.genes)[:10])  # первые 10 генов
for result in deletion_results:
    if result.flux < 0.01:
        print(f"Летальный: {result.gene.id}")
    else:
        print(f"{result.gene.id}: {result.flux:.4f}")

# 6. Сохранение результатов
print("\n=== Сохранение результатов ===")
fluxes_df = pd.DataFrame(solution.fluxes.items(), columns=["Reaction", "Flux"])
fluxes_df.to_csv("flux_results.csv", index=False)
print("Результаты сохранены в flux_results.csv")

# 7. Визуализация
print("\n=== Визуализация ===")
non_zero_fluxes = solution.fluxes[abs(solution.fluxes) > 0.01].sort_values()
plt.figure(figsize=(12, 6))
non_zero_fluxes.tail(20).plot(kind='barh')
plt.xlabel("Поток (ммоль/г/ч)")
plt.title("Топ-20 реакций по величине потока")
plt.tight_layout()
plt.savefig("flux_distribution.png")
print("График сохранен в flux_distribution.png")
```
**Результат:** Получен полный анализ метаболической сети с выявлением различий между аэробными и анаэробными условиями, идентификацией летальных генов и визуализацией распределения потоков.

---

## 12. Архитектура решения

Решение строится как модульное приложение, демонстрирующее различные сценарии анализа метаболических сетей.

**Структура проекта:**

1.  **Модели данных:** Объекты COBRApy (`Model`, `Reaction`, `Metabolite`, `Gene`) используются напрямую.
2.  **Анализаторы (Analyzers):** Модули, реализующие различные методы анализа (FBA, FVA, делеции, sampling).
3.  **Утилиты визуализации:** Функции для построения графиков потоков, тепловых карт и визуализации метаболических путей.
4.  **Слой ввода/вывода:** Функции для загрузки моделей из различных форматов и сохранения результатов.

```python
# --- analyzers.py ---
from cobra import Model
from cobra.flux_analysis import flux_variability_analysis, single_gene_deletion
import pandas as pd

class FBAAnalyzer:
    """Класс для проведения FBA-анализа"""
    
    def __init__(self, model: Model):
        self.model = model
    
    def optimize_growth(self, glucose_rate: float = -10.0, oxygen_rate: float = None):
        """Оптимизация роста в заданных условиях"""
        self.model.reactions.EX_glc__D_e.lower_bound = glucose_rate
        if oxygen_rate is not None:
            self.model.reactions.EX_o2_e.lower_bound = oxygen_rate
        self.model.objective = "Biomass_Ecoli_core"
        return self.model.optimize()
    
    def analyze_deletions(self, genes=None, reactions=None):
        """Анализ делеций генов или реакций"""
        if genes:
            return single_gene_deletion(self.model, genes)
        elif reactions:
            from cobra.flux_analysis import single_reaction_deletion
            return single_reaction_deletion(self.model, reactions)
        else:
            return None

# --- visualizers.py ---
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def plot_flux_distribution(solution, top_n=20, filename="fluxes.png"):
    """Визуализация распределения потоков"""
    non_zero = solution.fluxes[abs(solution.fluxes) > 0.01].sort_values()
    top_fluxes = non_zero.tail(top_n) if top_n else non_zero
    
    plt.figure(figsize=(10, max(6, top_n * 0.3)))
    top_fluxes.plot(kind='barh', color='steelblue')
    plt.xlabel("Поток (ммоль/г/ч)")
    plt.title(f"Топ-{top_n} реакций по величине потока")
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

def plot_growth_comparison(results, labels, filename="growth_comparison.png"):
    """Сравнение скоростей роста при разных условиях"""
    plt.figure(figsize=(8, 5))
    plt.bar(labels, results, color='steelblue')
    plt.ylabel("Скорость роста (1/ч)")
    plt.title("Сравнение скоростей роста")
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

# --- main.py ---
from cobra.io import load_model
import analyzers
import visualizers

def main():
    # Загрузка модели
    model = load_model("textbook")
    
    # Создание анализатора
    analyzer = analyzers.FBAAnalyzer(model)
    
    # Анализ в разных условиях
    aerobic = analyzer.optimize_growth(glucose_rate=-10.0, oxygen_rate=-20.0)
    anaerobic = analyzer.optimize_growth(glucose_rate=-10.0, oxygen_rate=0.0)
    
    # Визуализация сравнения
    visualizers.plot_growth_comparison(
        [aerobic.objective_value, anaerobic.objective_value],
        ["Аэробные", "Анаэробные"],
        "growth_comparison.png"
    )
    
    # Визуализация распределения потоков
    visualizers.plot_flux_distribution(aerobic, top_n=15, filename="aerobic_fluxes.png")
    visualizers.plot_flux_distribution(anaerobic, top_n=15, filename="anaerobic_fluxes.png")
    
    print("Анализ завершен. Результаты сохранены в изображениях.")

if __name__ == "__main__":
    main()
```

**Минимальный рабочий пример (все вместе)**
```python
from cobra.io import load_model

# Загрузка модели
model = load_model("textbook")

# Установка условий: глюкоза как единственный источник углерода
model.reactions.EX_glc__D_e.lower_bound = -10.0
model.reactions.EX_o2_e.lower_bound = -20.0

# Целевая функция — рост
model.objective = "Biomass_Ecoli_core"

# Оптимизация
solution = model.optimize()

print(f"Скорость роста: {solution.objective_value:.4f} 1/ч")
print(f"Поток через PGI: {solution.fluxes['PGI']:.4f}")
```
**Вывод:** Получается численное значение скорости роста и поток через ключевую реакцию гликолиза.

---

## 13. Этапы реализации

**Этап 1: Загрузка и исследование модели**
Научиться загружать встроенные модели и исследовать их структуру: количество реакций, метаболитов, генов, доступ к отдельным элементам.

**Этап 2: Настройка среды и проведение FBA**
Освоить изменение границ транспортных реакций для симуляции различных условий среды. Провести FBA и интерпретировать результаты.

**Этап 3: Анализ делеций**
Изучить методы анализа knockout генов и реакций. Идентифицировать летальные и нелетальные делеции.

**Этап 4: Расширенные методы анализа**
Освоить flux variability analysis, sampling, production envelope analysis. Понять, как эти методы дополняют FBA.

**Этап 5: Визуализация и сохранение**
Научиться экспортировать результаты в pandas DataFrame, сохранять их в CSV, строить графики потоков и сравнивать условия.

---

## 14. Возникшие сложности

**Сложности работы с COBRApy**

### 1. Понимание единиц измерения и масштабов
Потоки в COBRApy обычно выражаются в ммоль/г/ч (миллимоль на грамм биомассы в час). Начинающие могут неверно интерпретировать значения или задавать некорректные границы.
```python
# Границы должны быть разумными
reaction.lower_bound = -1000.0  # слишком высокое значение может привести к численной нестабильности
```

### 2. Проблемы с решателями
Не все решатели одинаково хорошо работают со всеми моделями. GLPK (встроенный) может быть медленным для больших моделей. Для iJO1366 (2583 реакции) рекомендуется использовать коммерческие решатели (CPLEX, Gurobi).
```python
# Выбор решателя
model.solver = "cplex"  # или "gurobi", "glpk"
```

### 3. Infeasible решения
Если модель не имеет решения (infeasible), это может быть вызвано противоречивыми ограничениями. Необходимо систематически проверять согласованность.
```python
from cobra.util import check_consistency

issues = check_consistency(model)
if issues:
    print(f"Проблемы: {issues}")
```

### 4. Работа с большими моделями
Модель iJO1366 содержит 2583 реакции и 1805 метаболитов. Загрузка и оптимизация такой модели требует больше времени и памяти, чем у core-модели.
```python
# Для больших моделей может потребоваться настройка толерантности
model.solver.configuration.tolerance_feasibility = 1e-7
```

**Ограничения проекта**

### 1. Зависимость от решателей
Для работы с большими моделями требуются коммерческие решатели (CPLEX, Gurobi), которые могут быть недоступны в open-source среде. GLPK бесплатен, но медленнее.

### 2. Отсутствие визуализации по умолчанию
COBRApy не включает встроенные инструменты визуализации. Пользователям необходимо самостоятельно использовать matplotlib, seaborn или специализированные инструменты (escher).

### 3. Ограниченная поддержка кинетических моделей
COBRApy ориентирована на анализ на основе ограничений (steady-state). Динамическое моделирование (кинетика) не поддерживается.

### 4. Сложность добавления новых реакций
Добавление новой реакции требует создания метаболитов с правильными формулами, зарядами и компартментами. Ошибки могут привести к нарушению масс-баланса.

---

## 15. Итоговая оценка инструмента

COBRApy — это **стандарт де-факто** для анализа метаболических сетей на основе ограничений в экосистеме Python. Она предоставляет мощный и гибкий инструментарий для исследования метаболизма организмов от бактерий до человека.

**Ключевые преимущества:**
- **Полнота функционала:** Поддерживает все основные методы анализа метаболических сетей: FBA, FVA, делеции, sampling, gap-filling.
- **Интеграция с научным стеком Python:** Бесшовная работа с pandas, numpy, matplotlib, RDKit, NetworkX .
- **Поддержка форматов:** Чтение и запись моделей в SBML, JSON, MATLAB.
- **Активное сообщество:** Регулярные обновления, документация и поддержка со стороны OpenCOBRA .

**Когда использовать:**
- Анализ метаболических сетей микроорганизмов.
- Предсказание эффектов генетических модификаций для метаболической инженерии.
- Исследование метаболических адаптаций к различным условиям среды.
- Образовательные цели для изучения метаболизма.

**Сравнение с альтернативами:**

| Инструмент | Сильные стороны | Слабые стороны |
| :--- | :--- | :--- |
| **COBRApy** | Интеграция с Python-экосистемой, открытый исходный код, поддержка SBML/JSON | Требует установки решателей, нет встроенной визуализации |
| **MATLAB COBRA Toolbox** | Богатый функционал, проверенный временем | Проприетарный (требует MATLAB), сложность интеграции с другими инструментами |
| **COMETS** | Динамическое моделирование пространственных систем | Узкая специализация, менее универсален |
| **Scipy/CVXOPT** | Полный контроль над оптимизацией | Требует ручного построения матриц и настройки решателей |



