# 🧬 Bioinformatics Practicals

## 📋 Overview

This repository contains bioinformatics computational biology practical solutions.

---

## 🧪 Q1: DNA Palindrome Analysis

### 🎯 Objective

Analyze DNA sequences to find palindromic substrings and calculate their GC content.

### 📝 Code Explanation:

#### 1️⃣ Input DNA Sequence

```python
dna_seq = "AGCTTTCGAATTCGA"
```

Starting DNA sequence to analyze.

---

#### 2️⃣ Find All Unique Palindromic Substrings 🔍

```python
def is_dna_palindrome(s):
    complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    return all(complement[s[i]] == s[-(i+1)] for i in range(len(s)))

unique_palindromes = set()
n = len(dna_seq)

for i in range(n):
    for j in range(i+2, n+1):   # Minimum length 2
        substring = dna_seq[i:j]
        if is_dna_palindrome(substring):
            unique_palindromes.add(substring)

unique_palindromes = sorted(unique_palindromes, key=lambda x: (-len(x), x))
```

**What it does:**

- `is_dna_palindrome()`: Checks if a DNA sequence is palindromic by comparing each base with its complement in reverse position
- Uses nested loops to extract all possible substrings (minimum length 2)
- Stores unique palindromes in a set
- Sorts by length (longest first) then alphabetically

---

#### 3️⃣ Calculate GC Content 📊

```python
def gc_content(seq):
    gc_count = seq.count('G') + seq.count('C')
    return round((gc_count / len(seq)) * 100, 2)

for p in unique_palindromes:
    print(f"{p:20} | {gc_content(p):5.2f}%")
```

**What it does:**

- Counts G and C nucleotides in the sequence
- Calculates percentage: `(G+C count / total length) × 100`
- Prints formatted table with palindrome and its GC%

---

#### 4️⃣ Find Longest Palindrome 🏆

```python
longest_length = max(len(p) for p in unique_palindromes)
longest_palindromes = [p for p in unique_palindromes if len(p) == longest_length]

print(f"\n Longest Palindromic SubString(s) (Length: {longest_length}):")
for p in longest_palindromes:
    print(f" {p} - GC: {gc_content(p)}")
```

**What it does:**

- Finds maximum palindrome length
- Filters palindromes matching that length
- Displays longest palindrome(s) with GC content

---

## 📈 Q2: Differential Gene Expression Analysis (DEA)

### 🎯 Objective

Compare gene expression levels between control and treatment samples to identify significantly changed genes.

### 📝 Code Explanation:

#### 1️⃣ Import Libraries & Load Data 📥

```python
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv('dea.csv')
df.head()
```

**What it does:**

- Imports required libraries for data analysis and visualization
- Reads CSV file containing gene expression data
- Displays first 5 rows

---

#### 2️⃣ Separate Control and Treatment Columns 🔀

```python
control = [c for c in df.columns if 'Control' in c]
treatment = [c for c in df.columns if 'Treatment' in c]

print("Control Columns: ", control)
print("Treatment Columns: ", treatment)
```

**What it does:**

- Uses list comprehension to filter column names
- Identifies all columns containing "Control" (25 samples)
- Identifies all columns containing "Treatment" (25 samples)

---

#### 3️⃣ Calculate Mean Expression Levels 🧮

```python
df['Control_Mean'] = df[control].mean(axis=1)
df['Treatment_Mean'] = df[treatment].mean(axis=1)
```

**What it does:**

- `mean(axis=1)`: Calculates row-wise average across all control samples
- Creates new column for control group average expression
- Creates new column for treatment group average expression

---

#### 4️⃣ Calculate Fold Change 📊

```python
epsilon = 1e-6
df['Fold_Change'] = (df['Treatment_Mean'] + epsilon) / (df['Control_Mean'] + epsilon)

df['log2FC'] = np.log2(df['Fold_Change'])
```

**What it does:**

- `epsilon`: Small value added to avoid division by zero
- `Fold_Change`: Ratio of treatment to control (how much gene expression changed)
- `log2FC`: Log2 transformation for better visualization and interpretation
  - log2FC = 1 means 2× increase
  - log2FC = -1 means 2× decrease

---

#### 5️⃣ Identify Up-regulated & Down-regulated Genes 🔺🔻

```python
up_regulated = df[df['Fold_Change'] >= 1.6]
down_regulated = df[df['Fold_Change'] < 0.67]

print(f"Up: {len(up_regulated)}, Down: {len(down_regulated)}")
print("\nTop 5 Up:", up_regulated[['Gene','Fold_Change']].head().values)
print("\nTop 5 Down:", down_regulated[['Gene','Fold_Change']].head().values)

display(up_regulated[["Gene", "Control_Mean", "Treatment_Mean", "Fold_Change", "log2FC"]])
display(down_regulated[["Gene", "Control_Mean", "Treatment_Mean", "Fold_Change", "log2FC"]])
```

**What it does:**

- Filters genes with FC ≥ 1.6 (60%+ increase) → **up-regulated** 📈
- Filters genes with FC < 0.67 (33%+ decrease) → **down-regulated** 📉
- Prints counts and top 5 genes from each category
- Displays full tables with gene details

---

#### 6️⃣ Visualize Fold Changes 📊

```python
plt.figure(figsize=(12, 6))
plt.bar(df.index[:30], df['Fold_Change'].sort_values(ascending=False)[:30])
plt.axhline(y=1.6, color='r', linestyle='--', label='Up-regulation Threshold (1.6)')
plt.axhline(y=0.67, color='b', linestyle='--', label='Down-regulation Threshold (0.67)')
plt.xticks(rotation=90)
plt.ylabel('Fold Change')
plt.legend()
plt.tight_layout()
plt.show()
```

**What it does:**

- Creates bar chart of top 30 genes by fold change
- Adds horizontal red line at 1.6 (up-regulation threshold)
- Adds horizontal blue line at 0.67 (down-regulation threshold)
- Rotates x-axis labels for readability
- Shows which genes changed significantly

### 🔬 Key Interpretations:

- **FC = 1.6** → Gene expression increased by 60%
- **FC = 0.67** → Gene expression decreased by 33%
- **FC = 1.0** → No change in expression
- **FC = 2.0** → Expression doubled
- **FC = 0.5** → Expression halved

---

## 📊 Data Files

- `dea.csv` - Differential expression analysis data (Q2)

---

## 👨‍💻 Author

Bioinformatics Student | IT 4133 Course

---

_Generated for academic purposes_ 📚
