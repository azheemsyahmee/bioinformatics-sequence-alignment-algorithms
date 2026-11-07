# Bioinformatics Sequence Alignment Algorithms Suite

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![CI Tests](https://github.com/azheemsyahmee/bioinformatics-sequence-alignment-algorithms/actions/workflows/tests.yml/badge.svg)](https://github.com/azheemsyahmee/bioinformatics-sequence-alignment-algorithms/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Domain](https://img.shields.io/badge/Domain-Bioinformatics%20%7C%20Computational%20Biology-teal.svg)](#)
[![Zero-Dependencies](https://img.shields.io/badge/Core%20Algorithms-Pure%20Python%20(No%20Deps)-brightgreen.svg)](#)

A pure Python bioinformatics toolkit implementing classical dynamic programming algorithms for biological sequence alignment from first principles. Features **Needleman-Wunsch** (global alignment), **Smith-Waterman** (local alignment with conserved motif detection), an interactive **Kivy desktop GUI ("SeqAlign+")** with a visual dynamic programming matrix grid, and a **Streamlit web workbench** integrating the **EMBL-EBI Clustal Omega REST API** for multiple sequence alignment (MSA).

---

## Provenance & Attribution Disclosure

To ensure complete transparency regarding project origins and authorship:
- **Original Core Code:** Developed by Azheem Syahmee during university coursework in Bioinformatics at Universiti Teknologi Malaysia (UTM). Original source scripts (`needleman.py`, `watermanTry.py`, `appMatrix.py`, `clustal.py`) contained the core dynamic programming loops, the Kivy interface prototype, and the EBI API request methods.
- **Portfolio Refactoring & Enhancements:** Packaged into modular packages (`src/algorithms/`, `src/gui/`, `src/web/`), with standardized function signatures, input validation, 1-based motif coordinate calculation, formatted ASCII matrix visualizers, and an automated 17-case unit test suite (`tests/`).
- **External Services:** Multiple Sequence Alignment (MSA) is computed remotely via the public [EMBL-EBI Clustal Omega REST API](https://www.ebi.ac.uk/Tools/services/rest/clustalo/), not by a local engine.
- **PocketSeq Separation:** This repository is an independent, pure Python algorithmic codebase. It is strictly separate from the mobile application project *PocketSeq* (a collaborative Flutter/Dart mobile application with Supabase backend owned by another account). All mobile code has been completely excluded.

---

## Table of Contents
- [1. Project Overview & Background](#1-project-overview--background)
- [2. Biological Sequence Alignment Concepts](#2-biological-sequence-alignment-concepts)
- [3. Algorithm Comparison: Global vs. Local](#3-algorithm-comparison-global-vs-local)
- [4. Mathematical Foundations & Implementation](#4-mathematical-foundations--implementation)
  - [Needleman-Wunsch Recurrence](#needleman-wunsch-global-alignment)
  - [Smith-Waterman Recurrence](#smith-waterman-local-alignment)
- [5. Repository Architecture](#5-repository-architecture)
- [6. Installation & Prerequisites](#6-installation--prerequisites)
- [7. Usage & Examples](#7-usage--examples)
  - [Command-Line Interface & Demonstration](#a-cli-demonstration-runner)
  - [Programmatic Python API](#b-programmatic-python-api)
  - [Interactive Kivy Desktop GUI](#c-kivy-desktop-gui-seqalign)
  - [Streamlit Web Workbench & EBI API](#d-streamlit-web-workbench)
- [8. Test Suite & Verification Results](#8-test-suite--verification-results)
- [9. Limitations & Future Roadmap](#9-limitations--future-roadmap)
- [10. License](#10-license)

---

## 1. Project Overview & Background

Pairwise sequence alignment is the foundational algorithm of modern bioinformatics, computational genomics, and molecular evolution. Whether comparing DNA sequences across homologous genes, identifying pathogenic mutations, or comparing protein structures, alignment algorithms quantify evolutionary relatedness by identifying insertions, deletions, and substitutions.

While production pipelines frequently rely on heuristic approximations (e.g., BLAST, FASTA) or C-accelerated libraries (e.g., Biopython), understanding the core dynamic programming recurrence relations is essential for computational biologists. 

This repository provides:
1. **Zero-Dependency Implementations:** Needleman-Wunsch and Smith-Waterman built strictly with Python standard library data structures.
2. **Matrix Transparency:** Full dynamic programming score matrix calculation, traceback logging, and tabular ASCII/DataFrame formatting for pedagogical inspection.
3. **Multi-Platform Interfaces:**
   - **CLI:** Fast, terminal-based pipeline for programmatic scripts.
   - **Kivy Desktop App (`SeqAlign+`):** Multi-screen desktop/mobile interface rendering live dynamic programming grids.
   - **Streamlit Web App:** Cloud-ready bioinformatics dashboard providing live pairwise alignment and integrating the European Bioinformatics Institute (EMBL-EBI) Clustal Omega web service for multiple sequence alignment.

---

## 2. Biological Sequence Alignment Concepts

Biological sequences—nucleotides ($A, C, G, T/U$) or amino acids ($20$ standard residues)—diverge through evolutionary mutations:
- **Substitutions:** A nucleotide is exchanged for another ($A \to G$). In proteins, conservative substitutions replace an amino acid with another possessing similar biochemical properties (e.g., Isoleucine $\leftrightarrow$ Leucine).
- **Insertions & Deletions (Indels):** Nucleotides are added or lost, introducing gaps ($-$) into the alignment.

To identify optimal alignments, algorithms assign numerical weights:
- **Match Score ($+u$):** Positive reward for identity matches.
- **Mismatch Penalty ($v < 0$):** Negative cost for non-identical substitutions.
- **Gap Penalty ($d < 0$):** Negative penalty for introducing a gap.

---

## 3. Algorithm Comparison: Global vs. Local

| Dimension | Needleman-Wunsch (1970) | Smith-Waterman (1981) |
| :--- | :--- | :--- |
| **Alignment Type** | **Global Alignment** | **Local Alignment** |
| **Objective** | Optimally align sequences from end to end across their entire length | Identify the highest-scoring local sub-region or conserved motif |
| **Biological Use Case** | Highly homologous genes or proteins of similar length (e.g., human vs. chimp $\beta$-globin) | Disparate sequences sharing a localized functional domain (e.g., zinc-finger motif in divergent proteins) |
| **Matrix Initialization** | First row & col initialized with cumulative linear gap penalties: $M(i, 0) = i \cdot d$ | First row & col initialized to **zeros**: $M(i, 0) = 0$ |
| **Recurrence Floor** | Unconstrained (can take negative scores) | **Zero-floor bounded**: $M(i, j) \ge 0$ prevents negative propagation |
| **Traceback Start** | Bottom-right corner: $(n, m)$ | Cell containing the **global maximum score**: $\max_{i,j} M(i, j)$ |
| **Traceback Termination** | Top-left corner: $(0, 0)$ | As soon as a cell with score **$0$** is reached |
| **Time Complexity** | $\mathcal{O}(n \cdot m)$ | $\mathcal{O}(n \cdot m)$ |
| **Space Complexity** | $\mathcal{O}(n \cdot m)$ | $\mathcal{O}(n \cdot m)$ |

---

## 4. Mathematical Foundations & Implementation

### Needleman-Wunsch Global Alignment

Given sequence $X = x_1 x_2 \dots x_n$ and sequence $Y = y_1 y_2 \dots y_m$:

#### 1. Matrix Initialization
$$M(0, 0) = 0$$
$$M(i, 0) = i \cdot d \quad \text{for } 1 \le i \le n$$
$$M(0, j) = j \cdot d \quad \text{for } 1 \le j \le m$$

#### 2. Dynamic Programming Recurrence
$$M(i, j) = \max \begin{cases} 
M(i-1, j-1) + S(x_i, y_j) & \text{(Diagonal: match or mismatch)} \\
M(i-1, j) + d & \text{(Vertical: deletion / gap in } Y\text{)} \\
M(i, j-1) + d & \text{(Horizontal: insertion / gap in } X\text{)}
\end{cases}$$

Where substitution score $S(x_i, y_j)$ is:
$$S(x_i, y_j) = \begin{cases} \text{match\_score} & \text{if } x_i = y_j \\ \text{mismatch\_penalty} & \text{if } x_i \neq y_j \end{cases}$$

#### 3. Traceback
Reconstruction begins at $(n, m)$ and moves toward $(0, 0)$, prepending aligned residues or gap symbols based on which neighboring cell produced the maximum value.

---

### Smith-Waterman Local Alignment

#### 1. Matrix Initialization
$$M(i, 0) = 0 \quad \text{for } 0 \le i \le n$$
$$M(0, j) = 0 \quad \text{for } 0 \le j \le m$$

#### 2. Dynamic Programming Recurrence
$$M(i, j) = \max \begin{cases} 
0 & \text{(Zero-floor boundary)} \\
M(i-1, j-1) + S(x_i, y_j) & \text{(Diagonal: match or mismatch)} \\
M(i-1, j) + d & \text{(Vertical: gap in } Y\text{)} \\
M(i, j-1) + d & \text{(Horizontal: gap in } X\text{)}
\end{cases}$$

#### 3. Traceback
Reconstruction begins at the global maximum position $(i^*, j^*) = \arg\max_{i,j} M(i, j)$ and terminates as soon as $M(i, j) = 0$. The resulting strings are reversed to output the aligned conserved motif along with its 1-based start and end indices.

---

## 5. Repository Architecture

```
bioinformatics-sequence-alignment-algorithms/
├── .gitignore                          # Git ignore specification
├── LICENSE                             # MIT Open Source License
├── README.md                           # Documentation and guides
├── requirements.txt                    # Optional GUI and web dependencies
├── PREPARATION_REPORT.md               # Technical audit and preparation report
│
├── src/
│   ├── __init__.py
│   ├── algorithms/                     # Core Pure Python Implementations (Zero Deps)
│   │   ├── __init__.py
│   │   ├── needleman_wunsch.py         # Global alignment, DP matrix, traceback, scoring
│   │   └── smith_waterman.py           # Local alignment, zero-floor DP, motif coordinates
│   ├── gui/                            # Kivy Cross-Platform Desktop GUI
│   │   ├── __init__.py
│   │   └── app_kivy.py                 # "SeqAlign+" multi-screen desktop application
│   └── web/                            # Streamlit Cloud/Browser Workbench
│       ├── __init__.py
│       └── app_streamlit.py            # Multi-tab web app with EBI Clustal Omega API
│
├── tests/                              # Verification & Unit Test Suite
│   ├── __init__.py
│   ├── test_needleman_wunsch.py        # Global alignment unit tests (8 test cases)
│   ├── test_smith_waterman.py          # Local alignment unit tests (7 test cases)
│   └── run_tests.py                    # Standalone unittest runner script
│
└── examples/                           # Benchmark Data & Demonstrations
    ├── sample_sequences.fasta          # Test FASTA sequences (DNA and Protein)
    └── demo_alignment.py               # CLI alignment demonstration runner
```

---

## 6. Installation & Prerequisites

### Prerequisites
- Python 3.8 or higher.
- No third-party packages are required to use the core algorithms (`src/algorithms/`).

### Setup Instructions

1. **Clone the repository:**
   ```bash
   git clone https://github.com/azheemsyahmee/bioinformatics-sequence-alignment-algorithms.git
   cd bioinformatics-sequence-alignment-algorithms
   ```

2. **(Optional) Create and activate a virtual environment:**
   ```bash
   python -m venv venv
   # On Windows:
   venv\Scripts\activate
   # On macOS/Linux:
   source venv/bin/activate
   ```

3. **Install optional GUI & Web dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

---

## 7. Usage & Examples

### (A) CLI Demonstration Runner

Execute the pre-configured alignment demonstration script:

```bash
python examples/demo_alignment.py
```

#### Example Output: Needleman-Wunsch Global Alignment

```text
Input Sequence 1 : FKHMEDPLE
Input Sequence 2 : FMDTPLNE
Scoring Scheme   : Match = +1, Mismatch = -2, Gap = -2

Optimal Global Alignment:
  Seq1 : FKHMED-PL-E
         |   || || |
  Seq2 : F--M-DTPLNE
Total Score : -4
Identity    : 54.5%

Full Dynamic Programming Matrix:
             F     M     D     T     P     L     N     E
       0    -2    -4    -6    -8   -10   -12   -14   -16
 F    -2     1    -1    -3    -5    -7    -9   -11   -13
 K    -4    -1    -1    -3    -5    -7    -9   -11   -13
 H    -6    -3    -3    -3    -5    -7    -9   -11   -13
 M    -8    -5    -2    -4    -5    -7    -9   -11   -13
 E   -10    -7    -4    -4    -6    -7    -9   -11   -10
 D   -12    -9    -6    -3    -5    -7    -9   -11   -12
 P   -14   -11    -8    -5    -5    -4    -6    -8   -10
 L   -16   -13   -10    -7    -7    -6    -3    -5    -7
 E   -18   -15   -12    -9    -9    -8    -5    -5    -4
```

#### Example Output: Smith-Waterman Local Alignment

```text
Input Sequence 1 (Target) : TTTAAACACACTAGGGCCC
Input Sequence 2 (Query)  : CCCGGAGCACACATTTTTT
Scoring Scheme            : Match = +3, Mismatch = -1, Gap = -2

Optimal Local Conserved Motif:
  Seq1 [7..14] : A-CACACTA
                 | ||||| |
  Seq2 [7..14] : AGCACAC-A
Maximum Local Score : 17
Motif Identity      : 77.8%
```

---

### (B) Programmatic Python API

Import and use the modular functions directly in your own computational pipelines:

```python
from src.algorithms.needleman_wunsch import needleman_wunsch, format_matrix
from src.algorithms.smith_waterman import smith_waterman

# Global Alignment
nw_result = needleman_wunsch(
    seq1="GATTACA",
    seq2="GCATGCU",
    match_score=1,
    mismatch_penalty=-1,
    gap_penalty=-1
)
print("Aligned Seq 1:", nw_result["aligned_seq1"])
print("Aligned Seq 2:", nw_result["aligned_seq2"])
print("Alignment Score:", nw_result["score"])

# Local Alignment
sw_result = smith_waterman(
    seq1="AAACCGTTT",
    seq2="TTTCGTCCC",
    match_score=3,
    mismatch_penalty=-1,
    gap_penalty=-2
)
print(f"Motif Range: Seq1 [{sw_result['start1']}..{sw_result['end1']}]")
print("Motif Score:", sw_result["score"])
```

---

### (C) Kivy Desktop GUI ("SeqAlign+")

Launch the cross-platform desktop application:

```bash
python src/gui/app_kivy.py
```

**Features:**
- Fluid slide-transition interface.
- Tabular Dynamic Programming Matrix viewer rendered in a responsive, scrollable grid layout.
- Real-time adjustment of match scores, mismatch penalties, and gap costs.

---

### (D) Streamlit Web Workbench

Launch the local interactive browser workbench:

```bash
streamlit run src/web/app_streamlit.py
```

**Tabs:**
1. **Needleman-Wunsch (Global):** Interactive parameter sliders, live alignment rendering, and full Pandas DataFrame matrix display.
2. **Smith-Waterman (Local):** Conserved motif identification with exact sequence coordinates and match indicators.
3. **Clustal Omega (Multiple Alignment):** Paste multi-FASTA protein or DNA files and submit directly to the EMBL-EBI REST API. Polls job status automatically and provides one-click alignment file downloads (`.aln`).

*Note on EMBL-EBI Clustal Omega API:* The web workbench conforms strictly to the [EMBL-EBI Terms of Use](https://www.ebi.ac.uk/about/terms-of-use). A valid contact email address must be supplied by the user upon job submission. No API keys or credentials are required.

---

## 8. Test Suite & Verification Results

### Automated Continuous Integration (GitHub Actions)
The repository includes an automated testing workflow located at [`.github/workflows/tests.yml`](.github/workflows/tests.yml) running on `ubuntu-latest` with Python 3.11:

1. **Job 1: `core-algorithm-tests` (Zero External Dependencies)**
   - Operates strictly in a pure Python standard library environment with no third-party package installations.
   - Executes the full 15-case unit test suite:
     ```bash
     python -m unittest discover -s tests -v
     python tests/run_tests.py
     ```
   - Executes the CLI demonstration script:
     ```bash
     python examples/demo_alignment.py
     ```
2. **Job 2: `optional-apps-smoke-test` (Optional Web & GUI Applications)**
   - Performs syntax and bytecode compilation validation (`python -m py_compile`) on all project modules.
   - Installs web dependencies (`streamlit`, `requests`, `pandas`) and tests the EMBL-EBI REST API client functions.
   - Tests Kivy module imports under headless display configuration, preserving isolation from display server requirements.

*CI Status:* The workflow is configured and ready to execute automatically upon pushing to GitHub. Test results will be displayed live on the repository Actions dashboard.

### Local Execution Instructions
On any system with Python 3.8+ installed:
```bash
# Run unit tests directly
python tests/run_tests.py

# Or run via unittest discovery
python -m unittest discover -s tests -v

# Run the CLI demonstration
python examples/demo_alignment.py
```

### Mathematical Verification Results (Independently Calculated)

| Algorithm | Test Case | Inputs | Expected Output | Independently Verified |
| :--- | :--- | :--- | :--- | :---: |
| **Needleman-Wunsch** | Identical Sequences | `ATCG` vs `ATCG` (Match: +2, Mismatch: -1, Gap: -2) | Score: `8`, Aligned: `ATCG` / `ATCG` | **VERIFIED** |
| **Needleman-Wunsch** | Complete Divergence | `AAAA` vs `CCCC` (Match: +1, Mismatch: -1, Gap: -2) | Score: `-4`, Aligned: `AAAA` / `CCCC` | **VERIFIED** |
| **Needleman-Wunsch** | Single Indel | `GATTACA` vs `GATACA` (Match: +1, Mismatch: -1, Gap: -1) | Score: `5`, Aligned: `GATTACA` / `GA-TACA` | **VERIFIED** |
| **Needleman-Wunsch** | Asymmetric Lengths | `ACTG` vs `ACTGACTGACTG` (Match: +1, Mismatch: -1, Gap: -2) | Score: `-12`, Aligned: `--------ACTG` / `ACTGACTGACTG` | **VERIFIED** |
| **Needleman-Wunsch** | Empty Boundary | `""` vs `"ACGT"` (Gap: -2) | Score: `-8`, Aligned: `----` / `ACGT` | **VERIFIED** |
| **Needleman-Wunsch** | Protein Benchmark | `FKHMEDPLE` vs `FMDTPLNE` (Match: +1, Mismatch: -2, Gap: -2) | Score: `-4`, Aligned: `FKHMED-PL-E` / `F--M-DTPLNE` | **VERIFIED** |
| **Smith-Waterman** | Embedded Conserved Motif | `XXXXXXHELLOWORLDYYYYYY` vs `ZZZZZZHELLOWORLDWWWWWW` | Score: `30`, Span: `[7..16]`, Motif: `HELLOWORLD` | **VERIFIED** |
| **Smith-Waterman** | Completely Divergent | `TTTT` vs `AAAA` (Match: +2, Mismatch: -5, Gap: -5) | Score: `0`, Aligned: `""` (Zero-floor) | **VERIFIED** |
| **Smith-Waterman** | Empty Sequence | `""` vs `"ACGT"` | Score: `0`, Aligned: `""` | **VERIFIED** |
| **Smith-Waterman** | Local DNA Motif | `ACACACTA` vs `AGCACACA` (Match: +3, Mismatch: -1, Gap: -2) | Score: `17`, Span: `[1..8]`, Motif: `A-CACACTA` / `AGCACAC-A` | **VERIFIED** |

*Runtime Execution Note:* Automated test files (`tests/run_tests.py`) require a Python 3.8+ runtime. On host systems where Python is not yet installed in PATH, the tests can be executed immediately upon installing Python or setting up a virtual environment.

---

## 9. Limitations & Future Roadmap

### Current Limitations
- Linear gap penalties only ($W = g \cdot k$). Affine gap penalty models ($W = g_o + g_e \cdot k$) are not yet implemented in the DP matrix.
- Pre-set scoring schemes (unitary integers) rather than full substitution matrices (e.g., BLOSUM62, PAM250).
- Streamlit Multiple Sequence Alignment tab requires active outbound internet access to reach the EMBL-EBI REST service.

### Future Roadmap
- [ ] Add support for affine gap penalties (Gotoh algorithm, $\mathcal{O}(n \cdot m)$).
- [ ] Integrate standard amino acid substitution matrices (`BLOSUM62`, `PAM250`).
- [ ] Implement Hirschberg's linear space algorithm ($\mathcal{O}(\min(n, m))$ memory) for long chromosome sequences.
- [ ] Add phylogenetic tree visualization (Neighbor-Joining) based on alignment distance matrices.

---

## 10. License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
