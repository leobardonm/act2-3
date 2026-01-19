# TC2038 - Advanced Algorithms Implementation
## Chapters 2-3: String Algorithms, Graph Algorithms, and Combined System

### Author Information
- **Course:** TC2038 - Analysis and Design of Advanced Algorithms
- **Assignment:** Activity 2-3

---

## 📁 Project Structure

```
act2-3/
├── main.cpp              # Main program with menu system
├── utils.h               # Utilities (benchmarking, data generation)
├── string_algorithms.h   # KMP, Z-Algorithm, Rabin-Karp
├── dna_algorithms.h      # Manacher, Suffix Array, LCS
├── graph_algorithms.h    # Max Flow, MST, 2nd MST, Edge Types, BCC
├── expression_parser.h   # Combined system (Expression Parser)
├── plot_benchmarks.py    # Python script for generating plots
├── README.md             # This file
└── TC2038-31.pdf         # Assignment specification
```

---

## 🔧 Compilation & Execution

# Compilar
g++ -std=c++17 -O2 -Wall -o algorithms main.cpp

# Ejecutar
./algorithms

# Generar gráficas (después de correr benchmarks)
python3 plot_benchmarks.py
---

## 📋 Menu Options

| Option | Description |
|--------|-------------|
| 1 | String Matching Demo (KMP, Z-Algorithm, Rabin-Karp) |
| 2 | DNA Algorithms Demo (Manacher, Suffix Array, LCS) |
| 3 | Security Log Analysis (Practical Application) |
| 4 | Genomics Application (Practical Application) |
| 5 | Graph Algorithms Demo (Max Flow, MST, 2nd MST, etc.) |
| 6 | Expression Parser (Polish/Infix, AST, Truth Tables) |
| 7 | String Algorithms Benchmark |
| 8 | DNA Algorithms Benchmark |
| 9 | Graph Algorithms Benchmark |
| 10 | Run All Benchmarks |
| 0 | Exit |

---

## 📊 Implemented Algorithms

### Problem 1: String Algorithms (Web Server Logs)
- **KMP (Knuth-Morris-Pratt):** O(n+m) - With π-array demonstration for patterns with nontrivial borders
- **Z-Algorithm:** O(n+m) - Pattern matching using P#T concatenation
- **Rabin-Karp:** O(n+m) average - With single and double hashing variants

### Problem 2: DNA Sequence Analysis
- **Manacher's Algorithm:** O(n) - Finding all palindromic substrings with correctness justification
- **Suffix Array + Kasai LCP:** O(n log n) - Efficient motif queries (4-12 bp)
- **Longest Common Substring:** Two implementations:
  - Via Suffix Array: O(n log n)
  - Via Dynamic Programming: O(nm)

### Graph Problems
1. **Maximum Flow:** Edmonds-Karp O(VE²) and Dinic O(V²E)
2. **Minimum Spanning Tree:** Kruskal O(E log E) and Prim O(E log V)
3. **Second MST:** Edge replacement strategy with necessity proof
4. **Edge Classification:** DFS/BFS edge types with mathematical proofs
5. **Biconnected Components:** Tarjan's algorithm O(V+E)

### Combined System: Expression Parser
- Detects Polish (prefix) vs Infix notation using KMP/Z-Algorithm
- Builds Abstract Syntax Tree (AST)
- Evaluates logical expressions with truth tables
- Supports operators: NOT (~), AND (&), OR (|), XOR (^)
- Variables: P, Q, R

---

## 📈 Benchmarking

The program includes empirical benchmarking for all algorithms:
- Tests multiple input sizes
- Averages over 10 trials
- Exports results to CSV files
- Python script generates publication-ready plots

### Output Files
- `string_benchmark.csv` - String algorithm timings
- `dna_benchmark.csv` - DNA algorithm timings
- `graph_benchmark.csv` - Graph algorithm timings
- `*_benchmark_plot.png` - Visualization plots

---

## 🔬 Complexity Summary

| Algorithm | Time | Space | Application |
|-----------|------|-------|-------------|
| KMP | O(n+m) | O(m) | Pattern matching |
| Z-Algorithm | O(n+m) | O(n) | Pattern matching |
| Rabin-Karp | O(n+m) avg | O(m) | Multi-pattern matching |
| Manacher | O(n) | O(n) | Palindrome detection |
| Suffix Array | O(n log n) | O(n) | String indexing |
| LCS (SA) | O(n log n) | O(n) | Sequence comparison |
| LCS (DP) | O(nm) | O(nm) | Sequence comparison |
| Edmonds-Karp | O(VE²) | O(V²) | Network flow |
| Dinic | O(V²E) | O(V²) | Network flow |
| Kruskal | O(E log E) | O(V) | MST |
| Prim | O(E log V) | O(V) | MST |
| Second MST | O(E log E + V²) | O(V+E) | Redundant MST |
| Biconnected | O(V+E) | O(V+E) | Graph connectivity |
| Expression Parser | O(n) | O(n) | Logical evaluation |

---

## 💡 Key Insights

### Cybersecurity Application
The string matching algorithms are essential for:
- Intrusion detection systems (IDS)
- Log analysis for attack patterns
- Malware signature detection
- Real-time network traffic analysis

### Genomics Application
DNA analysis algorithms enable:
- Finding palindromic sequences (restriction enzyme sites)
- Genome indexing for rapid queries
- Sequence alignment and comparison
- Motif discovery in regulatory regions

### Expression Parser Design
The combined system demonstrates:
- Integration of string matching for notation detection
- Recursive descent parsing for AST construction
- Complete truth table generation for verification

---

## 📝 Test Cases (Expression Parser)

### Polish Notation (Prefix)
1. `& P Q` → Valid, AST: AND(P, Q)
2. `| & P Q ~ R` → Valid, AST: OR(AND(P, Q), NOT(R))
3. `& P` → Invalid (missing operand)
4. `P Q` → Invalid (no operator)

### Infix Notation
1. `P & Q` → Valid, AST: AND(P, Q)
2. `(P | Q) & ~R` → Valid, AST: AND(OR(P, Q), NOT(R))
3. `P & & Q` → Invalid (consecutive operators)
4. `(P & Q` → Invalid (unmatched parenthesis)

---

## 📚 References

- Cormen, T. H., et al. "Introduction to Algorithms" (CLRS)
- Gusfield, D. "Algorithms on Strings, Trees, and Sequences"
- Sedgewick, R. "Algorithms in C++"
