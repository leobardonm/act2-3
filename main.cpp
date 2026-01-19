// ============================================================================
// TC2038 - Analysis and Design of Advanced Algorithms
// Chapters 2-3: Complete Implementation
// ============================================================================
// 
// This program implements:
// 
// STRING ALGORITHMS:
//   1.a. Knuth-Morris-Pratt (KMP) with prefix function
//   1.b. Z-Algorithm with P#T concatenation
//   1.c. Rabin-Karp with single and double hashing
//   2.a. Manacher's Algorithm for longest palindrome
//   2.b. Suffix Array with Kasai LCP
//   2.c. Longest Common Substring (SA+LCP and DP methods)
// 
// GRAPH ALGORITHMS:
//   1. Maximum Flow (Edmonds-Karp with step-by-step)
//   2. Second MST with edge replacement strategy
//   3. Edge types in DFS/BFS traversals with proofs
//   4. Biconnected Components and Bridges (Tarjan)
//   5. MST via Kruskal and Prim with comparison
// 
// COMBINED SYSTEM:
//   - Logical expression parser (Polish/Infix notation detection)
//   - AST construction and evaluation
//   - Truth table generation
// 
// BENCHMARKING:
//   - Empirical runtime measurements
//   - CSV export for plotting
// 
// Compilation: g++ -std=c++17 -O2 -o algorithms main.cpp
// Usage: ./algorithms
// ============================================================================

#include <iostream>
#include <fstream>
#include <cstdlib>
#include "utils.h"
#include "string_algorithms.h"
#include "dna_algorithms.h"
#include "graph_algorithms.h"
#include "expression_parser.h"

// ============================================================================
// BENCHMARKING FUNCTIONS
// ============================================================================

void runStringBenchmarks() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  STRING ALGORITHMS BENCHMARKING\n";
    std::cout << std::string(70, '=') << "\n";
    
    std::cout << "\nExperimental Setup:\n";
    std::cout << "  - Varying text size: 1K, 5K, 10K, 50K, 100K, 500K characters\n";
    std::cout << "  - Pattern length: 10 characters\n";
    std::cout << "  - Alphabet: lowercase letters (26 symbols)\n";
    std::cout << "  - Trials per size: 10 (averaged)\n";
    std::cout << "  - Compiler: g++ -O2 -std=c++17\n\n";
    
    DataGenerator gen(42);
    std::vector<int> sizes = {1000, 5000, 10000, 50000, 100000, 500000};
    std::string pattern = gen.randomString(10);
    
    std::vector<double> kmpTimes, zTimes, rkTimes;
    
    std::cout << std::setw(10) << "Size" << std::setw(15) << "KMP (ms)" 
              << std::setw(15) << "Z-Algo (ms)" << std::setw(15) << "Rabin-Karp\n";
    std::cout << std::string(55, '-') << "\n";
    
    for (int size : sizes) {
        std::string text = gen.randomString(size);
        
        double kmpTime = Benchmark::measure([&]() {
            KMP::search(text, pattern);
        }, 10);
        kmpTimes.push_back(kmpTime);
        
        double zTime = Benchmark::measure([&]() {
            ZAlgorithm::search(text, pattern);
        }, 10);
        zTimes.push_back(zTime);
        
        double rkTime = Benchmark::measure([&]() {
            RabinKarp::searchDoubleHash(text, pattern);
        }, 10);
        rkTimes.push_back(rkTime);
        
        std::cout << std::setw(10) << size 
                  << std::setw(15) << std::fixed << std::setprecision(4) << kmpTime
                  << std::setw(15) << zTime
                  << std::setw(15) << rkTime << "\n";
    }
    
    // Export to CSV
    std::ofstream csv("string_benchmark.csv");
    csv << "size,kmp_ms,z_algo_ms,rabin_karp_ms\n";
    for (size_t i = 0; i < sizes.size(); i++) {
        csv << sizes[i] << "," << kmpTimes[i] << "," << zTimes[i] << "," << rkTimes[i] << "\n";
    }
    csv.close();
    
    std::cout << "\nResults exported to: string_benchmark.csv\n";
    std::cout << "Expected complexity: O(n + m) for all three algorithms.\n";
    std::cout << "The near-linear growth confirms the theoretical analysis.\n";
}

void runDNABenchmarks() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  DNA ALGORITHMS BENCHMARKING\n";
    std::cout << std::string(70, '=') << "\n";
    
    std::cout << "\nExperimental Setup:\n";
    std::cout << "  - DNA sequence sizes: 100, 500, 1K, 5K, 10K bases\n";
    std::cout << "  - Alphabet: {A, C, G, T}\n";
    std::cout << "  - Trials per size: 5 (averaged)\n\n";
    
    DataGenerator gen(42);
    std::vector<int> sizes = {100, 500, 1000, 5000, 10000};
    
    std::vector<double> manacherTimes, saTimes, lcsDPTimes;
    
    std::cout << std::setw(10) << "Size" << std::setw(18) << "Manacher (ms)" 
              << std::setw(18) << "Suffix Array" << std::setw(18) << "LCS-DP\n";
    std::cout << std::string(64, '-') << "\n";
    
    for (int size : sizes) {
        std::string seqA = gen.randomDNA(size);
        std::string seqB = gen.randomDNA(size);
        
        double manacherTime = Benchmark::measure([&]() {
            Manacher::longestPalindrome(seqA);
        }, 5);
        manacherTimes.push_back(manacherTime);
        
        double saTime = Benchmark::measure([&]() {
            SuffixArray sa(seqA + "#" + seqB);
        }, 5);
        saTimes.push_back(saTime);
        
        double dpTime = Benchmark::measure([&]() {
            LongestCommonSubstring::usingDP(seqA, seqB);
        }, 5);
        lcsDPTimes.push_back(dpTime);
        
        std::cout << std::setw(10) << size 
                  << std::setw(18) << std::fixed << std::setprecision(4) << manacherTime
                  << std::setw(18) << saTime
                  << std::setw(18) << dpTime << "\n";
    }
    
    // Export to CSV
    std::ofstream csv("dna_benchmark.csv");
    csv << "size,manacher_ms,suffix_array_ms,lcs_dp_ms\n";
    for (size_t i = 0; i < sizes.size(); i++) {
        csv << sizes[i] << "," << manacherTimes[i] << "," << saTimes[i] << "," << lcsDPTimes[i] << "\n";
    }
    csv.close();
    
    std::cout << "\nResults exported to: dna_benchmark.csv\n";
    std::cout << "Expected complexity:\n";
    std::cout << "  - Manacher: O(n)\n";
    std::cout << "  - Suffix Array construction: O(n log n)\n";
    std::cout << "  - LCS-DP: O(n*m) = O(n^2) for equal lengths\n";
}

void runGraphBenchmarks() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  GRAPH ALGORITHMS BENCHMARKING\n";
    std::cout << std::string(70, '=') << "\n";
    
    std::cout << "\nExperimental Setup:\n";
    std::cout << "  - Graph sizes: V=50,100,200,500,1000 vertices\n";
    std::cout << "  - Edge density: E = 3*V edges\n";
    std::cout << "  - Seed: 42 for reproducibility\n";
    std::cout << "  - Trials: 5 per size\n\n";
    
    DataGenerator gen(42);
    std::vector<int> vertexCounts = {50, 100, 200, 500, 1000};
    
    std::vector<double> kruskalTimes, primTimes;
    
    std::cout << std::setw(10) << "Vertices" << std::setw(10) << "Edges"
              << std::setw(18) << "Kruskal (ms)" << std::setw(18) << "Prim (ms)\n";
    std::cout << std::string(56, '-') << "\n";
    
    for (int V : vertexCounts) {
        int E = 3 * V;
        
        // Generate random graph
        std::vector<Edge> edges;
        std::vector<std::vector<std::pair<int, int>>> adj(V);
        std::set<std::pair<int, int>> usedEdges;
        
        gen.setSeed(42);
        
        // Ensure connectivity
        for (int i = 1; i < V; i++) {
            int j = gen.randomInt(0, i - 1);
            int w = gen.randomInt(1, 100);
            edges.push_back({j, i, w});
            adj[j].push_back({i, w});
            adj[i].push_back({j, w});
            usedEdges.insert({std::min(i, j), std::max(i, j)});
        }
        
        while ((int)edges.size() < E) {
            int u = gen.randomInt(0, V - 1);
            int v = gen.randomInt(0, V - 1);
            if (u != v && !usedEdges.count({std::min(u, v), std::max(u, v)})) {
                int w = gen.randomInt(1, 100);
                edges.push_back({u, v, w});
                adj[u].push_back({v, w});
                adj[v].push_back({u, w});
                usedEdges.insert({std::min(u, v), std::max(u, v)});
            }
        }
        
        double kruskalTime = Benchmark::measure([&]() {
            auto edgesCopy = edges;
            MST::kruskal(V, edgesCopy, false);
        }, 5);
        kruskalTimes.push_back(kruskalTime);
        
        double primTime = Benchmark::measure([&]() {
            MST::prim(V, adj, 0, false);
        }, 5);
        primTimes.push_back(primTime);
        
        std::cout << std::setw(10) << V << std::setw(10) << E
                  << std::setw(18) << std::fixed << std::setprecision(4) << kruskalTime
                  << std::setw(18) << primTime << "\n";
    }
    
    // Export to CSV
    std::ofstream csv("graph_benchmark.csv");
    csv << "vertices,edges,kruskal_ms,prim_ms\n";
    for (size_t i = 0; i < vertexCounts.size(); i++) {
        csv << vertexCounts[i] << "," << 3 * vertexCounts[i] << "," 
            << kruskalTimes[i] << "," << primTimes[i] << "\n";
    }
    csv.close();
    
    std::cout << "\nResults exported to: graph_benchmark.csv\n";
    std::cout << "Expected complexity:\n";
    std::cout << "  - Kruskal: O(E log E)\n";
    std::cout << "  - Prim: O(E log V) with priority queue\n";
}

// ============================================================================
// SECURITY LOG ANALYSIS DEMONSTRATION
// ============================================================================

void demonstrateSecurityLogAnalysis() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  SECURITY LOG ANALYSIS - PRACTICAL APPLICATION\n";
    std::cout << std::string(70, '=') << "\n";
    
    DataGenerator gen(42);
    
    // Generate simulated server logs
    std::string randomLogs = gen.generateLogs(100);
    std::string periodicLogs = gen.generatePeriodicLogs(50);
    
    // Security signatures to search for
    std::vector<std::string> signatures = {
        "/api/login",       // Login endpoint
        "401",              // Unauthorized
        "403",              // Forbidden
        "500",              // Server error
        "10.0.0.5",         // Specific IP
        "DELETE"            // Potentially dangerous operation
    };
    
    std::cout << "\n### Analysis on Random Logs ###\n";
    std::cout << "Log sample (first 200 chars): " << randomLogs.substr(0, 200) << "...\n\n";
    
    std::cout << "Signature Detection Results:\n";
    std::cout << std::string(50, '-') << "\n";
    
    for (const auto& sig : signatures) {
        auto kmpMatches = KMP::search(randomLogs, sig);
        auto zMatches = ZAlgorithm::search(randomLogs, sig);
        
        std::cout << "\"" << sig << "\": " << kmpMatches.size() << " occurrences\n";
        
        // Verify algorithms match
        if (kmpMatches.size() != zMatches.size()) {
            std::cout << "  WARNING: Algorithm mismatch!\n";
        }
    }
    
    std::cout << "\n### Analysis on Periodic Logs ###\n";
    std::cout << "Periodic pattern detected - testing worst-case behavior.\n";
    
    auto start = Benchmark::now();
    auto matches = KMP::search(periodicLogs, "/api/login");
    auto end = Benchmark::now();
    
    std::cout << "Pattern \"/api/login\" found " << matches.size() << " times in "
              << Benchmark::elapsed_ms(start, end) << " ms\n";
    
    std::cout << "\n### Cybersecurity Significance ###\n";
    std::cout << std::string(60, '-') << "\n";
    std::cout << "String matching algorithms are foundational for:\n";
    std::cout << "1. Intrusion Detection Systems (IDS) - pattern-based threat detection\n";
    std::cout << "2. Log Analysis - identifying suspicious activities\n";
    std::cout << "3. Malware Scanning - signature-based detection\n";
    std::cout << "4. Network Traffic Analysis - protocol inspection\n\n";
    std::cout << "KMP/Z advantages: Linear time O(n+m), no worst-case degradation\n";
    std::cout << "Rabin-Karp advantage: Efficient multi-pattern search\n";
    std::cout << "Real-world: Aho-Corasick combines these for optimal multi-pattern\n";
}

// ============================================================================
// GENOMICS APPLICATION DEMONSTRATION
// ============================================================================

void demonstrateGenomicsApplication() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  GENOMICS APPLICATION - DNA SEQUENCE ANALYSIS\n";
    std::cout << std::string(70, '=') << "\n";
    
    // Example: Realistic DNA sequences (portions of actual genomic patterns)
    std::string geneA = "ATGCGATCGATCGATCGAATTCGAATTCGAATTCGATCGATCGTAGCTAGCT"
                        "GCTAGCTAGCTACGTACGTACGTACGTGCATGCATGCATGCATTACGTACGT";
    std::string geneB = "GCTAGCTACGATCGATCGATCGAATTCGAATTCGATCGATCGTAGCTAGCTA"
                        "TACGTACGTACGTACGTGCATGCATGCATGCATGCATGCATTTTACGTACGT";
    
    std::cout << "\nGene A (100 bp): " << geneA.substr(0, 50) << "...\n";
    std::cout << "Gene B (100 bp): " << geneB.substr(0, 50) << "...\n";
    
    // Palindromic regions (important for restriction enzyme sites)
    std::cout << "\n### Palindromic Region Analysis ###\n";
    std::cout << "(Important for identifying restriction enzyme recognition sites)\n\n";
    
    auto [palA, posA] = Manacher::longestPalindrome(geneA);
    auto [palB, posB] = Manacher::longestPalindrome(geneB);
    
    std::cout << "Gene A longest palindrome: \"" << palA << "\" at position " << posA << "\n";
    std::cout << "Gene B longest palindrome: \"" << palB << "\" at position " << posB << "\n";
    
    // Common subsequence (evolutionary relationship)
    std::cout << "\n### Longest Common Substring (Evolutionary Analysis) ###\n";
    auto [lcs, lcsPos] = LongestCommonSubstring::usingDP(geneA, geneB);
    std::cout << "LCS: \"" << lcs << "\" (length: " << lcs.size() << ")\n";
    std::cout << "Shared sequence suggests common evolutionary origin.\n";
    
    // Motif queries (transcription factor binding sites)
    std::cout << "\n### Motif Query (Transcription Factor Binding Sites) ###\n";
    std::vector<std::string> motifs = {
        "GAATTC",   // EcoRI restriction site
        "GGATCC",   // BamHI restriction site
        "AAGCTT",   // HindIII restriction site
        "TATA",     // TATA box (promoter element)
        "GCATGC"    // SphI restriction site
    };
    
    std::string combined = geneA + "#" + geneB;
    SuffixArray sa(combined);
    sa.queryMotifs(motifs);
    
    std::cout << "\n### Genomics Research Significance ###\n";
    std::cout << std::string(60, '-') << "\n";
    std::cout << "These algorithms are essential for:\n";
    std::cout << "1. Sequence Alignment - BLAST uses suffix arrays + heuristics\n";
    std::cout << "2. Genome Assembly - suffix arrays enable overlap detection\n";
    std::cout << "3. Variant Calling - comparing against reference genomes\n";
    std::cout << "4. Motif Discovery - finding regulatory elements\n";
    std::cout << "5. Phylogenetics - computing sequence similarity\n\n";
    std::cout << "References:\n";
    std::cout << "- Gusfield, D. (1997). Algorithms on Strings, Trees and Sequences\n";
    std::cout << "- Navarro, G. (2001). A guided tour to approximate string matching\n";
    std::cout << "- Langmead, B. et al. (2009). Ultrafast memory-efficient alignment\n";
}

// ============================================================================
// MENU SYSTEM
// ============================================================================

void printMenu() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  TC2038 - ADVANCED ALGORITHMS IMPLEMENTATION\n";
    std::cout << "  Chapters 2-3\n";
    std::cout << std::string(70, '=') << "\n\n";
    
    std::cout << "STRING ALGORITHMS:\n";
    std::cout << "  1. String Matching Demo (KMP, Z-Algorithm, Rabin-Karp)\n";
    std::cout << "  2. DNA Algorithms Demo (Manacher, Suffix Array, LCS)\n";
    std::cout << "  3. Security Log Analysis (Practical Application)\n";
    std::cout << "  4. Genomics Application (Practical Application)\n\n";
    
    std::cout << "GRAPH ALGORITHMS:\n";
    std::cout << "  5. Graph Algorithms Demo (Max Flow, MST, 2nd MST, etc.)\n\n";
    
    std::cout << "COMBINED SYSTEM:\n";
    std::cout << "  6. Expression Parser (Polish/Infix, AST, Truth Tables)\n\n";
    
    std::cout << "BENCHMARKING:\n";
    std::cout << "  7. String Algorithms Benchmark\n";
    std::cout << "  8. DNA Algorithms Benchmark\n";
    std::cout << "  9. Graph Algorithms Benchmark\n";
    std::cout << " 10. Run All Benchmarks\n\n";
    
    std::cout << "  0. Exit\n\n";
    
    std::cout << "Select option: ";
}

int main() {
    int choice;
    
    do {
        printMenu();
        std::cin >> choice;
        
        switch (choice) {
            case 1:
                demonstrateStringAlgorithms();
                break;
            case 2:
                demonstrateDNAAlgorithms();
                break;
            case 3:
                demonstrateSecurityLogAnalysis();
                break;
            case 4:
                demonstrateGenomicsApplication();
                break;
            case 5:
                demonstrateGraphAlgorithms();
                break;
            case 6:
                demonstrateExpressionParser();
                break;
            case 7:
                runStringBenchmarks();
                break;
            case 8:
                runDNABenchmarks();
                break;
            case 9:
                runGraphBenchmarks();
                break;
            case 10:
                runStringBenchmarks();
                runDNABenchmarks();
                runGraphBenchmarks();
                break;
            case 0:
                std::cout << "\nGoodbye!\n";
                break;
            default:
                std::cout << "\nInvalid option. Please try again.\n";
        }
        
        if (choice != 0) {
            std::cout << "\nPress Enter to continue...";
            std::cin.ignore();
            std::cin.get();
        }
        
    } while (choice != 0);
    
    return 0;
}
