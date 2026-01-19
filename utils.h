#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <random>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <functional>


class Benchmark {
public:
    using TimePoint = std::chrono::high_resolution_clock::time_point;
    
    static TimePoint now() {
        return std::chrono::high_resolution_clock::now();
    }
    
    static double elapsed_ms(TimePoint start, TimePoint end) {
        return std::chrono::duration<double, std::milli>(end - start).count();
    }
    
    static double elapsed_us(TimePoint start, TimePoint end) {
        return std::chrono::duration<double, std::micro>(end - start).count();
    }
    
    // Run a function multiple times and return average time in milliseconds
    template<typename Func>
    static double measure(Func f, int trials = 10) {
        double total = 0;
        for (int i = 0; i < trials; i++) {
            auto start = now();
            f();
            auto end = now();
            total += elapsed_ms(start, end);
        }
        return total / trials;
    }
    
    // Export benchmark results to CSV for plotting
    static void exportCSV(const std::string& filename,
                         const std::vector<int>& sizes,
                         const std::vector<double>& times,
                         const std::string& algorithm) {
        std::ofstream file(filename, std::ios::app);
        file << "# Algorithm: " << algorithm << "\n";
        file << "size,time_ms\n";
        for (size_t i = 0; i < sizes.size(); i++) {
            file << sizes[i] << "," << times[i] << "\n";
        }
        file << "\n";
        file.close();
    }
};

// ============================================================================
// DATA GENERATORS
// ============================================================================

class DataGenerator {
private:
    std::mt19937 rng;
    
public:
    DataGenerator(unsigned int seed = 42) : rng(seed) {}
    
    void setSeed(unsigned int seed) { rng.seed(seed); }
    
    // Generate random string with given alphabet
    std::string randomString(int length, const std::string& alphabet = "abcdefghijklmnopqrstuvwxyz") {
        std::uniform_int_distribution<int> dist(0, alphabet.size() - 1);
        std::string result(length, ' ');
        for (int i = 0; i < length; i++) {
            result[i] = alphabet[dist(rng)];
        }
        return result;
    }
    
    // Generate DNA sequence
    std::string randomDNA(int length) {
        return randomString(length, "ACGT");
    }
    
    // Generate periodic string (for testing edge cases)
    std::string periodicString(int length, const std::string& pattern) {
        std::string result;
        while ((int)result.size() < length) {
            result += pattern;
        }
        return result.substr(0, length);
    }
    
    // Generate simulated server log entry
    std::string generateLogEntry() {
        std::vector<std::string> ips = {"192.168.1.1", "10.0.0.5", "172.16.0.100", "8.8.8.8"};
        std::vector<std::string> methods = {"GET", "POST", "PUT", "DELETE"};
        std::vector<std::string> paths = {"/api/login", "/api/users", "/admin/dashboard", "/api/data"};
        std::vector<std::string> statuses = {"200", "401", "403", "404", "500"};
        
        std::uniform_int_distribution<int> ipDist(0, ips.size() - 1);
        std::uniform_int_distribution<int> methodDist(0, methods.size() - 1);
        std::uniform_int_distribution<int> pathDist(0, paths.size() - 1);
        std::uniform_int_distribution<int> statusDist(0, statuses.size() - 1);
        
        return ips[ipDist(rng)] + " - - [2025/09/08:10:30:45] \"" + 
               methods[methodDist(rng)] + " " + paths[pathDist(rng)] + " HTTP/1.1\" " +
               statuses[statusDist(rng)] + " " + std::to_string(1000 + rng() % 9000);
    }
    
    // Generate multiple log entries
    std::string generateLogs(int count) {
        std::string logs;
        for (int i = 0; i < count; i++) {
            logs += generateLogEntry() + "\n";
        }
        return logs;
    }
    
    // Generate highly periodic logs (for testing)
    std::string generatePeriodicLogs(int count) {
        std::string baseLog = "192.168.1.1 - - [2025/09/08:10:30:45] \"GET /api/login HTTP/1.1\" 200 1234\n";
        return periodicString(count * baseLog.size(), baseLog);
    }
    
    // Generate random integer
    int randomInt(int min, int max) {
        std::uniform_int_distribution<int> dist(min, max);
        return dist(rng);
    }
    
    // Generate random weighted graph as adjacency list with weights
    std::vector<std::vector<std::pair<int, int>>> randomWeightedGraph(int V, int E, int maxWeight = 100) {
        std::vector<std::vector<std::pair<int, int>>> adj(V);
        std::vector<std::pair<int, int>> edges;
        
        // First ensure connectivity with a spanning tree
        for (int i = 1; i < V; i++) {
            int j = randomInt(0, i - 1);
            int w = randomInt(1, maxWeight);
            adj[i].push_back({j, w});
            adj[j].push_back({i, w});
            edges.push_back({i, j});
        }
        
        // Add remaining edges
        int remaining = E - (V - 1);
        while (remaining > 0) {
            int u = randomInt(0, V - 1);
            int v = randomInt(0, V - 1);
            if (u != v) {
                bool exists = false;
                for (auto& p : adj[u]) {
                    if (p.first == v) { exists = true; break; }
                }
                if (!exists) {
                    int w = randomInt(1, maxWeight);
                    adj[u].push_back({v, w});
                    adj[v].push_back({u, w});
                    remaining--;
                }
            }
        }
        return adj;
    }
};

// ============================================================================
// PRINT UTILITIES
// ============================================================================

template<typename T>
void printVector(const std::vector<T>& v, const std::string& name = "") {
    if (!name.empty()) std::cout << name << ": ";
    std::cout << "[";
    for (size_t i = 0; i < v.size(); i++) {
        std::cout << v[i];
        if (i < v.size() - 1) std::cout << ", ";
    }
    std::cout << "]\n";
}

template<typename T>
void printMatrix(const std::vector<std::vector<T>>& m, const std::string& name = "") {
    if (!name.empty()) std::cout << name << ":\n";
    for (const auto& row : m) {
        for (const auto& val : row) {
            std::cout << std::setw(4) << val << " ";
        }
        std::cout << "\n";
    }
}

void printSeparator(const std::string& title = "") {
    std::cout << "\n" << std::string(60, '=') << "\n";
    if (!title.empty()) {
        std::cout << "  " << title << "\n";
        std::cout << std::string(60, '=') << "\n";
    }
}

#endif // UTILS_H
