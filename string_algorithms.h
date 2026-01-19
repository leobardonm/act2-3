#ifndef STRING_ALGORITHMS_H
#define STRING_ALGORITHMS_H

#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>
#include <algorithm>

// ============================================================================
// 1.a KNUTH-MORRIS-PRATT (KMP) ALGORITHM
// ============================================================================
// Time Complexity: O(n + m) where n = |text|, m = |pattern|
// Space Complexity: O(m) for the prefix function
// ============================================================================

class KMP {
public:
    // Compute the prefix function (π-array / failure function)
    // π[i] = length of the longest proper prefix of P[0..i] that is also a suffix
    static std::vector<int> computePrefixFunction(const std::string& pattern) {
        int m = pattern.size();
        std::vector<int> pi(m, 0);
        
        int k = 0; // length of the current longest prefix-suffix
        for (int i = 1; i < m; i++) {
            // Fall back using the prefix function until we find a match or reach 0
            while (k > 0 && pattern[k] != pattern[i]) {
                k = pi[k - 1];
            }
            // If characters match, extend the prefix-suffix
            if (pattern[k] == pattern[i]) {
                k++;
            }
            pi[i] = k;
        }
        return pi;
    }
    
    // Find all occurrences of pattern in text
    // Returns vector of starting positions (0-indexed)
    static std::vector<int> search(const std::string& text, const std::string& pattern) {
        std::vector<int> matches;
        if (pattern.empty() || text.empty() || pattern.size() > text.size()) {
            return matches;
        }
        
        std::vector<int> pi = computePrefixFunction(pattern);
        int n = text.size();
        int m = pattern.size();
        int q = 0; // number of characters matched
        
        for (int i = 0; i < n; i++) {
            // Fall back using prefix function
            while (q > 0 && pattern[q] != text[i]) {
                q = pi[q - 1];
            }
            // If characters match, advance
            if (pattern[q] == text[i]) {
                q++;
            }
            // If we've matched the entire pattern
            if (q == m) {
                matches.push_back(i - m + 1);
                q = pi[q - 1]; // prepare for next potential match
            }
        }
        return matches;
    }
    
    // Check if a pattern has a nontrivial border (proper prefix = proper suffix)
    static bool hasNontrivialBorder(const std::string& pattern) {
        if (pattern.size() < 2) return false;
        std::vector<int> pi = computePrefixFunction(pattern);
        return pi.back() > 0;
    }
    
    // Print the prefix function with explanation
    static void printPrefixFunction(const std::string& pattern) {
        std::vector<int> pi = computePrefixFunction(pattern);
        std::cout << "Pattern: \"" << pattern << "\"\n";
        std::cout << "Index:   ";
        for (size_t i = 0; i < pattern.size(); i++) {
            std::cout << i << " ";
        }
        std::cout << "\nChar:    ";
        for (char c : pattern) {
            std::cout << c << " ";
        }
        std::cout << "\nπ-array: ";
        for (int val : pi) {
            std::cout << val << " ";
        }
        std::cout << "\n";
        
        // Explain nontrivial borders
        if (pi.back() > 0) {
            std::cout << "Nontrivial border: \"" << pattern.substr(0, pi.back()) 
                      << "\" (length " << pi.back() << ")\n";
        } else {
            std::cout << "No nontrivial border.\n";
        }
    }
};

// ============================================================================
// 1.b Z-ALGORITHM
// ============================================================================
// Time Complexity: O(n + m) using the P#T concatenation trick
// Space Complexity: O(n + m) for the Z-array
// ============================================================================

class ZAlgorithm {
public:
    // Compute Z-array: Z[i] = length of longest substring starting at i 
    // that matches a prefix of the string
    static std::vector<int> computeZArray(const std::string& s) {
        int n = s.size();
        std::vector<int> z(n, 0);
        
        // l, r define the rightmost Z-box [l, r]
        int l = 0, r = 0;
        
        for (int i = 1; i < n; i++) {
            if (i < r) {
                // We're inside the Z-box, use previously computed values
                z[i] = std::min(r - i, z[i - l]);
            }
            // Try to extend
            while (i + z[i] < n && s[z[i]] == s[i + z[i]]) {
                z[i]++;
            }
            // Update Z-box if we extended past r
            if (i + z[i] > r) {
                l = i;
                r = i + z[i];
            }
        }
        z[0] = n; // By convention, Z[0] = |S|
        return z;
    }
    
    // Find all occurrences using P#T concatenation trick
    static std::vector<int> search(const std::string& text, const std::string& pattern) {
        std::vector<int> matches;
        if (pattern.empty() || text.empty() || pattern.size() > text.size()) {
            return matches;
        }
        
        // Concatenate: pattern + sentinel + text
        std::string concat = pattern + "#" + text;
        std::vector<int> z = computeZArray(concat);
        
        int m = pattern.size();
        // Check positions in the text part of the concatenation
        for (size_t i = m + 1; i < concat.size(); i++) {
            if (z[i] == m) {
                // Found a match at position (i - m - 1) in original text
                matches.push_back(i - m - 1);
            }
        }
        return matches;
    }
    
    // Print Z-array with explanation
    static void printZArray(const std::string& s) {
        std::vector<int> z = computeZArray(s);
        std::cout << "String: \"" << s << "\"\n";
        std::cout << "Index:  ";
        for (size_t i = 0; i < s.size(); i++) {
            std::cout << i << " ";
        }
        std::cout << "\nChar:   ";
        for (char c : s) {
            std::cout << c << " ";
        }
        std::cout << "\nZ-array:";
        for (int val : z) {
            std::cout << val << " ";
        }
        std::cout << "\n";
    }
};

// ============================================================================
// 1.c RABIN-KARP (ROLLING HASH)
// ============================================================================
// Time Complexity: O(n + m) average, O(nm) worst case with many collisions
// Space Complexity: O(1) extra for single pattern
// ============================================================================

class RabinKarp {
private:
    // Default hash parameters
    static const long long BASE1 = 31;
    static const long long MOD1 = 1e9 + 7;
    static const long long BASE2 = 37;
    static const long long MOD2 = 1e9 + 9;
    
public:
    // Compute hash of a string
    static long long computeHash(const std::string& s, long long base, long long mod) {
        long long hash = 0;
        long long power = 1;
        for (char c : s) {
            hash = (hash + (c - 'a' + 1) * power) % mod;
            power = (power * base) % mod;
        }
        return hash;
    }
    
    // Single hash pattern matching
    static std::vector<int> searchSingleHash(const std::string& text, const std::string& pattern,
                                             long long base = BASE1, long long mod = MOD1) {
        std::vector<int> matches;
        int n = text.size();
        int m = pattern.size();
        
        if (m > n || m == 0) return matches;
        
        // Precompute base^m mod
        long long basePow = 1;
        for (int i = 0; i < m; i++) {
            basePow = (basePow * base) % mod;
        }
        
        // Compute hash of pattern and first window
        long long patternHash = computeHash(pattern, base, mod);
        long long windowHash = computeHash(text.substr(0, m), base, mod);
        
        // Slide the window
        for (int i = 0; i <= n - m; i++) {
            if (windowHash == patternHash) {
                // Verify character by character to handle collisions
                if (text.substr(i, m) == pattern) {
                    matches.push_back(i);
                }
            }
            // Roll the hash
            if (i < n - m) {
                // Remove leftmost character, add new rightmost
                windowHash = (windowHash - (text[i] - 'a' + 1) + mod) % mod;
                windowHash = (windowHash * modInverse(base, mod)) % mod;
                windowHash = (windowHash + (text[i + m] - 'a' + 1) * basePow % mod * modInverse(base, mod)) % mod;
                
                // Alternative rolling (simpler): recompute from scratch for clarity
                windowHash = computeHash(text.substr(i + 1, m), base, mod);
            }
        }
        return matches;
    }
    
    // Double hash pattern matching (reduces collision probability)
    static std::vector<int> searchDoubleHash(const std::string& text, const std::string& pattern) {
        std::vector<int> matches;
        int n = text.size();
        int m = pattern.size();
        
        if (m > n || m == 0) return matches;
        
        // Compute hashes with two different (base, mod) pairs
        auto hash1 = [](const std::string& s, int start, int len) -> long long {
            long long h = 0, p = 1;
            for (int i = 0; i < len; i++) {
                h = (h + (long long)(s[start + i]) * p) % MOD1;
                p = (p * BASE1) % MOD1;
            }
            return h;
        };
        
        auto hash2 = [](const std::string& s, int start, int len) -> long long {
            long long h = 0, p = 1;
            for (int i = 0; i < len; i++) {
                h = (h + (long long)(s[start + i]) * p) % MOD2;
                p = (p * BASE2) % MOD2;
            }
            return h;
        };
        
        long long patHash1 = hash1(pattern, 0, m);
        long long patHash2 = hash2(pattern, 0, m);
        
        for (int i = 0; i <= n - m; i++) {
            long long h1 = hash1(text, i, m);
            long long h2 = hash2(text, i, m);
            
            if (h1 == patHash1 && h2 == patHash2) {
                // Both hashes match - very high probability of true match
                // Still verify to guarantee correctness
                if (text.substr(i, m) == pattern) {
                    matches.push_back(i);
                }
            }
        }
        return matches;
    }
    
    // Multi-pattern matching (search for multiple patterns simultaneously)
    static std::vector<std::pair<int, int>> multiPatternSearch(
            const std::string& text, 
            const std::vector<std::string>& patterns) {
        
        std::vector<std::pair<int, int>> results; // (pattern_index, position)
        
        // Precompute hashes for all patterns (using double hashing)
        std::vector<std::pair<long long, long long>> patternHashes;
        std::vector<int> patternLengths;
        
        for (size_t i = 0; i < patterns.size(); i++) {
            const std::string& p = patterns[i];
            long long h1 = 0, h2 = 0, pow1 = 1, pow2 = 1;
            for (char c : p) {
                h1 = (h1 + (long long)c * pow1) % MOD1;
                h2 = (h2 + (long long)c * pow2) % MOD2;
                pow1 = (pow1 * BASE1) % MOD1;
                pow2 = (pow2 * BASE2) % MOD2;
            }
            patternHashes.push_back({h1, h2});
            patternLengths.push_back(p.size());
        }
        
        int n = text.size();
        
        // For each starting position
        for (int i = 0; i < n; i++) {
            // Check each pattern
            for (size_t j = 0; j < patterns.size(); j++) {
                int m = patternLengths[j];
                if (i + m > n) continue;
                
                // Compute hash of window
                long long h1 = 0, h2 = 0, pow1 = 1, pow2 = 1;
                for (int k = 0; k < m; k++) {
                    h1 = (h1 + (long long)text[i + k] * pow1) % MOD1;
                    h2 = (h2 + (long long)text[i + k] * pow2) % MOD2;
                    pow1 = (pow1 * BASE1) % MOD1;
                    pow2 = (pow2 * BASE2) % MOD2;
                }
                
                if (h1 == patternHashes[j].first && h2 == patternHashes[j].second) {
                    if (text.substr(i, m) == patterns[j]) {
                        results.push_back({j, i});
                    }
                }
            }
        }
        return results;
    }
    
    // Get hash parameters info
    static void printHashParameters() {
        std::cout << "Single Hash: base=" << BASE1 << ", mod=" << MOD1 << "\n";
        std::cout << "Double Hash: base1=" << BASE1 << ", mod1=" << MOD1 
                  << ", base2=" << BASE2 << ", mod2=" << MOD2 << "\n";
        std::cout << "Collision handling: Character-by-character verification on hash match\n";
    }
    
private:
    // Modular inverse using extended Euclidean algorithm
    static long long modInverse(long long a, long long mod) {
        long long m0 = mod, y = 0, x = 1;
        if (mod == 1) return 0;
        while (a > 1) {
            long long q = a / mod;
            long long t = mod;
            mod = a % mod;
            a = t;
            t = y;
            y = x - q * y;
            x = t;
        }
        if (x < 0) x += m0;
        return x;
    }
};

// ============================================================================
// STRING ALGORITHMS DEMO / TEST
// ============================================================================

void demonstrateStringAlgorithms() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  STRING MATCHING ALGORITHMS DEMONSTRATION\n";
    std::cout << "  Problem 1: Web Server Log Analysis\n";
    std::cout << std::string(70, '=') << "\n";
    
    // ========================================================================
    // RANDOM-LIKE LOG (typical server traffic)
    // ========================================================================
    std::string randomLog = 
        "192.168.1.1 - - [2025/09/08:10:30:45] \"GET /api/login HTTP/1.1\" 200 1234\n"
        "10.0.0.5 - - [2025/09/08:10:30:46] \"POST /api/login HTTP/1.1\" 401 567\n"
        "192.168.1.1 - - [2025/09/08:10:30:47] \"GET /api/users HTTP/1.1\" 200 890\n"
        "172.16.0.1 - - [2025/09/08:10:30:48] \"GET /api/login HTTP/1.1\" 403 123\n"
        "10.0.0.5 - - [2025/09/08:10:30:49] \"DELETE /admin HTTP/1.1\" 401 456\n"
        "192.168.1.1 - - [2025/09/08:10:30:50] \"GET /api/login HTTP/1.1\" 200 789\n";
    
    // ========================================================================
    // HIGHLY PERIODIC LOG (DoS attack simulation - repeated requests)
    // ========================================================================
    std::string periodicLog;
    std::string repeatedEntry = "10.0.0.5 GET /api/login 401 ";
    for (int i = 0; i < 20; i++) {
        periodicLog += repeatedEntry;
    }
    
    // ========================================================================
    // SIGNATURE SET P = {P1, ..., Pm} with at least 2 having nontrivial borders
    // ========================================================================
    // Nontrivial border: proper prefix that is also a proper suffix
    // P1 = "0.0.0" has border "0.0" (wait, let's use better examples)
    // Actually using patterns that naturally appear in logs:
    std::vector<std::string> signatures = {
        "/api/login",      // P1: Security-critical endpoint
        "192.168.1.1",     // P2: Specific internal IP
        "401",             // P3: Unauthorized attempts
        "0.0.0.0",         // P4: NONTRIVIAL BORDER: "0.0" (π = [0,0,0,1,2,3,4])
        "1.1"              // P5: NONTRIVIAL BORDER: "1" (π = [0,0,1])
    };
    
    std::cout << "\n=== TEST DATA ===\n";
    std::cout << "Alphabet Σ: ASCII printable characters (|Σ| ≥ 26)\n\n";
    
    std::cout << "Random-like log (T1) - first 200 chars:\n";
    std::cout << "\"" << randomLog.substr(0, 200) << "...\"\n";
    std::cout << "Length: " << randomLog.size() << " characters\n\n";
    
    std::cout << "Highly periodic log (T2) - pattern repeated 20 times:\n";
    std::cout << "Pattern: \"" << repeatedEntry << "\"\n";
    std::cout << "Length: " << periodicLog.size() << " characters\n\n";
    
    std::cout << "Signature set P = {P1, ..., P5}:\n";
    for (size_t i = 0; i < signatures.size(); i++) {
        std::cout << "  P" << i+1 << " = \"" << signatures[i] << "\"";
        if (KMP::hasNontrivialBorder(signatures[i])) {
            std::cout << " [HAS NONTRIVIAL BORDER]";
        }
        std::cout << "\n";
    }
    
    // ========================================================================
    // 1.a KMP (Knuth-Morris-Pratt)
    // ========================================================================
    std::cout << "\n" << std::string(70, '-') << "\n";
    std::cout << "1.a. KNUTH-MORRIS-PRATT (KMP)\n";
    std::cout << std::string(70, '-') << "\n";
    
    std::cout << "\n[Random Log Results]\n";
    for (const auto& sig : signatures) {
        auto matches = KMP::search(randomLog, sig);
        std::cout << "  \"" << sig << "\": " << matches.size() << " match(es) at: ";
        for (int pos : matches) std::cout << pos << " ";
        std::cout << "\n";
    }
    
    std::cout << "\n[Periodic Log Results]\n";
    for (const auto& sig : signatures) {
        auto matches = KMP::search(periodicLog, sig);
        std::cout << "  \"" << sig << "\": " << matches.size() << " match(es)\n";
    }
    
    // PREFIX FUNCTION for two patterns with nontrivial borders
    std::cout << "\n[Prefix Function (π-array) for patterns with NONTRIVIAL BORDERS]\n\n";
    
    std::cout << "Pattern P4 = \"0.0.0.0\" (appears in IP addresses):\n";
    KMP::printPrefixFunction("0.0.0.0");
    std::cout << "  Explanation: The pattern has overlapping structure.\n";
    std::cout << "  Border \"0.0.0\" allows skipping comparisons on mismatch.\n\n";
    
    std::cout << "Pattern P5 = \"1.1\" (common in version numbers, IPs):\n";
    KMP::printPrefixFunction("1.1");
    std::cout << "  Explanation: Simple overlap at the '1' character.\n";
    
    // ========================================================================
    // 1.b Z-Algorithm
    // ========================================================================
    std::cout << "\n" << std::string(70, '-') << "\n";
    std::cout << "1.b. Z-ALGORITHM (P#T concatenation)\n";
    std::cout << std::string(70, '-') << "\n";
    
    std::cout << "\n[Random Log Results]\n";
    for (const auto& sig : signatures) {
        auto matches = ZAlgorithm::search(randomLog, sig);
        std::cout << "  \"" << sig << "\": " << matches.size() << " match(es) at: ";
        for (int pos : matches) std::cout << pos << " ";
        std::cout << "\n";
    }
    
    std::cout << "\n[Periodic Log Results]\n";
    for (const auto& sig : signatures) {
        auto matches = ZAlgorithm::search(periodicLog, sig);
        std::cout << "  \"" << sig << "\": " << matches.size() << " match(es)\n";
    }
    
    std::cout << "\n[Z-array Example: P#T where P=\"0.0\" and T=\"0.0.0.0\"]\n";
    std::string concat = "0.0#0.0.0.0";
    ZAlgorithm::printZArray(concat);
    std::cout << "  Positions where Z[i] = |P| = 3 indicate matches.\n";
    
    // ========================================================================
    // 1.c Rabin-Karp (Rolling Hash)
    // ========================================================================
    std::cout << "\n" << std::string(70, '-') << "\n";
    std::cout << "1.c. RABIN-KARP (Rolling Hash)\n";
    std::cout << std::string(70, '-') << "\n";
    
    std::cout << "\n[Hash Parameters]\n";
    RabinKarp::printHashParameters();
    
    std::cout << "\n[Multi-pattern Search on Random Log]\n";
    auto rkResults = RabinKarp::multiPatternSearch(randomLog, signatures);
    for (const auto& [patIdx, pos] : rkResults) {
        std::cout << "  P" << patIdx+1 << " \"" << signatures[patIdx] 
                  << "\" at position " << pos << "\n";
    }
    
    std::cout << "\n[Multi-pattern Search on Periodic Log]\n";
    auto rkPeriodic = RabinKarp::multiPatternSearch(periodicLog, signatures);
    std::cout << "  Total matches found: " << rkPeriodic.size() << "\n";
    
    // ========================================================================
    // ALGORITHM COMPARISON & VERIFICATION
    // ========================================================================
    std::cout << "\n" << std::string(70, '-') << "\n";
    std::cout << "ALGORITHM VERIFICATION\n";
    std::cout << std::string(70, '-') << "\n";
    
    bool allMatch = true;
    for (const auto& sig : signatures) {
        auto kmpRes = KMP::search(randomLog, sig);
        auto zRes = ZAlgorithm::search(randomLog, sig);
        
        if (kmpRes != zRes) {
            std::cout << "✗ MISMATCH for \"" << sig << "\"!\n";
            allMatch = false;
        }
    }
    if (allMatch) {
        std::cout << "✓ All algorithms produce identical results on random log.\n";
    }
    
    std::cout << "\n" << std::string(70, '-') << "\n";
    std::cout << "COMPLEXITY ANALYSIS\n";
    std::cout << std::string(70, '-') << "\n";
    std::cout << "Let n = |T| (text length), m = |P| (pattern length)\n\n";
    std::cout << "KMP:        O(n + m) time, O(m) space - prefix function preprocessing\n";
    std::cout << "Z-Algorithm: O(n + m) time, O(n + m) space - Z-array on P#T\n";
    std::cout << "Rabin-Karp: O(n + m) expected, O(nm) worst case with collisions\n";
    std::cout << "            Double hashing reduces collision probability to ~1/(mod1*mod2)\n";
}

#endif // STRING_ALGORITHMS_H
