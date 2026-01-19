#ifndef DNA_ALGORITHMS_H
#define DNA_ALGORITHMS_H

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <unordered_map>

// ============================================================================
// 2.a MANACHER'S ALGORITHM - Longest Palindromic Substring
// ============================================================================
// Time Complexity: O(n)
// Space Complexity: O(n)
// 
// CORRECTNESS JUSTIFICATION (Radius Propagation):
// The algorithm maintains a "rightmost palindrome" [l, r] seen so far.
// For each position i:
//   1. If i < r, we mirror position i' = l + (r - i) inside [l, r]
//   2. The radius at i is at least min(p[i'], r - i) due to symmetry
//   3. We then try to extend beyond this known bound
// This avoids redundant comparisons since characters within [l, r] 
// have already been matched during expansion of the center (l+r)/2.
// Each character is visited at most twice: once for expansion, once 
// for mirror copying. Hence O(n) total.
// ============================================================================

class Manacher {
public:
    // Transform string: "abc" -> "^#a#b#c#$"
    // This handles both odd and even length palindromes uniformly
    static std::string transform(const std::string& s) {
        std::string t = "^#";
        for (char c : s) {
            t += c;
            t += '#';
        }
        t += '$';
        return t;
    }
    
    // Compute palindrome radii array
    // p[i] = radius of palindrome centered at i in transformed string
    static std::vector<int> computeRadii(const std::string& t) {
        int n = t.size();
        std::vector<int> p(n, 0);
        int center = 0, right = 0;
        
        for (int i = 1; i < n - 1; i++) {
            // Mirror index
            int mirror = 2 * center - i;
            
            // If within right boundary, use mirror's radius as starting point
            if (i < right) {
                p[i] = std::min(right - i, p[mirror]);
            }
            
            // Try to expand palindrome centered at i
            while (t[i + 1 + p[i]] == t[i - 1 - p[i]]) {
                p[i]++;
            }
            
            // If expanded beyond right, update center and right
            if (i + p[i] > right) {
                center = i;
                right = i + p[i];
            }
        }
        return p;
    }
    
    // Find the longest palindromic substring
    static std::pair<std::string, int> longestPalindrome(const std::string& s) {
        if (s.empty()) return {"", -1};
        
        std::string t = transform(s);
        std::vector<int> p = computeRadii(t);
        
        // Find max radius and its center
        int maxLen = 0, centerIdx = 0;
        for (int i = 1; i < (int)t.size() - 1; i++) {
            if (p[i] > maxLen) {
                maxLen = p[i];
                centerIdx = i;
            }
        }
        
        // Convert back to original string position
        // In transformed string, position i corresponds to (i-1)/2 in original
        int start = (centerIdx - maxLen) / 2;
        std::string palindrome = s.substr(start, maxLen);
        
        return {palindrome, start};
    }
    
    // Find all maximal palindromic substrings (for detailed analysis)
    static std::vector<std::tuple<std::string, int, int>> allMaximalPalindromes(
            const std::string& s, int minLength = 2) {
        std::vector<std::tuple<std::string, int, int>> result;
        std::string t = transform(s);
        std::vector<int> p = computeRadii(t);
        
        for (int i = 1; i < (int)t.size() - 1; i++) {
            if (p[i] >= minLength) {
                int start = (i - p[i]) / 2;
                std::string pal = s.substr(start, p[i]);
                result.push_back({pal, start, p[i]});
            }
        }
        return result;
    }
    
    static void printRadiiExplanation(const std::string& s) {
        std::cout << "Original: \"" << s << "\"\n";
        std::string t = transform(s);
        std::cout << "Transformed: \"" << t << "\"\n";
        std::vector<int> p = computeRadii(t);
        
        std::cout << "Index:  ";
        for (size_t i = 0; i < t.size(); i++) std::cout << i % 10 << " ";
        std::cout << "\nChar:   ";
        for (char c : t) std::cout << c << " ";
        std::cout << "\nRadius: ";
        for (int r : p) std::cout << r << " ";
        std::cout << "\n";
    }
};

// ============================================================================
// 2.b SUFFIX ARRAY WITH KASAI LCP
// ============================================================================
// Suffix Array Construction: O(n log n) using prefix doubling
// Kasai LCP: O(n)
// Space: O(n)
// ============================================================================

class SuffixArray {
public:
    std::string text;
    std::vector<int> sa;    // suffix array
    std::vector<int> lcp;   // LCP array (Kasai's algorithm)
    std::vector<int> rank_; // rank array (inverse of SA)
    
    SuffixArray() {}
    
    SuffixArray(const std::string& s) : text(s) {
        buildSuffixArray();
        buildLCP();
    }
    
    // Build suffix array using prefix doubling (O(n log n))
    void buildSuffixArray() {
        int n = text.size();
        if (n == 0) return;
        
        sa.resize(n);
        rank_.resize(n);
        std::vector<int> tmp(n);
        
        // Initialize with single character ranks
        for (int i = 0; i < n; i++) {
            sa[i] = i;
            rank_[i] = text[i];
        }
        
        // Prefix doubling
        for (int k = 1; k < n; k *= 2) {
            // Sort by (rank[i], rank[i+k])
            auto cmp = [&](int a, int b) {
                if (rank_[a] != rank_[b]) return rank_[a] < rank_[b];
                int ra = (a + k < n) ? rank_[a + k] : -1;
                int rb = (b + k < n) ? rank_[b + k] : -1;
                return ra < rb;
            };
            std::sort(sa.begin(), sa.end(), cmp);
            
            // Compute new ranks
            tmp[sa[0]] = 0;
            for (int i = 1; i < n; i++) {
                tmp[sa[i]] = tmp[sa[i-1]] + (cmp(sa[i-1], sa[i]) ? 1 : 0);
            }
            rank_ = tmp;
            
            // Early termination if all ranks are unique
            if (rank_[sa[n-1]] == n - 1) break;
        }
    }
    
    // Build LCP array using Kasai's algorithm (O(n))
    void buildLCP() {
        int n = text.size();
        if (n == 0) return;
        
        lcp.resize(n, 0);
        
        int k = 0;
        for (int i = 0; i < n; i++) {
            if (rank_[i] == 0) {
                k = 0;
                continue;
            }
            int j = sa[rank_[i] - 1];
            while (i + k < n && j + k < n && text[i + k] == text[j + k]) {
                k++;
            }
            lcp[rank_[i]] = k;
            if (k > 0) k--;
        }
    }
    
    // Check if pattern exists (binary search on SA)
    bool contains(const std::string& pattern) const {
        int lo = 0, hi = sa.size() - 1;
        while (lo <= hi) {
            int mid = (lo + hi) / 2;
            std::string suffix = text.substr(sa[mid], pattern.size());
            if (suffix == pattern) return true;
            if (suffix < pattern) lo = mid + 1;
            else hi = mid - 1;
        }
        return false;
    }
    
    // Find all occurrences of pattern
    std::vector<int> findAll(const std::string& pattern) const {
        std::vector<int> positions;
        int n = sa.size();
        int m = pattern.size();
        
        // Find lower bound
        int lo = 0, hi = n - 1;
        int left = n;
        while (lo <= hi) {
            int mid = (lo + hi) / 2;
            std::string suffix = text.substr(sa[mid], std::min((int)text.size() - sa[mid], m));
            if (suffix >= pattern) {
                left = mid;
                hi = mid - 1;
            } else {
                lo = mid + 1;
            }
        }
        
        // Find upper bound
        lo = 0; hi = n - 1;
        int right = -1;
        while (lo <= hi) {
            int mid = (lo + hi) / 2;
            std::string suffix = text.substr(sa[mid], std::min((int)text.size() - sa[mid], m));
            if (suffix <= pattern && suffix.substr(0, pattern.size()) == pattern) {
                right = mid;
                lo = mid + 1;
            } else if (suffix < pattern) {
                lo = mid + 1;
            } else {
                hi = mid - 1;
            }
        }
        
        // Collect all matches
        for (int i = left; i < n; i++) {
            if (text.compare(sa[i], m, pattern) == 0) {
                positions.push_back(sa[i]);
            } else if (text.compare(sa[i], std::min(m, (int)text.size() - sa[i]), pattern.substr(0, std::min(m, (int)text.size() - sa[i]))) > 0) {
                break;
            }
        }
        return positions;
    }
    
    // Query multiple motifs
    void queryMotifs(const std::vector<std::string>& motifs) const {
        std::cout << "Motif membership queries:\n";
        for (const auto& motif : motifs) {
            auto positions = findAll(motif);
            std::cout << "  \"" << motif << "\": ";
            if (positions.empty()) {
                std::cout << "NOT FOUND\n";
            } else {
                std::cout << "FOUND at positions: ";
                for (int p : positions) std::cout << p << " ";
                std::cout << "\n";
            }
        }
    }
    
    void print() const {
        std::cout << "Suffix Array for \"" << text.substr(0, std::min(50, (int)text.size())) 
                  << (text.size() > 50 ? "..." : "") << "\":\n";
        std::cout << std::string(60, '-') << "\n";
        std::cout << "i\tSA[i]\tLCP[i]\tSuffix\n";
        std::cout << std::string(60, '-') << "\n";
        for (size_t i = 0; i < sa.size() && i < 20; i++) {
            std::string suffix = text.substr(sa[i]);
            if (suffix.size() > 30) suffix = suffix.substr(0, 30) + "...";
            std::cout << i << "\t" << sa[i] << "\t" << lcp[i] << "\t\"" << suffix << "\"\n";
        }
        if (sa.size() > 20) std::cout << "... (" << sa.size() - 20 << " more entries)\n";
    }
};

// ============================================================================
// 2.c LONGEST COMMON SUBSTRING (LCS)
// ============================================================================
// Method 1: Using Suffix Array + LCP on A#B (O((n+m) log(n+m)))
// Method 2: Dynamic Programming baseline (O(n*m))
// ============================================================================

class LongestCommonSubstring {
public:
    // Method 1: Using Suffix Array with LCP
    // Build SA for A#B, then find max LCP where adjacent suffixes come from different strings
    static std::pair<std::string, int> usingSuffixArray(const std::string& A, const std::string& B) {
        if (A.empty() || B.empty()) return {"", -1};
        
        // Use a sentinel character not in alphabet
        char sentinel = '#';  // Assuming # not in DNA alphabet
        std::string combined = A + sentinel + B;
        
        SuffixArray sa(combined);
        
        int lenA = A.size();
        int maxLen = 0;
        int startPos = -1;
        
        // Find max LCP between suffixes from different original strings
        for (size_t i = 1; i < sa.sa.size(); i++) {
            int pos1 = sa.sa[i - 1];
            int pos2 = sa.sa[i];
            
            // Check if suffixes are from different strings
            bool from1 = pos1 < lenA;
            bool from2 = pos2 < lenA;
            
            if (from1 != from2) {  // From different strings
                // LCP should not cross the sentinel
                int lcpVal = sa.lcp[i];
                int maxPossible1 = (pos1 < lenA) ? (lenA - pos1) : ((int)combined.size() - pos1);
                int maxPossible2 = (pos2 < lenA) ? (lenA - pos2) : ((int)combined.size() - pos2);
                lcpVal = std::min(lcpVal, std::min(maxPossible1, maxPossible2));
                
                if (lcpVal > maxLen) {
                    maxLen = lcpVal;
                    startPos = from1 ? pos1 : pos2 - lenA - 1;
                }
            }
        }
        
        if (maxLen == 0) return {"", -1};
        
        std::string lcs = (startPos < lenA) ? A.substr(startPos, maxLen) 
                                            : B.substr(startPos, maxLen);
        return {lcs, startPos};
    }
    
    // Method 2: Dynamic Programming O(|A| * |B|)
    // dp[i][j] = length of longest common suffix ending at A[i-1] and B[j-1]
    static std::pair<std::string, int> usingDP(const std::string& A, const std::string& B) {
        if (A.empty() || B.empty()) return {"", -1};
        
        int n = A.size();
        int m = B.size();
        
        // Space-optimized: only need previous row
        std::vector<int> prev(m + 1, 0), curr(m + 1, 0);
        
        int maxLen = 0;
        int endPosA = 0;
        
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= m; j++) {
                if (A[i - 1] == B[j - 1]) {
                    curr[j] = prev[j - 1] + 1;
                    if (curr[j] > maxLen) {
                        maxLen = curr[j];
                        endPosA = i;
                    }
                } else {
                    curr[j] = 0;
                }
            }
            std::swap(prev, curr);
            std::fill(curr.begin(), curr.end(), 0);
        }
        
        if (maxLen == 0) return {"", -1};
        
        int startPos = endPosA - maxLen;
        return {A.substr(startPos, maxLen), startPos};
    }
    
    // Full DP with backtracking (returns all LCS if multiple exist)
    static std::vector<std::string> allLCS_DP(const std::string& A, const std::string& B) {
        int n = A.size();
        int m = B.size();
        
        std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1, 0));
        int maxLen = 0;
        
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= m; j++) {
                if (A[i - 1] == B[j - 1]) {
                    dp[i][j] = dp[i - 1][j - 1] + 1;
                    maxLen = std::max(maxLen, dp[i][j]);
                }
            }
        }
        
        // Collect all LCS
        std::unordered_map<std::string, bool> seen;
        std::vector<std::string> results;
        
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= m; j++) {
                if (dp[i][j] == maxLen) {
                    std::string lcs = A.substr(i - maxLen, maxLen);
                    if (!seen[lcs]) {
                        seen[lcs] = true;
                        results.push_back(lcs);
                    }
                }
            }
        }
        return results;
    }
};

// ============================================================================
// DNA ALGORITHMS DEMONSTRATION
// ============================================================================

void demonstrateDNAAlgorithms() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  DNA SEQUENCE ANALYSIS ALGORITHMS\n";
    std::cout << "  Problem 2: Genomic String Analysis\n";
    std::cout << std::string(70, '=') << "\n";
    
    // ========================================================================
    // INPUT: Two DNA strings A, B ∈ {A,C,G,T}*
    // ========================================================================
    std::string seqA = "ACGTACGTGACGTACGTAATCGAATCGACGTACGTACGTTATAACGTACGT";
    std::string seqB = "GCTAGCTAGCTACGTACGTTTTACGTACGTAGCTAGCTATAAGCTAGCTAT";
    
    std::cout << "\n=== INPUT DNA SEQUENCES ===\n";
    std::cout << "Alphabet: Σ = {A, C, G, T}\n\n";
    std::cout << "Sequence A (|A| = " << seqA.size() << "):\n";
    std::cout << "  " << seqA << "\n\n";
    std::cout << "Sequence B (|B| = " << seqB.size() << "):\n";
    std::cout << "  " << seqB << "\n";
    
    // ========================================================================
    // 2.a MANACHER'S ALGORITHM - Longest Palindromic Substring
    // ========================================================================
    std::cout << "\n" << std::string(70, '-') << "\n";
    std::cout << "2.a. MANACHER'S ALGORITHM - Longest Palindromic Regions\n";
    std::cout << std::string(70, '-') << "\n";
    
    auto [palA, posA] = Manacher::longestPalindrome(seqA);
    std::cout << "\nSequence A:\n";
    std::cout << "  Longest palindrome: \"" << palA << "\"\n";
    std::cout << "  Position: " << posA << ", Length: " << palA.size() << "\n";
    
    auto [palB, posB] = Manacher::longestPalindrome(seqB);
    std::cout << "\nSequence B:\n";
    std::cout << "  Longest palindrome: \"" << palB << "\"\n";
    std::cout << "  Position: " << posB << ", Length: " << palB.size() << "\n";
    
    std::cout << "\n[Radius Propagation Example]\n";
    std::cout << "String: \"ACGTGCA\" (palindrome ACGTGCA with center at T)\n";
    Manacher::printRadiiExplanation("ACGTGCA");
    
    std::cout << "\n[CORRECTNESS JUSTIFICATION - Radius Propagation]\n";
    std::cout << "The algorithm maintains the rightmost palindrome boundary [l, r].\n";
    std::cout << "For each position i:\n";
    std::cout << "  1. If i < r: Use mirror position i' = 2*center - i\n";
    std::cout << "     Initial radius = min(p[i'], r - i) due to symmetry\n";
    std::cout << "  2. Expand beyond known bounds if possible\n";
    std::cout << "  3. Update [l, r] if expansion extends past r\n";
    std::cout << "\n";
    std::cout << "PROOF OF O(n):\n";
    std::cout << "  - Right boundary r only moves right (never left)\n";
    std::cout << "  - Each character is involved in at most 2 comparisons:\n";
    std::cout << "    once during expansion, once through mirroring\n";
    std::cout << "  - Total comparisons ≤ 2n → O(n) time complexity\n";
    
    // ========================================================================
    // 2.b SUFFIX ARRAY with KASAI LCP - Motif Queries
    // ========================================================================
    std::cout << "\n" << std::string(70, '-') << "\n";
    std::cout << "2.b. SUFFIX ARRAY + KASAI LCP - Exact Substring Queries\n";
    std::cout << std::string(70, '-') << "\n";
    
    // Build suffix array for A#B (# is unique sentinel)
    std::string combined = seqA + "#" + seqB;
    SuffixArray sa(combined);
    
    std::cout << "\nBuilt Suffix Array for A#B with unique sentinel '#'\n";
    std::cout << "Combined length: " << combined.size() << "\n";
    sa.print();
    
    // MOTIF SET M = {M1, ..., Mk} with lengths 4-12
    std::cout << "\n[MOTIF SET M - Exact Substring Membership Queries]\n";
    std::cout << "Motif lengths: 4-12 bp (as specified)\n\n";
    
    std::vector<std::string> motifs = {
        "ACGT",         // M1: 4 bp - common tetranucleotide
        "GAATTC",       // M2: 6 bp - EcoRI restriction site
        "GCTA",         // M3: 4 bp - common pattern
        "AATCGAATCG",   // M4: 10 bp - longer motif
        "TACGTACGT",    // M5: 9 bp - repeated unit
        "TATA",         // M6: 4 bp - TATA box element
        "TTTTACGT",     // M7: 8 bp 
        "XXXXXX"        // M8: 6 bp - should NOT be found
    };
    
    std::cout << "Motif membership queries (using binary search on SA):\n";
    std::cout << std::string(50, '-') << "\n";
    for (size_t i = 0; i < motifs.size(); i++) {
        const auto& motif = motifs[i];
        auto positions = sa.findAll(motif);
        std::cout << "M" << i+1 << " \"" << motif << "\" (len=" << motif.size() << "): ";
        if (positions.empty()) {
            std::cout << "NOT FOUND\n";
        } else {
            std::cout << "FOUND at " << positions.size() << " position(s): ";
            for (size_t j = 0; j < std::min(positions.size(), (size_t)5); j++) {
                std::cout << positions[j];
                if (j < positions.size() - 1 && j < 4) std::cout << ", ";
            }
            if (positions.size() > 5) std::cout << "...";
            std::cout << "\n";
        }
    }
    
    std::cout << "\nSuffix Array Query Complexity: O(m log n) per query\n";
    std::cout << "  where m = motif length, n = |A| + |B| + 1\n";
    
    // ========================================================================
    // 2.c LONGEST COMMON SUBSTRING (LCSsubstr)
    // ========================================================================
    std::cout << "\n" << std::string(70, '-') << "\n";
    std::cout << "2.c. LONGEST COMMON SUBSTRING (LCSsubstr)\n";
    std::cout << std::string(70, '-') << "\n";
    
    std::cout << "\n[Method 1: Suffix Array + LCP]\n";
    std::cout << "Approach: Find max LCP where adjacent suffixes come from different strings\n";
    auto [lcsSA, posSA] = LongestCommonSubstring::usingSuffixArray(seqA, seqB);
    std::cout << "Result: \"" << lcsSA << "\"\n";
    std::cout << "Length: " << lcsSA.size() << " bp\n";
    std::cout << "Complexity: O((|A|+|B|) log(|A|+|B|)) for SA construction + O(|A|+|B|) LCP\n";
    
    std::cout << "\n[Method 2: Dynamic Programming - O(|A|×|B|) baseline]\n";
    std::cout << "Approach: dp[i][j] = length of LCS ending at A[i-1] and B[j-1]\n";
    auto [lcsDP, posDP] = LongestCommonSubstring::usingDP(seqA, seqB);
    std::cout << "Result: \"" << lcsDP << "\"\n";
    std::cout << "Length: " << lcsDP.size() << " bp\n";
    std::cout << "Complexity: O(|A|×|B|) time, O(min(|A|,|B|)) space (optimized)\n";
    
    // Comparison
    std::cout << "\n[METHODOLOGICAL COMPARISON]\n";
    if (lcsSA.size() == lcsDP.size()) {
        std::cout << "✓ Both methods found LCS of same length: " << lcsSA.size() << "\n";
    }
    
    auto allLCS = LongestCommonSubstring::allLCS_DP(seqA, seqB);
    std::cout << "\nAll longest common substrings found:\n";
    for (size_t i = 0; i < allLCS.size(); i++) {
        std::cout << "  " << i+1 << ". \"" << allLCS[i] << "\"\n";
    }
    
    std::cout << "\n[COMPLEXITY SUMMARY]\n";
    std::cout << "Suffix Array method: O(n log n) construction + O(n) LCP + O(n) scan\n";
    std::cout << "DP method:          O(|A|×|B|) time, better for short sequences\n";
    std::cout << "SA preferred when: |A|,|B| are large and SA reused for multiple queries\n";
}

#endif // DNA_ALGORITHMS_H
