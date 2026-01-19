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
    std::cout << std::string(70, '=') << "\n";
    
    // Sample DNA sequences
    std::string seqA = "ACGTACGTGACGTACGTAATCGAATCGACGTACGTACGT";
    std::string seqB = "GCTAGCTAGCTACGTACGTTTTACGTACGTAGCTAGCTA";
    
    std::cout << "\nSequence A: " << seqA << " (length: " << seqA.size() << ")\n";
    std::cout << "Sequence B: " << seqB << " (length: " << seqB.size() << ")\n";
    
    // 2.a Manacher's Algorithm - Longest Palindrome
    std::cout << "\n--- 2.a Manacher's Algorithm (Longest Palindrome) ---\n";
    
    auto [palA, posA] = Manacher::longestPalindrome(seqA);
    std::cout << "Sequence A longest palindrome: \"" << palA << "\" at position " << posA 
              << " (length: " << palA.size() << ")\n";
    
    auto [palB, posB] = Manacher::longestPalindrome(seqB);
    std::cout << "Sequence B longest palindrome: \"" << palB << "\" at position " << posB 
              << " (length: " << palB.size() << ")\n";
    
    // Detailed radius propagation example
    std::cout << "\nRadius propagation example (small string):\n";
    Manacher::printRadiiExplanation("ACGTGCA");
    
    std::cout << "\nCorrectness Justification:\n";
    std::cout << "- The algorithm exploits symmetry: if position i is within a known\n";
    std::cout << "  palindrome [l,r], its mirror i' = 2*center - i gives initial radius.\n";
    std::cout << "- Expansion only occurs past previously checked boundaries.\n";
    std::cout << "- Each character participates in at most 2 comparisons → O(n) total.\n";
    
    // 2.b Suffix Array with Kasai LCP
    std::cout << "\n--- 2.b Suffix Array + Kasai LCP ---\n";
    
    std::string combined = seqA + "#" + seqB;
    SuffixArray sa(combined);
    
    std::cout << "Built suffix array for A#B (length: " << combined.size() << ")\n";
    sa.print();
    
    // Motif queries
    std::vector<std::string> motifs = {"ACGT", "GCTA", "TTTT", "AATCG", "XXXX", "TACGTACGT"};
    std::cout << "\nExact motif membership queries:\n";
    sa.queryMotifs(motifs);
    
    // 2.c Longest Common Substring
    std::cout << "\n--- 2.c Longest Common Substring ---\n";
    
    std::cout << "\nMethod 1: Suffix Array + LCP\n";
    auto [lcsSA, posSA] = LongestCommonSubstring::usingSuffixArray(seqA, seqB);
    std::cout << "LCS: \"" << lcsSA << "\" (length: " << lcsSA.size() << ")\n";
    
    std::cout << "\nMethod 2: Dynamic Programming O(|A|*|B|)\n";
    auto [lcsDP, posDP] = LongestCommonSubstring::usingDP(seqA, seqB);
    std::cout << "LCS: \"" << lcsDP << "\" (length: " << lcsDP.size() << ")\n";
    
    // Verify both methods give same result
    if (lcsSA == lcsDP) {
        std::cout << "\n✓ Both methods produce the same LCS.\n";
    } else {
        std::cout << "\n! Methods produced different results (both valid if same length).\n";
    }
    
    // Find all LCS (if multiple)
    auto allLCS = LongestCommonSubstring::allLCS_DP(seqA, seqB);
    if (allLCS.size() > 1) {
        std::cout << "All LCS found: ";
        for (const auto& s : allLCS) std::cout << "\"" << s << "\" ";
        std::cout << "\n";
    }
}

#endif // DNA_ALGORITHMS_H
