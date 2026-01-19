#ifndef EXPRESSION_PARSER_H
#define EXPRESSION_PARSER_H

#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <stack>
#include <memory>
#include <cctype>
#include <algorithm>
#include <iomanip>
#include "string_algorithms.h"

// ============================================================================
// COMBINED SYSTEM: LOGICAL EXPRESSION PARSER
// ============================================================================
// Features:
// 1.a. Accepts logical expressions with up to 3 propositional variables
//      Operators: NOT, AND, OR, XOR with parentheses support
// 1.b. Detects Polish (prefix) vs Infix notation using string matching (KMP/Z)
// 1.c. Parses to AST and evaluates; produces truth tables
// ============================================================================

// Forward declarations
class ASTNode;
using ASTPtr = std::shared_ptr<ASTNode>;

// ============================================================================
// TOKENIZER
// ============================================================================

enum class TokenType {
    VAR,        // Variable (P, Q, R, etc.)
    NOT,        // NOT operator
    AND,        // AND operator
    OR,         // OR operator
    XOR,        // XOR operator
    LPAREN,     // (
    RPAREN,     // )
    END         // End of input
};

struct Token {
    TokenType type;
    std::string value;
    
    std::string toString() const {
        switch(type) {
            case TokenType::VAR: return "VAR(" + value + ")";
            case TokenType::NOT: return "NOT";
            case TokenType::AND: return "AND";
            case TokenType::OR: return "OR";
            case TokenType::XOR: return "XOR";
            case TokenType::LPAREN: return "LPAREN";
            case TokenType::RPAREN: return "RPAREN";
            case TokenType::END: return "END";
        }
        return "UNKNOWN";
    }
};

class Tokenizer {
public:
    static std::vector<Token> tokenize(const std::string& expr) {
        std::vector<Token> tokens;
        std::string upper = toUpper(expr);
        size_t i = 0;
        
        while (i < upper.size()) {
            if (std::isspace(upper[i])) {
                i++;
                continue;
            }
            
            if (upper[i] == '(') {
                tokens.push_back({TokenType::LPAREN, "("});
                i++;
            } else if (upper[i] == ')') {
                tokens.push_back({TokenType::RPAREN, ")"});
                i++;
            } else if (upper.substr(i, 3) == "NOT" || upper[i] == '!') {
                tokens.push_back({TokenType::NOT, "NOT"});
                i += (upper[i] == '!' ? 1 : 3);
            } else if (upper.substr(i, 3) == "AND" || upper[i] == '&') {
                tokens.push_back({TokenType::AND, "AND"});
                i += (upper[i] == '&' ? 1 : 3);
            } else if (upper.substr(i, 2) == "OR" || upper[i] == '|') {
                tokens.push_back({TokenType::OR, "OR"});
                i += (upper[i] == '|' ? 1 : 2);
            } else if (upper.substr(i, 3) == "XOR" || upper[i] == '^') {
                tokens.push_back({TokenType::XOR, "XOR"});
                i += (upper[i] == '^' ? 1 : 3);
            } else if (std::isalpha(upper[i])) {
                // Variable name (single letter)
                tokens.push_back({TokenType::VAR, std::string(1, upper[i])});
                i++;
            } else {
                // Skip unknown characters
                i++;
            }
        }
        
        tokens.push_back({TokenType::END, ""});
        return tokens;
    }
    
private:
    static std::string toUpper(const std::string& s) {
        std::string result = s;
        for (char& c : result) c = std::toupper(c);
        return result;
    }
};

// ============================================================================
// NOTATION DETECTION (Using KMP/Z-Algorithm)
// ============================================================================

class NotationDetector {
public:
    enum Notation { INFIX, PREFIX, INVALID };
    
    // Detect notation using string pattern matching
    // Polish (prefix): operator appears before operands, e.g., "AND P Q"
    // Infix: operator appears between operands, e.g., "P AND Q"
    static Notation detect(const std::string& expr) {
        std::string normalized = normalizeSpaces(expr);
        
        // Patterns that indicate prefix notation
        // Prefix: expression starts with operator followed by space
        std::vector<std::string> prefixPatterns = {
            "NOT ", "AND ", "OR ", "XOR "
        };
        
        // Patterns that indicate infix notation
        // Infix: variable/paren followed by operator
        std::vector<std::string> infixPatterns = {
            " AND ", " OR ", " XOR "
        };
        
        bool hasPrefix = false;
        bool hasInfix = false;
        
        // Use KMP to search for patterns
        for (const auto& pattern : prefixPatterns) {
            auto matches = KMP::search(normalized, pattern);
            // Check if pattern appears at start (after optional parentheses)
            for (int pos : matches) {
                if (pos == 0 || (pos > 0 && normalized[pos-1] == '(')) {
                    hasPrefix = true;
                    break;
                }
            }
        }
        
        for (const auto& pattern : infixPatterns) {
            auto matches = ZAlgorithm::search(normalized, pattern);
            if (!matches.empty()) {
                hasInfix = true;
                break;
            }
        }
        
        // Additional heuristics
        std::vector<Token> tokens = Tokenizer::tokenize(expr);
        
        // Check first non-paren token
        int firstIdx = 0;
        while (firstIdx < (int)tokens.size() && tokens[firstIdx].type == TokenType::LPAREN) {
            firstIdx++;
        }
        
        if (firstIdx < (int)tokens.size()) {
            TokenType first = tokens[firstIdx].type;
            
            // If first token (after parens) is a binary operator → likely PREFIX
            if (first == TokenType::AND || first == TokenType::OR || first == TokenType::XOR) {
                return PREFIX;
            }
            
            // If first token is variable, check what follows
            if (first == TokenType::VAR) {
                // Look for binary operator after variable (or after NOT VAR)
                int nextIdx = firstIdx + 1;
                if (nextIdx < (int)tokens.size() - 1) {
                    TokenType next = tokens[nextIdx].type;
                    if (next == TokenType::AND || next == TokenType::OR || next == TokenType::XOR) {
                        return INFIX;
                    }
                }
            }
            
            // If first token is NOT, could be either
            if (first == TokenType::NOT) {
                // "NOT P AND Q" could be prefix or infix
                // Check if there's "VAR OP VAR" pattern for infix
                for (size_t i = 0; i < tokens.size() - 2; i++) {
                    if (tokens[i].type == TokenType::VAR && 
                        (tokens[i+1].type == TokenType::AND || 
                         tokens[i+1].type == TokenType::OR ||
                         tokens[i+1].type == TokenType::XOR) &&
                        tokens[i+2].type == TokenType::VAR) {
                        return INFIX;
                    }
                }
                // Otherwise assume prefix
                return PREFIX;
            }
        }
        
        // Default based on pattern matching results
        if (hasInfix && !hasPrefix) return INFIX;
        if (hasPrefix && !hasInfix) return PREFIX;
        if (hasInfix) return INFIX;  // Preference for infix in ambiguous cases
        
        return INVALID;
    }
    
    static std::string notationName(Notation n) {
        switch(n) {
            case INFIX: return "INFIX (normal notation)";
            case PREFIX: return "PREFIX (Polish notation)";
            case INVALID: return "INVALID/UNKNOWN";
        }
        return "UNKNOWN";
    }
    
private:
    static std::string normalizeSpaces(const std::string& s) {
        std::string result;
        bool lastSpace = true;
        for (char c : s) {
            if (std::isspace(c)) {
                if (!lastSpace) {
                    result += ' ';
                    lastSpace = true;
                }
            } else {
                result += std::toupper(c);
                lastSpace = false;
            }
        }
        if (!result.empty() && result.back() == ' ') {
            result.pop_back();
        }
        return result;
    }
};

// ============================================================================
// ABSTRACT SYNTAX TREE (AST)
// ============================================================================

class ASTNode {
public:
    virtual ~ASTNode() = default;
    virtual bool evaluate(const std::unordered_map<char, bool>& vars) const = 0;
    virtual void print(int depth = 0) const = 0;
    virtual std::string toString() const = 0;
};

class VarNode : public ASTNode {
    char name;
public:
    VarNode(char n) : name(n) {}
    
    bool evaluate(const std::unordered_map<char, bool>& vars) const override {
        auto it = vars.find(name);
        if (it != vars.end()) return it->second;
        return false;
    }
    
    void print(int depth = 0) const override {
        std::cout << std::string(depth * 2, ' ') << name << "\n";
    }
    
    std::string toString() const override {
        return std::string(1, name);
    }
};

class NotNode : public ASTNode {
    ASTPtr child;
public:
    NotNode(ASTPtr c) : child(c) {}
    
    bool evaluate(const std::unordered_map<char, bool>& vars) const override {
        return !child->evaluate(vars);
    }
    
    void print(int depth = 0) const override {
        std::cout << std::string(depth * 2, ' ') << "NOT\n";
        child->print(depth + 1);
    }
    
    std::string toString() const override {
        return "NOT(" + child->toString() + ")";
    }
};

class BinaryNode : public ASTNode {
protected:
    ASTPtr left, right;
    std::string op;
public:
    BinaryNode(ASTPtr l, ASTPtr r, const std::string& o) : left(l), right(r), op(o) {}
    
    void print(int depth = 0) const override {
        std::cout << std::string(depth * 2, ' ') << op << "\n";
        left->print(depth + 1);
        right->print(depth + 1);
    }
    
    std::string toString() const override {
        return "(" + left->toString() + " " + op + " " + right->toString() + ")";
    }
};

class AndNode : public BinaryNode {
public:
    AndNode(ASTPtr l, ASTPtr r) : BinaryNode(l, r, "AND") {}
    bool evaluate(const std::unordered_map<char, bool>& vars) const override {
        return left->evaluate(vars) && right->evaluate(vars);
    }
};

class OrNode : public BinaryNode {
public:
    OrNode(ASTPtr l, ASTPtr r) : BinaryNode(l, r, "OR") {}
    bool evaluate(const std::unordered_map<char, bool>& vars) const override {
        return left->evaluate(vars) || right->evaluate(vars);
    }
};

class XorNode : public BinaryNode {
public:
    XorNode(ASTPtr l, ASTPtr r) : BinaryNode(l, r, "XOR") {}
    bool evaluate(const std::unordered_map<char, bool>& vars) const override {
        return left->evaluate(vars) != right->evaluate(vars);
    }
};

// ============================================================================
// PARSERS
// ============================================================================

// Infix Parser (using precedence: NOT > AND > XOR > OR)
class InfixParser {
    std::vector<Token> tokens;
    size_t pos;
    
    Token current() const {
        return pos < tokens.size() ? tokens[pos] : Token{TokenType::END, ""};
    }
    
    Token consume() {
        return tokens[pos++];
    }
    
    bool match(TokenType type) {
        if (current().type == type) {
            consume();
            return true;
        }
        return false;
    }
    
    // Grammar:
    // expr    -> orExpr
    // orExpr  -> xorExpr (OR xorExpr)*
    // xorExpr -> andExpr (XOR andExpr)*
    // andExpr -> notExpr (AND notExpr)*
    // notExpr -> NOT notExpr | primary
    // primary -> VAR | LPAREN expr RPAREN
    
    ASTPtr parseOr() {
        ASTPtr left = parseXor();
        while (match(TokenType::OR)) {
            ASTPtr right = parseXor();
            left = std::make_shared<OrNode>(left, right);
        }
        return left;
    }
    
    ASTPtr parseXor() {
        ASTPtr left = parseAnd();
        while (match(TokenType::XOR)) {
            ASTPtr right = parseAnd();
            left = std::make_shared<XorNode>(left, right);
        }
        return left;
    }
    
    ASTPtr parseAnd() {
        ASTPtr left = parseNot();
        while (match(TokenType::AND)) {
            ASTPtr right = parseNot();
            left = std::make_shared<AndNode>(left, right);
        }
        return left;
    }
    
    ASTPtr parseNot() {
        if (match(TokenType::NOT)) {
            return std::make_shared<NotNode>(parseNot());
        }
        return parsePrimary();
    }
    
    ASTPtr parsePrimary() {
        if (match(TokenType::LPAREN)) {
            ASTPtr expr = parseOr();
            match(TokenType::RPAREN);  // Consume closing paren
            return expr;
        }
        if (current().type == TokenType::VAR) {
            char name = current().value[0];
            consume();
            return std::make_shared<VarNode>(name);
        }
        // Error: unexpected token
        return nullptr;
    }
    
public:
    ASTPtr parse(const std::string& expr) {
        tokens = Tokenizer::tokenize(expr);
        pos = 0;
        return parseOr();
    }
};

// Prefix Parser (Polish notation)
class PrefixParser {
    std::vector<Token> tokens;
    size_t pos;
    
    Token current() const {
        return pos < tokens.size() ? tokens[pos] : Token{TokenType::END, ""};
    }
    
    Token consume() {
        return tokens[pos++];
    }
    
    // Grammar:
    // expr -> VAR | NOT expr | AND expr expr | OR expr expr | XOR expr expr | LPAREN expr RPAREN
    
    ASTPtr parseExpr() {
        Token tok = current();
        
        if (tok.type == TokenType::VAR) {
            consume();
            return std::make_shared<VarNode>(tok.value[0]);
        }
        
        if (tok.type == TokenType::NOT) {
            consume();
            return std::make_shared<NotNode>(parseExpr());
        }
        
        if (tok.type == TokenType::AND) {
            consume();
            ASTPtr left = parseExpr();
            ASTPtr right = parseExpr();
            return std::make_shared<AndNode>(left, right);
        }
        
        if (tok.type == TokenType::OR) {
            consume();
            ASTPtr left = parseExpr();
            ASTPtr right = parseExpr();
            return std::make_shared<OrNode>(left, right);
        }
        
        if (tok.type == TokenType::XOR) {
            consume();
            ASTPtr left = parseExpr();
            ASTPtr right = parseExpr();
            return std::make_shared<XorNode>(left, right);
        }
        
        if (tok.type == TokenType::LPAREN) {
            consume();
            ASTPtr expr = parseExpr();
            if (current().type == TokenType::RPAREN) consume();
            return expr;
        }
        
        return nullptr;
    }
    
public:
    ASTPtr parse(const std::string& expr) {
        tokens = Tokenizer::tokenize(expr);
        pos = 0;
        return parseExpr();
    }
};

// ============================================================================
// EXPRESSION EVALUATOR WITH TRUTH TABLE
// ============================================================================

class ExpressionEvaluator {
public:
    // Parse and evaluate expression, generating truth table
    static void evaluate(const std::string& expr, 
                        const std::vector<char>& varNames,
                        bool verbose = true) {
        
        if (verbose) {
            std::cout << "\n" << std::string(60, '-') << "\n";
            std::cout << "Expression: " << expr << "\n";
        }
        
        // Detect notation
        auto notation = NotationDetector::detect(expr);
        if (verbose) {
            std::cout << "Detected notation: " << NotationDetector::notationName(notation) << "\n";
        }
        
        // Parse to AST
        ASTPtr ast;
        if (notation == NotationDetector::PREFIX) {
            PrefixParser parser;
            ast = parser.parse(expr);
        } else {
            InfixParser parser;
            ast = parser.parse(expr);
        }
        
        if (!ast) {
            std::cout << "ERROR: Failed to parse expression!\n";
            return;
        }
        
        if (verbose) {
            std::cout << "\nAbstract Syntax Tree:\n";
            ast->print(1);
            std::cout << "\nCanonical form: " << ast->toString() << "\n";
        }
        
        // Generate truth table
        if (verbose) {
            std::cout << "\nTruth Table:\n";
            printTruthTable(ast, varNames);
        }
    }
    
    static void printTruthTable(ASTPtr ast, const std::vector<char>& varNames) {
        int n = varNames.size();
        int rows = 1 << n;
        
        // Header
        for (char v : varNames) {
            std::cout << std::setw(3) << v << " ";
        }
        std::cout << "| Result\n";
        std::cout << std::string(4 * n + 9, '-') << "\n";
        
        // Rows
        for (int i = 0; i < rows; i++) {
            std::unordered_map<char, bool> vars;
            for (int j = 0; j < n; j++) {
                bool val = (i >> (n - 1 - j)) & 1;
                vars[varNames[j]] = val;
                std::cout << std::setw(3) << val << " ";
            }
            
            bool result = ast->evaluate(vars);
            std::cout << "|   " << result << "\n";
        }
    }
    
    // Evaluate single assignment
    static bool evaluateSingle(ASTPtr ast, const std::unordered_map<char, bool>& vars) {
        return ast->evaluate(vars);
    }
};

// ============================================================================
// TEST CASES
// ============================================================================

void runExpressionParserTests() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  LOGICAL EXPRESSION PARSER - TEST CASES\n";
    std::cout << std::string(70, '=') << "\n";
    
    std::vector<char> vars = {'P', 'Q', 'R'};
    
    // INFIX NOTATION TEST CASES
    std::cout << "\n### INFIX NOTATION TEST CASES ###\n";
    
    std::vector<std::string> infixTests = {
        "P AND Q",                    // Valid: simple AND
        "(P OR Q) AND R",             // Valid: with parentheses
        "NOT P AND (Q XOR R)",        // Valid: mixed operators
        "P OR Q OR R",                // Valid: chained OR
        "NOT (P AND Q)",              // Valid: NOT with paren
        "P XOR NOT Q",                // Valid: XOR with NOT
        "((P))",                      // Valid: nested parens
        "AND P Q"                     // Invalid for infix (detected as prefix)
    };
    
    std::cout << "\n--- Valid Infix Expressions ---\n";
    for (size_t i = 0; i < infixTests.size() - 1; i++) {
        std::cout << "\nTest case " << i + 1 << ":\n";
        ExpressionEvaluator::evaluate(infixTests[i], vars);
    }
    
    std::cout << "\n--- Invalid/Ambiguous Infix Expression ---\n";
    std::cout << "\nTest case (should detect as prefix):\n";
    ExpressionEvaluator::evaluate(infixTests.back(), vars);
    
    // PREFIX NOTATION TEST CASES
    std::cout << "\n\n### PREFIX (POLISH) NOTATION TEST CASES ###\n";
    
    std::vector<std::string> prefixTests = {
        "AND P Q",                    // Valid: simple AND
        "OR AND P Q R",               // Valid: nested
        "NOT AND P Q",                // Valid: NOT of AND
        "XOR P NOT Q",                // Valid: XOR with NOT
        "AND OR P Q NOT R",           // Valid: complex
        "NOT NOT P",                  // Valid: double NOT
        "P AND Q",                    // Invalid for prefix (detected as infix)
    };
    
    std::cout << "\n--- Valid Prefix Expressions ---\n";
    for (size_t i = 0; i < prefixTests.size() - 1; i++) {
        std::cout << "\nTest case " << i + 1 << ":\n";
        ExpressionEvaluator::evaluate(prefixTests[i], vars);
    }
    
    std::cout << "\n--- Invalid/Ambiguous Prefix Expression ---\n";
    std::cout << "\nTest case (should detect as infix):\n";
    ExpressionEvaluator::evaluate(prefixTests.back(), vars);
    
    // AST EXAMPLES
    std::cout << "\n\n### AST EXAMPLES ###\n";
    
    std::cout << "\n--- Infix Example 1: (P AND Q) OR (NOT R) ---\n";
    ExpressionEvaluator::evaluate("(P AND Q) OR (NOT R)", vars);
    
    std::cout << "\n--- Infix Example 2: P XOR (Q AND R) ---\n";
    ExpressionEvaluator::evaluate("P XOR (Q AND R)", vars);
    
    std::cout << "\n--- Prefix Example 1: OR AND P Q NOT R ---\n";
    ExpressionEvaluator::evaluate("OR AND P Q NOT R", vars);
    
    std::cout << "\n--- Prefix Example 2: XOR P AND Q R ---\n";
    ExpressionEvaluator::evaluate("XOR P AND Q R", vars);
    
    // Complexity Discussion
    std::cout << "\n\n### COMPLEXITY ANALYSIS ###\n";
    std::cout << std::string(60, '-') << "\n";
    std::cout << "PARSING COMPLEXITY:\n";
    std::cout << "  Tokenization: O(n) where n = expression length\n";
    std::cout << "  Notation detection: O(n*k) using KMP/Z, k = pattern count\n";
    std::cout << "  Infix parsing: O(n) - each token visited once\n";
    std::cout << "  Prefix parsing: O(n) - each token visited once\n";
    std::cout << "  Total: O(n)\n\n";
    std::cout << "EVALUATION COMPLEXITY:\n";
    std::cout << "  Single evaluation: O(m) where m = AST nodes\n";
    std::cout << "  Truth table (k vars): O(2^k * m)\n\n";
    std::cout << "CORRECTNESS:\n";
    std::cout << "  Infix parser uses recursive descent with proper precedence:\n";
    std::cout << "    NOT > AND > XOR > OR (highest to lowest)\n";
    std::cout << "  Prefix parser directly follows Polish notation semantics.\n";
    std::cout << "  AST evaluation uses post-order DFS traversal (evaluate\n";
    std::cout << "  children before parent), guaranteeing correct results.\n";
}

void demonstrateExpressionParser() {
    runExpressionParserTests();
}

#endif // EXPRESSION_PARSER_H
