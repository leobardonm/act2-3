#ifndef GRAPH_ALGORITHMS_H
#define GRAPH_ALGORITHMS_H

#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>
#include <climits>
#include <functional>
#include <set>
#include <map>
#include <iomanip>

// ============================================================================
// GRAPH UTILITY STRUCTURES
// ============================================================================

struct Edge {
    int u, v, weight;
    bool operator<(const Edge& other) const {
        return weight < other.weight;
    }
};

// ============================================================================
// 1. MAXIMUM FLOW - EDMONDS-KARP ALGORITHM
// ============================================================================
// Time Complexity: O(V * E^2)
// Uses BFS to find shortest augmenting paths (in terms of edges)
// ============================================================================

class EdmondsKarp {
private:
    int n;
    std::vector<std::vector<int>> capacity;
    std::vector<std::vector<int>> flow;
    std::vector<std::vector<int>> adj;
    
public:
    EdmondsKarp(int vertices) : n(vertices) {
        capacity.assign(n, std::vector<int>(n, 0));
        flow.assign(n, std::vector<int>(n, 0));
        adj.resize(n);
    }
    
    void addEdge(int u, int v, int cap) {
        capacity[u][v] += cap;
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    
    // BFS to find augmenting path
    bool bfs(int source, int sink, std::vector<int>& parent) {
        std::fill(parent.begin(), parent.end(), -1);
        parent[source] = source;
        std::queue<std::pair<int, int>> q;
        q.push({source, INT_MAX});
        
        while (!q.empty()) {
            int cur = q.front().first;
            int curFlow = q.front().second;
            q.pop();
            
            for (int next : adj[cur]) {
                int residual = capacity[cur][next] - flow[cur][next];
                if (parent[next] == -1 && residual > 0) {
                    parent[next] = cur;
                    int newFlow = std::min(curFlow, residual);
                    if (next == sink) return true;
                    q.push({next, newFlow});
                }
            }
        }
        return false;
    }
    
    // Compute maximum flow with step-by-step output
    int maxFlow(int source, int sink, bool verbose = true) {
        int totalFlow = 0;
        std::vector<int> parent(n);
        int iteration = 0;
        
        if (verbose) {
            std::cout << "\n--- Edmonds-Karp Algorithm Execution ---\n";
            std::cout << "Source: " << source << ", Sink: " << sink << "\n\n";
        }
        
        while (bfs(source, sink, parent)) {
            iteration++;
            
            // Find bottleneck (min residual on path)
            int pathFlow = INT_MAX;
            std::vector<int> path;
            for (int v = sink; v != source; v = parent[v]) {
                int u = parent[v];
                pathFlow = std::min(pathFlow, capacity[u][v] - flow[u][v]);
                path.push_back(v);
            }
            path.push_back(source);
            std::reverse(path.begin(), path.end());
            
            if (verbose) {
                std::cout << "Iteration " << iteration << ":\n";
                std::cout << "  Augmenting path: ";
                for (size_t i = 0; i < path.size(); i++) {
                    std::cout << path[i];
                    if (i < path.size() - 1) std::cout << " -> ";
                }
                std::cout << "\n  Bottleneck flow: " << pathFlow << "\n";
            }
            
            // Update flow along path
            for (int v = sink; v != source; v = parent[v]) {
                int u = parent[v];
                flow[u][v] += pathFlow;
                flow[v][u] -= pathFlow;
            }
            
            totalFlow += pathFlow;
            
            if (verbose) {
                std::cout << "  Total flow so far: " << totalFlow << "\n\n";
            }
        }
        
        if (verbose) {
            printResidualGraph();
            printFinalFlow();
        }
        
        return totalFlow;
    }
    
    void printResidualGraph() {
        std::cout << "Residual Graph (capacity - flow):\n";
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                int residual = capacity[i][j] - flow[i][j];
                if (capacity[i][j] > 0 || residual > 0) {
                    std::cout << "  " << i << " -> " << j << ": residual = " << residual << "\n";
                }
            }
        }
        std::cout << "\n";
    }
    
    void printFinalFlow() {
        std::cout << "Final Flow per Edge:\n";
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (capacity[i][j] > 0 && flow[i][j] > 0) {
                    std::cout << "  " << i << " -> " << j << ": " << flow[i][j] 
                              << "/" << capacity[i][j] << "\n";
                }
            }
        }
    }
    
    std::vector<std::vector<int>>& getFlow() { return flow; }
    std::vector<std::vector<int>>& getCapacity() { return capacity; }
};

// ============================================================================
// DINIC'S ALGORITHM (Alternative for comparison)
// ============================================================================
// Time Complexity: O(V^2 * E)
// Uses level graph and blocking flows for better performance on dense graphs
// ============================================================================

class Dinic {
private:
    int n;
    std::vector<std::vector<int>> capacity;
    std::vector<std::vector<int>> adj;
    std::vector<int> level, iter;
    
public:
    Dinic(int vertices) : n(vertices) {
        capacity.assign(n, std::vector<int>(n, 0));
        adj.resize(n);
        level.resize(n);
        iter.resize(n);
    }
    
    void addEdge(int u, int v, int cap) {
        capacity[u][v] += cap;
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    
    bool bfs(int source, int sink) {
        std::fill(level.begin(), level.end(), -1);
        level[source] = 0;
        std::queue<int> q;
        q.push(source);
        
        while (!q.empty()) {
            int v = q.front(); q.pop();
            for (int to : adj[v]) {
                if (capacity[v][to] > 0 && level[to] < 0) {
                    level[to] = level[v] + 1;
                    q.push(to);
                }
            }
        }
        return level[sink] >= 0;
    }
    
    int dfs(int v, int sink, int pushed) {
        if (v == sink || pushed == 0) return pushed;
        for (int& cid = iter[v]; cid < (int)adj[v].size(); cid++) {
            int to = adj[v][cid];
            if (level[v] + 1 != level[to] || capacity[v][to] <= 0)
                continue;
            int tr = dfs(to, sink, std::min(pushed, capacity[v][to]));
            if (tr > 0) {
                capacity[v][to] -= tr;
                capacity[to][v] += tr;
                return tr;
            }
        }
        return 0;
    }
    
    int maxFlow(int source, int sink) {
        int flow = 0;
        while (bfs(source, sink)) {
            std::fill(iter.begin(), iter.end(), 0);
            while (int pushed = dfs(source, sink, INT_MAX)) {
                flow += pushed;
            }
        }
        return flow;
    }
};

// ============================================================================
// 2. MINIMUM SPANNING TREE - KRUSKAL & PRIM
// ============================================================================

// Union-Find (Disjoint Set Union) for Kruskal's
class UnionFind {
public:
    std::vector<int> parent, rank_;
    
    UnionFind(int n) : parent(n), rank_(n, 0) {
        for (int i = 0; i < n; i++) parent[i] = i;
    }
    
    int find(int x) {
        if (parent[x] != x) parent[x] = find(parent[x]);
        return parent[x];
    }
    
    bool unite(int x, int y) {
        int px = find(x), py = find(y);
        if (px == py) return false;
        if (rank_[px] < rank_[py]) std::swap(px, py);
        parent[py] = px;
        if (rank_[px] == rank_[py]) rank_[px]++;
        return true;
    }
};

class MST {
public:
    // Kruskal's Algorithm: O(E log E)
    static std::pair<int, std::vector<Edge>> kruskal(int V, std::vector<Edge>& edges, bool verbose = true) {
        std::sort(edges.begin(), edges.end());
        UnionFind uf(V);
        
        std::vector<Edge> mstEdges;
        int totalWeight = 0;
        
        if (verbose) std::cout << "\n--- Kruskal's Algorithm ---\n";
        
        for (const Edge& e : edges) {
            if (uf.unite(e.u, e.v)) {
                mstEdges.push_back(e);
                totalWeight += e.weight;
                if (verbose) {
                    std::cout << "  Added edge: " << e.u << " -- " << e.v 
                              << " (weight: " << e.weight << ")\n";
                }
                if ((int)mstEdges.size() == V - 1) break;
            }
        }
        
        if (verbose) {
            std::cout << "Total MST weight: " << totalWeight << "\n";
        }
        
        return {totalWeight, mstEdges};
    }
    
    // Prim's Algorithm: O(E log V) with priority queue
    static std::pair<int, std::vector<Edge>> prim(
            int V, 
            const std::vector<std::vector<std::pair<int, int>>>& adj,
            int start = 0,
            bool verbose = true) {
        
        std::vector<bool> inMST(V, false);
        std::vector<Edge> mstEdges;
        int totalWeight = 0;
        
        // Priority queue: (weight, vertex, parent)
        std::priority_queue<std::tuple<int, int, int>,
                          std::vector<std::tuple<int, int, int>>,
                          std::greater<>> pq;
        
        pq.push({0, start, -1});
        
        if (verbose) std::cout << "\n--- Prim's Algorithm ---\n";
        
        while (!pq.empty() && (int)mstEdges.size() < V - 1) {
            auto [w, u, parent] = pq.top();
            pq.pop();
            
            if (inMST[u]) continue;
            inMST[u] = true;
            
            if (parent != -1) {
                mstEdges.push_back({parent, u, w});
                totalWeight += w;
                if (verbose) {
                    std::cout << "  Added edge: " << parent << " -- " << u 
                              << " (weight: " << w << ")\n";
                }
            }
            
            for (auto [v, weight] : adj[u]) {
                if (!inMST[v]) {
                    pq.push({weight, v, u});
                }
            }
        }
        
        if (verbose) {
            std::cout << "Total MST weight: " << totalWeight << "\n";
        }
        
        return {totalWeight, mstEdges};
    }
};

// ============================================================================
// 2. SECOND MINIMUM SPANNING TREE (2nd MST)
// ============================================================================
// Strategy: For each edge in MST, try replacing it with a non-MST edge
// that creates a cycle, then remove the max weight edge on that cycle
// Time Complexity: O(V^2) or O(E log V) with optimization
// ============================================================================

class SecondMST {
public:
    // Find 2nd MST given original MST edges and all graph edges
    static std::pair<int, std::vector<Edge>> compute(
            int V,
            const std::vector<Edge>& allEdges,
            const std::vector<Edge>& mstEdges,
            bool verbose = true) {
        
        int mstWeight = 0;
        for (const auto& e : mstEdges) mstWeight += e.weight;
        
        // Build MST adjacency for path queries
        std::vector<std::vector<std::pair<int, int>>> mstAdj(V);
        std::set<std::pair<int, int>> mstEdgeSet;
        
        for (const auto& e : mstEdges) {
            mstAdj[e.u].push_back({e.v, e.weight});
            mstAdj[e.v].push_back({e.u, e.weight});
            mstEdgeSet.insert({std::min(e.u, e.v), std::max(e.u, e.v)});
        }
        
        int secondMSTWeight = INT_MAX;
        Edge bestSwapIn, bestSwapOut;
        
        if (verbose) {
            std::cout << "\n--- Second MST Computation ---\n";
            std::cout << "Original MST weight: " << mstWeight << "\n";
            std::cout << "Trying edge replacements:\n";
        }
        
        // For each non-MST edge, find the max edge on the MST path between its endpoints
        for (const Edge& e : allEdges) {
            auto edgeKey = std::make_pair(std::min(e.u, e.v), std::max(e.u, e.v));
            if (mstEdgeSet.count(edgeKey)) continue;  // Skip MST edges
            
            // Find max weight edge on path from e.u to e.v in MST
            auto [maxWeight, maxEdge] = findMaxOnPath(V, mstAdj, e.u, e.v);
            
            if (maxWeight == -1) continue;  // No path (shouldn't happen in connected MST)
            
            int newWeight = mstWeight - maxWeight + e.weight;
            
            if (verbose && e.weight <= maxWeight + 5) {  // Only show promising swaps
                std::cout << "  Swap out (" << maxEdge.u << "-" << maxEdge.v << ", w=" << maxWeight 
                          << ") for (" << e.u << "-" << e.v << ", w=" << e.weight << "): "
                          << "new weight = " << newWeight << "\n";
            }
            
            if (newWeight < secondMSTWeight && newWeight > mstWeight) {
                secondMSTWeight = newWeight;
                bestSwapIn = e;
                bestSwapOut = maxEdge;
            }
        }
        
        // Construct 2nd MST edges
        std::vector<Edge> secondMSTEdges;
        for (const auto& e : mstEdges) {
            if (!(e.u == bestSwapOut.u && e.v == bestSwapOut.v) &&
                !(e.u == bestSwapOut.v && e.v == bestSwapOut.u)) {
                secondMSTEdges.push_back(e);
            }
        }
        secondMSTEdges.push_back(bestSwapIn);
        
        if (verbose) {
            std::cout << "\nBest swap: Remove (" << bestSwapOut.u << "-" << bestSwapOut.v 
                      << "), Add (" << bestSwapIn.u << "-" << bestSwapIn.v << ")\n";
            std::cout << "Second MST weight: " << secondMSTWeight << "\n";
        }
        
        return {secondMSTWeight, secondMSTEdges};
    }
    
    // BFS to find maximum weight edge on path between u and v in MST
    static std::pair<int, Edge> findMaxOnPath(
            int V,
            const std::vector<std::vector<std::pair<int, int>>>& adj,
            int start, int end) {
        
        std::vector<int> parent(V, -1);
        std::vector<int> parentWeight(V, 0);
        std::queue<int> q;
        
        q.push(start);
        parent[start] = start;
        
        while (!q.empty()) {
            int u = q.front(); q.pop();
            if (u == end) break;
            
            for (auto [v, w] : adj[u]) {
                if (parent[v] == -1) {
                    parent[v] = u;
                    parentWeight[v] = w;
                    q.push(v);
                }
            }
        }
        
        if (parent[end] == -1) return {-1, {-1, -1, -1}};
        
        // Trace back and find max
        int maxWeight = 0;
        Edge maxEdge;
        int cur = end;
        while (cur != start) {
            if (parentWeight[cur] > maxWeight) {
                maxWeight = parentWeight[cur];
                maxEdge = {parent[cur], cur, parentWeight[cur]};
            }
            cur = parent[cur];
        }
        
        return {maxWeight, maxEdge};
    }
    
    // Discussion: Is it necessary to compute MST first?
    static void discussNecessity() {
        std::cout << "\n--- Is computing MST first necessary? ---\n";
        std::cout << "ANSWER: Yes, it is necessary in most efficient algorithms.\n\n";
        std::cout << "JUSTIFICATION:\n";
        std::cout << "1. The 2nd MST differs from MST by exactly one edge swap.\n";
        std::cout << "2. To find the optimal swap, we need to know which edges are in MST.\n";
        std::cout << "3. For each non-MST edge (u,v), we must find the max-weight edge\n";
        std::cout << "   on the unique path from u to v in the MST.\n\n";
        std::cout << "COUNTEREXAMPLE attempt:\n";
        std::cout << "Consider trying to find 2nd MST directly by 'almost-Kruskal':\n";
        std::cout << "  - Skip the first edge that would complete a cycle.\n";
        std::cout << "  This FAILS because the optimal swap might involve edges\n";
        std::cout << "  added earlier in Kruskal's process, not the most recent one.\n\n";
        std::cout << "COMPLEXITY: O(E log E + V^2) or O(E log E + E log V) with LCA preprocessing.\n";
    }
};

// ============================================================================
// 3. EDGE TYPES IN TRAVERSALS (DFS & BFS)
// ============================================================================

class EdgeClassification {
public:
    enum EdgeType { TREE, BACK, FORWARD, CROSS };
    
    static std::string edgeTypeName(EdgeType t) {
        switch(t) {
            case TREE: return "TREE";
            case BACK: return "BACK";
            case FORWARD: return "FORWARD";
            case CROSS: return "CROSS";
        }
        return "UNKNOWN";
    }
    
    // DFS edge classification
    static void dfsClassify(int V, const std::vector<std::vector<int>>& adj, 
                           bool directed, bool verbose = true) {
        std::vector<int> discovery(V, -1), finish(V, -1);
        std::vector<int> color(V, 0);  // 0=white, 1=gray, 2=black
        int time = 0;
        
        std::vector<std::tuple<int, int, EdgeType>> edges;
        
        std::function<void(int, int)> dfs = [&](int u, int parent) {
            color[u] = 1;
            discovery[u] = time++;
            
            for (int v : adj[u]) {
                if (!directed && v == parent) continue;  // Skip parent edge in undirected
                
                if (color[v] == 0) {  // White - unvisited
                    edges.push_back({u, v, TREE});
                    dfs(v, u);
                } else if (color[v] == 1) {  // Gray - ancestor
                    edges.push_back({u, v, BACK});
                } else {  // Black - finished
                    if (discovery[u] < discovery[v]) {
                        edges.push_back({u, v, FORWARD});
                    } else {
                        edges.push_back({u, v, CROSS});
                    }
                }
            }
            
            color[u] = 2;
            finish[u] = time++;
        };
        
        for (int i = 0; i < V; i++) {
            if (color[i] == 0) dfs(i, -1);
        }
        
        if (verbose) {
            std::cout << "\n--- DFS Edge Classification (" 
                      << (directed ? "Directed" : "Undirected") << " Graph) ---\n";
            std::cout << "Vertex:    ";
            for (int i = 0; i < V; i++) std::cout << std::setw(3) << i;
            std::cout << "\nDiscovery: ";
            for (int i = 0; i < V; i++) std::cout << std::setw(3) << discovery[i];
            std::cout << "\nFinish:    ";
            for (int i = 0; i < V; i++) std::cout << std::setw(3) << finish[i];
            std::cout << "\n\nEdge classifications:\n";
            for (auto& [u, v, type] : edges) {
                std::cout << "  " << u << " -> " << v << ": " << edgeTypeName(type) << "\n";
            }
        }
    }
    
    // BFS edge classification
    static void bfsClassify(int V, const std::vector<std::vector<int>>& adj,
                           bool directed, int start = 0, bool verbose = true) {
        std::vector<int> level(V, -1);
        std::vector<std::tuple<int, int, EdgeType>> edges;
        
        std::queue<int> q;
        q.push(start);
        level[start] = 0;
        
        std::set<std::pair<int,int>> visited;
        
        while (!q.empty()) {
            int u = q.front(); q.pop();
            
            for (int v : adj[u]) {
                auto edge = directed ? std::make_pair(u, v) 
                                    : std::make_pair(std::min(u,v), std::max(u,v));
                
                if (!directed && visited.count(edge)) continue;
                visited.insert(edge);
                
                if (level[v] == -1) {  // Unvisited
                    level[v] = level[u] + 1;
                    edges.push_back({u, v, TREE});
                    q.push(v);
                } else if (directed && level[v] <= level[u]) {
                    // In directed BFS: back edge if v is ancestor
                    if (level[v] < level[u]) {
                        edges.push_back({u, v, BACK});
                    } else {
                        edges.push_back({u, v, CROSS});
                    }
                } else {
                    edges.push_back({u, v, CROSS});
                }
            }
        }
        
        if (verbose) {
            std::cout << "\n--- BFS Edge Classification ("
                      << (directed ? "Directed" : "Undirected") << " Graph) ---\n";
            std::cout << "Vertex: ";
            for (int i = 0; i < V; i++) std::cout << std::setw(3) << i;
            std::cout << "\nLevel:  ";
            for (int i = 0; i < V; i++) std::cout << std::setw(3) << level[i];
            std::cout << "\n\nEdge classifications:\n";
            for (auto& [u, v, type] : edges) {
                std::cout << "  " << u << " -> " << v << ": " << edgeTypeName(type) << "\n";
            }
        }
    }
    
    // Prove claims about edge types
    static void proveClaims() {
        std::cout << "\n" << std::string(70, '=') << "\n";
        std::cout << "  PROOFS: EDGE TYPES IN DFS AND BFS\n";
        std::cout << std::string(70, '=') << "\n";
        
        std::cout << "\n[CLAIM 1] DFS on undirected graph: only TREE and BACK edges.\n";
        std::cout << "PROOF:\n";
        std::cout << "- When DFS visits edge (u,v) where v is already visited:\n";
        std::cout << "  * v must be gray (in current DFS path) → BACK edge\n";
        std::cout << "  * v cannot be black: if v finished before u started exploring (u,v),\n";
        std::cout << "    then v would have explored (v,u) and u would already be visited.\n";
        std::cout << "- No FORWARD: in undirected graphs, forward = back from other direction.\n";
        std::cout << "- No CROSS: cross edges require two independent branches, but in\n";
        std::cout << "  undirected graphs, any such connection makes them non-independent.\n";
        
        std::cout << "\n[CLAIM 2] BFS on undirected graph: only TREE and CROSS edges.\n";
        std::cout << "PROOF:\n";
        std::cout << "- TREE edges: first edge that discovers a vertex.\n";
        std::cout << "- Non-tree edges connect vertices at same or adjacent levels.\n";
        std::cout << "- No BACK: if (u,v) with level[v] < level[u]-1, then v would have\n";
        std::cout << "  discovered u at level[v]+1, contradiction.\n";
        std::cout << "- CROSS: edges between vertices at same level or adjacent levels.\n";
        
        std::cout << "\n[CLAIM 3] BFS on directed graph: TREE, BACK, CROSS, but no FORWARD.\n";
        std::cout << "PROOF:\n";
        std::cout << "- TREE: first edge discovering a vertex.\n";
        std::cout << "- BACK: edge (u,v) where v is ancestor of u (level[v] < level[u]).\n";
        std::cout << "- CROSS: edge (u,v) where v at same or different branch, level[v] <= level[u]+1.\n";
        std::cout << "- No FORWARD: if (u,v) is forward, level[v] > level[u]+1.\n";
        std::cout << "  But BFS processes vertices level by level, so when u is processed,\n";
        std::cout << "  any unvisited descendant v would be at level[u]+1, not further.\n";
    }
};

// ============================================================================
// 4. BICONNECTED COMPONENTS AND BRIDGES (TARJAN)
// ============================================================================

class BiconnectedComponents {
public:
    int V;
    std::vector<std::vector<int>> adj;
    std::vector<int> disc, low;
    std::vector<bool> visited;
    std::stack<std::pair<int, int>> edgeStack;
    std::vector<std::vector<std::pair<int, int>>> components;
    std::vector<std::pair<int, int>> bridges;
    std::vector<int> articulationPoints;
    int timer;
    
    BiconnectedComponents(int vertices) : V(vertices), adj(vertices), 
        disc(vertices, -1), low(vertices, -1), visited(vertices, false), timer(0) {}
    
    void addEdge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    
    void dfs(int u, int parent) {
        disc[u] = low[u] = timer++;
        visited[u] = true;
        int children = 0;
        bool isArticulation = false;
        
        for (int v : adj[u]) {
            if (!visited[v]) {
                children++;
                edgeStack.push({u, v});
                dfs(v, u);
                
                low[u] = std::min(low[u], low[v]);
                
                // Check for articulation point
                if ((parent == -1 && children > 1) || 
                    (parent != -1 && low[v] >= disc[u])) {
                    isArticulation = true;
                    
                    // Pop component
                    std::vector<std::pair<int, int>> component;
                    while (!edgeStack.empty()) {
                        auto edge = edgeStack.top();
                        edgeStack.pop();
                        component.push_back(edge);
                        if ((edge.first == u && edge.second == v) ||
                            (edge.first == v && edge.second == u)) break;
                    }
                    components.push_back(component);
                }
                
                // Check for bridge
                if (low[v] > disc[u]) {
                    bridges.push_back({u, v});
                }
            } else if (v != parent && disc[v] < disc[u]) {
                // Back edge
                edgeStack.push({u, v});
                low[u] = std::min(low[u], disc[v]);
            }
        }
        
        if (isArticulation) {
            articulationPoints.push_back(u);
        }
    }
    
    void findAll() {
        for (int i = 0; i < V; i++) {
            if (!visited[i]) {
                dfs(i, -1);
                // Remaining edges form a component
                if (!edgeStack.empty()) {
                    std::vector<std::pair<int, int>> component;
                    while (!edgeStack.empty()) {
                        component.push_back(edgeStack.top());
                        edgeStack.pop();
                    }
                    components.push_back(component);
                }
            }
        }
        
        // Remove duplicate articulation points
        std::sort(articulationPoints.begin(), articulationPoints.end());
        articulationPoints.erase(
            std::unique(articulationPoints.begin(), articulationPoints.end()),
            articulationPoints.end());
    }
    
    void printResults() {
        std::cout << "\n--- Biconnected Components (Tarjan) ---\n";
        
        std::cout << "\nDiscovery and Low arrays:\n";
        std::cout << "Vertex: ";
        for (int i = 0; i < V; i++) std::cout << std::setw(3) << i;
        std::cout << "\nDisc:   ";
        for (int i = 0; i < V; i++) std::cout << std::setw(3) << disc[i];
        std::cout << "\nLow:    ";
        for (int i = 0; i < V; i++) std::cout << std::setw(3) << low[i];
        std::cout << "\n";
        
        std::cout << "\nBiconnected Components (" << components.size() << " total):\n";
        for (size_t i = 0; i < components.size(); i++) {
            std::cout << "  Component " << i + 1 << ": ";
            std::set<int> vertices;
            for (auto [u, v] : components[i]) {
                vertices.insert(u);
                vertices.insert(v);
                std::cout << "(" << u << "-" << v << ") ";
            }
            std::cout << " | Vertices: {";
            bool first = true;
            for (int v : vertices) {
                if (!first) std::cout << ",";
                std::cout << v;
                first = false;
            }
            std::cout << "}\n";
        }
        
        std::cout << "\nBridges (" << bridges.size() << " total): ";
        for (auto [u, v] : bridges) {
            std::cout << "(" << u << "-" << v << ") ";
        }
        std::cout << "\n";
        
        std::cout << "\nArticulation Points: ";
        for (int v : articulationPoints) {
            std::cout << v << " ";
        }
        std::cout << "\n";
    }
};

// ============================================================================
// GRAPH ALGORITHMS DEMONSTRATION
// ============================================================================

void demonstrateGraphAlgorithms() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  GRAPH ALGORITHMS DEMONSTRATION\n";
    std::cout << std::string(70, '=') << "\n";
    
    // 1. Maximum Flow (Edmonds-Karp)
    std::cout << "\n### 1. MAXIMUM FLOW (Edmonds-Karp) ###\n";
    std::cout << "Graph: 6 vertices, directed edges with capacities\n";
    std::cout << "       0 --10--> 1 --10--> 3\n";
    std::cout << "       |         |         |\n";
    std::cout << "      10        2/4       10\n";
    std::cout << "       v         v         v\n";
    std::cout << "       2 --10--> 4 --10--> 5\n";
    
    EdmondsKarp ek(6);
    ek.addEdge(0, 1, 10);
    ek.addEdge(0, 2, 10);
    ek.addEdge(1, 3, 10);
    ek.addEdge(1, 4, 2);
    ek.addEdge(2, 4, 10);
    ek.addEdge(3, 5, 10);
    ek.addEdge(4, 3, 4);
    ek.addEdge(4, 5, 10);
    
    int maxFlowValue = ek.maxFlow(0, 5);
    std::cout << "\n*** MAXIMUM FLOW: " << maxFlowValue << " ***\n";
    
    std::cout << "\nComplexity Comparison:\n";
    std::cout << "  Edmonds-Karp: O(VE^2) - uses BFS for shortest augmenting paths\n";
    std::cout << "  Dinic:        O(V^2E) - uses level graphs and blocking flows\n";
    std::cout << "  Ford-Fulkerson: O(Ef) where f=max flow - can be slow with irrational capacities\n";
    
    // 2. Second MST
    std::cout << "\n\n### 2. SECOND MINIMUM SPANNING TREE ###\n";
    
    // Small example graph
    std::vector<Edge> edges = {
        {0, 1, 1}, {0, 2, 2}, {1, 2, 3}, {1, 3, 4}, {2, 3, 5}, {2, 4, 6}, {3, 4, 7}
    };
    std::vector<std::vector<std::pair<int, int>>> adj(5);
    for (auto& e : edges) {
        adj[e.u].push_back({e.v, e.weight});
        adj[e.v].push_back({e.u, e.weight});
    }
    
    auto [mstWeight, mstEdges] = MST::kruskal(5, edges);
    SecondMST::compute(5, edges, mstEdges);
    SecondMST::discussNecessity();
    
    // 3. Edge Classification
    std::cout << "\n\n### 3. EDGE TYPES IN TRAVERSALS ###\n";
    EdgeClassification::proveClaims();
    
    // Example: DFS on undirected
    std::cout << "\nExample - DFS on undirected graph:\n";
    std::cout << "Graph: 0-1, 1-2, 2-0, 2-3\n";
    std::vector<std::vector<int>> undirAdj = {{1, 2}, {0, 2}, {1, 0, 3}, {2}};
    EdgeClassification::dfsClassify(4, undirAdj, false);
    
    // Example: BFS on directed
    std::cout << "\nExample - BFS on directed graph:\n";
    std::cout << "Graph: 0->1, 0->2, 1->2, 2->0, 2->3, 3->3\n";
    std::vector<std::vector<int>> dirAdj = {{1, 2}, {2}, {0, 3}, {3}};
    EdgeClassification::bfsClassify(4, dirAdj, true, 0);
    
    // 4. Biconnected Components
    std::cout << "\n\n### 4. BICONNECTED COMPONENTS AND BRIDGES ###\n";
    
    // Graph with 5+ biconnected components including a bridge
    BiconnectedComponents bc(12);
    // Component 1: 0-1-2-0 (triangle)
    bc.addEdge(0, 1); bc.addEdge(1, 2); bc.addEdge(2, 0);
    // Bridge: 2-3
    bc.addEdge(2, 3);
    // Component 2: 3-4-5-3 (triangle)
    bc.addEdge(3, 4); bc.addEdge(4, 5); bc.addEdge(5, 3);
    // Bridge: 5-6
    bc.addEdge(5, 6);
    // Component 3: 6-7-8-9-6 (4-cycle with diagonal)
    bc.addEdge(6, 7); bc.addEdge(7, 8); bc.addEdge(8, 9); bc.addEdge(9, 6); bc.addEdge(6, 8);
    // Bridge: 9-10
    bc.addEdge(9, 10);
    // Component 4: 10-11 (bridge itself as 2-vertex component)
    bc.addEdge(10, 11);
    
    std::cout << "Graph structure:\n";
    std::cout << "  Triangle(0,1,2) --bridge(2,3)-- Triangle(3,4,5) --bridge(5,6)--\n";
    std::cout << "  4-cycle+diag(6,7,8,9) --bridge(9,10)-- Edge(10,11)\n";
    
    bc.findAll();
    bc.printResults();
    
    // 5. MST Comparison (Kruskal vs Prim)
    std::cout << "\n\n### 5. MST: KRUSKAL vs PRIM ###\n";
    std::cout << "Generating random connected graph: V=10, E=30, seed=42\n";
    
    // Use fixed seed for reproducibility
    std::mt19937 rng(42);
    std::uniform_int_distribution<int> weightDist(1, 100);
    
    int V = 10, E = 30;
    std::vector<Edge> randomEdges;
    std::vector<std::vector<std::pair<int, int>>> randomAdj(V);
    std::set<std::pair<int, int>> usedEdges;
    
    // Ensure connectivity
    for (int i = 1; i < V; i++) {
        int j = rng() % i;
        int w = weightDist(rng);
        randomEdges.push_back({j, i, w});
        randomAdj[j].push_back({i, w});
        randomAdj[i].push_back({j, w});
        usedEdges.insert({std::min(i, j), std::max(i, j)});
    }
    
    // Add remaining edges
    while ((int)randomEdges.size() < E) {
        int u = rng() % V, v = rng() % V;
        if (u != v && !usedEdges.count({std::min(u, v), std::max(u, v)})) {
            int w = weightDist(rng);
            randomEdges.push_back({u, v, w});
            randomAdj[u].push_back({v, w});
            randomAdj[v].push_back({u, w});
            usedEdges.insert({std::min(u, v), std::max(u, v)});
        }
    }
    
    std::cout << "\nEdges: ";
    for (auto& e : randomEdges) {
        std::cout << "(" << e.u << "-" << e.v << ":" << e.weight << ") ";
    }
    std::cout << "\n";
    
    auto [kruskalWeight, kruskalEdges] = MST::kruskal(V, randomEdges);
    auto [primWeight, primEdges] = MST::prim(V, randomAdj);
    
    std::cout << "\nComparison:\n";
    std::cout << "  Kruskal MST weight: " << kruskalWeight << "\n";
    std::cout << "  Prim MST weight:    " << primWeight << "\n";
    std::cout << "  Weights match: " << (kruskalWeight == primWeight ? "YES ✓" : "NO ✗") << "\n";
    
    std::cout << "\nNOTA BENE: Dijkstra's algorithm finds shortest paths from a source,\n";
    std::cout << "NOT minimum spanning trees. Prim's uses a similar greedy heap approach\n";
    std::cout << "but optimizes total tree weight, not individual path lengths.\n";
}

#endif // GRAPH_ALGORITHMS_H
