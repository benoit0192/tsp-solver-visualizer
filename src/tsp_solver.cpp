#include <cmath>
#include <queue>
#include <stack>
#include <chrono>
#include <vector>
#include <fstream>
#include <climits>
#include <iomanip>
#include <iostream>
#include <filesystem>
#include <unordered_map>
#include <boost/format.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/max_cardinality_matching.hpp>

#include "utils.hpp"


/** BEGIN Forward declarations **/
bool
isValidFilePath(
    const std::string& filepath
);
bool
isValidDirectoryPath(
    const std::string& dirPath);
/** END Forward declarations **/

std::vector<std::vector<float>>
getDistanceMatrix(
    std::vector<std::pair<float,float>>& nodes
)
{
    int n=nodes.size();
    std::vector<std::vector<float>> distMat(n, std::vector<float>(n));
    for (int i=0; i<n; ++i) {
        auto& [x1,y1] = nodes[i];
        for (int j=0; j<n; ++j) {
            auto& [x2,y2] = nodes[j];
            distMat[i][j] = std::sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
        }
    }
    return distMat;
}

void
printDistMat(
    std::vector<std::vector<float>>& distMat,
    int precision=2,
    int width=10
)
{
    int n=distMat.size();
    for (int i=0; i<n; ++i) {
        for (int j=0; j<n; ++j) {
            std::cout << std::fixed
                      << std::setprecision(precision)
                      << std::setw(width)
                      << distMat[i][j];
        }
        std::cout << '\n';
    }
}

std::vector<std::pair<int,int>>
findMSTVertices(
    std::vector<std::vector<float>>& distMat
)
{
    /**
     * Find the Minimum Spawning Tree based on Prim's algorithm
     */
    int n=distMat.size();
    if(n==0) return {};
    std::vector<int> vis_nodes(n,0);
    vis_nodes[0]=1;
    auto comp = [](const auto& p1, const auto& p2){
        return p1.first > p2.first;
    };
    std::priority_queue<std::pair<float,std::pair<int,int>>,
                        std::vector<std::pair<float,std::pair<int,int>>>,
                        decltype(comp)> pq;
    /* Init the queue */
    for (int v=1; v<n; ++v) {
        if (distMat[0][v]>0)
            pq.push({distMat[0][v],{0,v}});
    }

    std::vector<std::pair<int,int>> edges;
    float totWeight = 0;
    while (!pq.empty()) {
        auto  [weight,e] = pq.top(); pq.pop();
        auto& [u,v] = e;
        if (vis_nodes[v])   continue;
        vis_nodes[v]=1;
        edges.push_back({u,v});
        totWeight += distMat[u][v];
        for (int w=0; w<n; ++w) {
            if (!vis_nodes[w] && distMat[v][w]>0)
                pq.push({distMat[v][w],{v,w}});
        }
    }
    if ((int)edges.size() != (n-1)) {
        throw std::runtime_error(
            (boost::format("Unexpected number of edges: expected %1%, got %2%.")
                          % (n-1) % edges.size()).str());
    }
    return edges;
}

void
printMSTEdges(
    std::vector<std::pair<int,int>>& edges
)
{
    std::string str_edges   = "  ";
    std::string str_sepator = " - ";
    for (auto& e: edges) {
         str_edges +=  "(" + std::to_string(e.first) + ","
                           + std::to_string(e.second) + ")" + str_sepator;
    }
    auto i = str_edges.rfind(str_sepator);
    if (i!=str_edges.npos)   str_edges = str_edges.substr(0,i);
    std::cout << str_edges << '\n';
}

std::vector<int>
findOddDegreeVertices(
    std::vector<std::pair<int,int>>& edges
)
{
    std::unordered_map<int,int> degMap;
    for (auto& [u,v]: edges) {
        degMap[u]++; degMap[v]++;
    }
    std::vector<int> vertices;
    for (auto& [u,deg]: degMap) {
        if (deg&1)   vertices.push_back(u);
    }
    if (vertices.size()%2 != 0) {
        throw std::runtime_error(
            (boost::format("Expected number of vertices to be even: got %1%.")
                          % vertices.size()).str());
    }
    return vertices;
}

void
printVector(
    std::vector<int>& vertices
)
{
    std::string str_vertices  = "  ";
    std::string str_separator = ", ";
    for (int v: vertices) {
        str_vertices += std::to_string(v) + str_separator;
    }
    auto i=str_vertices.rfind(str_separator);
    if (i != str_vertices.npos) {
        str_vertices = str_vertices.substr(0,i);
    }
    std::cout << str_vertices << '\n';
}

float
getDijkstraDistance(
    int start,
    int goal,
    std::vector<std::vector<float>>& weightsMat
)
{
    int n=weightsMat.size();
    if (n==0) {
        throw std::runtime_error("Weights matrix expected to be non null.");
    }

    auto comp = [](const auto& p1, const auto& p2){
        return p1.first > p2.first;
    };
    std::priority_queue<std::pair<float,int>,
                        std::vector<std::pair<float,int>>,
                        decltype(comp)> pq;
    pq.push({0.0,start});

    std::vector<float> nodes_weight(n,std::numeric_limits<float>::max());
    nodes_weight[start]=0.0;

    while (!pq.empty()) {
        auto [w,u] = pq.top(); pq.pop();
        if (u==goal) return w;
        for (int v=0; v<(int)weightsMat[u].size(); ++v) {
            if (u!=v &&
                weightsMat[u][v]>0 &&
                w+weightsMat[u][v]<nodes_weight[v])
            {
                nodes_weight[v] = w + weightsMat[u][v];
                pq.push({nodes_weight[v], v});
            }
        }
    }
    return std::numeric_limits<float>::max();
}

struct pair_hash {
    template <class T1, class T2>
    size_t operator()(const std::pair<T1, T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ h2;
    }
};

std::unordered_map<std::pair<int,int>,float, pair_hash>
findShortestPairDistances(
    std::vector<int>& vertices,
    std::vector<std::vector<float>> weightsMat
)
{
    int n=vertices.size();
    std::unordered_map<std::pair<int, int>, float, pair_hash> distMap;
    for (int i=0; i<n-1; ++i) {
        for (int j=i+1; j<n; ++j) {
            float d = getDijkstraDistance(vertices[i], vertices[j], weightsMat);
            distMap[{vertices[i],vertices[j]}] = d;
        }
    }
    return distMap;
}

void
printDistMap(
    std::unordered_map<std::pair<int,int>,float, pair_hash> distMap
)
{
    for (auto& [p, w]: distMap) {
        auto& [u,v] = p;
        std::cout << "  ";
        std::cout << "(" + std::to_string(u) + "," +
                           std::to_string(v) + "): " + std::to_string(w) + '\n';
    }
}

std::unordered_map<int,int>
findPairMatching(
    std::unordered_map<std::pair<int,int>,float, pair_hash>& distMap
)
{
    typedef boost::adjacency_list<boost::vecS,
                                  boost::vecS,
                                  boost::undirectedS,
                                  boost::no_property,
                                  boost::property<boost::edge_weight_t,
                                  float>> Graph;
    typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;

    Graph g;
    for (auto& [p, w] : distMap) {
        auto& [u,v] = p;
        boost::add_edge(u, v, w, g);
    }
    std::vector<Vertex> mate(boost::num_vertices(g),
                             boost::graph_traits<Graph>::null_vertex());

    bool success = boost::checked_edmonds_maximum_cardinality_matching(
                                                                    g,&mate[0]);
    if (!success) {
        std::cerr << "ERROR computing maximum cardinality matching." << '\n';
    }
    std::unordered_map<int,int> matching;
    for (std::size_t u = 0; u < mate.size(); ++u) {
        Vertex v = mate[u];
        if (v != boost::graph_traits<Graph>::null_vertex() && u < v) {
            matching[u] = v;
        }
    }
    return matching;
}

void
printMatchingPair(
    std::unordered_map<int,int>& matching
)
{
    for (auto& [u,v] : matching) {
        std::cout << "  ";
        std::cout << "(" << u << "," << v << ")" << '\n';
    }
}

std::vector<std::pair<int,int>>
mergeEdges(
    std::vector<std::pair<int,int>>& mstEdges,
    std::unordered_map<int,int>& perfectMatching
)
{
    auto makeUniform = [](const int u, const int v){
        return (u < v) ? std::make_pair(u,v) :
                         std::make_pair(v,u);
    };
    std::set<std::pair<int,int>> edgeSet;
    for (auto& e: mstEdges)
        edgeSet.insert(makeUniform(e.first,e.second));
    for (auto& [u,v]: perfectMatching)
        edgeSet.insert(makeUniform(u,v));
    return std::vector<std::pair<int,int>>(edgeSet.begin(), edgeSet.end());
}

std::vector<int>
findEulerianCircuit(
    std::vector<std::pair<int,int>>& edges
)
{
    /**
     * Find the Eulerian circuit based on Hierholzer's algorithm
     */
    if (edges.empty())   return {};

    std::unordered_map<int, std::vector<int>> adjList;
    for (auto& [u,v]: edges) {
        adjList[u].push_back(v);
        adjList[v].push_back(u);
    }

    std::vector<int> circuit;
    std::stack<int> stack;
    int start = edges[0].first;
    stack.push(start);

    while (!stack.empty()) {
        int u = stack.top();

        if (!adjList[u].empty()) {
            int v = adjList[u].back();
            adjList[u].pop_back();

            auto& neighbors = adjList[v];
            neighbors.erase(std::remove(neighbors.begin(), neighbors.end(), u),
                            neighbors.end());

            stack.push(v);
        }
        else {
            circuit.push_back(u);
            stack.pop();
        }
    }
    return circuit;
}

std::vector<int>
convertEulerianToHamiltonian(
    std::vector<int>& circuit
)
{
    std::unordered_set<int> vis;
    std::vector<int> hCircuit;
    for (int& x: circuit) {
        if (vis.contains(x))
            continue;
        vis.insert(x);
        hCircuit.push_back(x);
    }
    return hCircuit;
}

std::vector<int>
resolveTSP(
    std::vector<std::pair<float,float>>& nodes,
    bool verbose=false
)
{
    std::vector<std::vector<float>> distMat = getDistanceMatrix(nodes);

    std::vector<std::pair<int,int>> mstEdges = findMSTVertices(distMat);

    std::vector<int> oddDegVertices = findOddDegreeVertices(mstEdges);

    using um_t = std::unordered_map<std::pair<int,int>,float, pair_hash>;
    um_t oddDegDistMap = findShortestPairDistances(oddDegVertices, distMat);

    std::unordered_map<int,int> matching = findPairMatching(oddDegDistMap);

    std::vector<std::pair<int,int>> edges = mergeEdges(mstEdges, matching);

    std::vector<int> circuit = findEulerianCircuit(edges);

    std::vector<int> hCircuit = convertEulerianToHamiltonian(circuit);

    if (verbose) {
        std::cout << "Distance matrix:" << '\n';
        printDistMat(distMat);
        std::cout << "MST edges:" << '\n';
        printMSTEdges(mstEdges);
        std::cout << "Odd-degree vertices:" << '\n';
        printVector(oddDegVertices);
        std::cout << "Odd-degree vertices distance map:" << '\n';
        printDistMap(oddDegDistMap);
        std::cout << "Pair matching:" << '\n';
        printMatchingPair(matching);
        std::cout << "Merged edges:" << '\n';
        printMSTEdges(edges);
        std::cout << "Eulerian circuit:" << '\n';
        printVector(circuit);
        std::cout << "Hamiltonian circuit:" << '\n';
        printVector(hCircuit);
    }
    return hCircuit;
}

std::vector<std::pair<float, float>>
readInputFromFile(
    const std::string& filename
)
{
    std::vector<std::pair<float, float>> nodes;
    std::ifstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error(
                (boost::format("Failed to open file: %1%.") % filename).str());
    }
    std::string line;
    while (std::getline(file, line)) {
        if (line.contains(',')) {
            int cIx = line.find(',');
            std::string xStr = line.substr(0,cIx);
            std::string yStr = line.substr(cIx+1);
            if (xStr.empty() || yStr.empty()) {
                throw std::runtime_error(
                    (boost::format("Failed to parse entry: %1%.") % line).str());
            }
            float x = std::stof(xStr);
            float y = std::stof(yStr);
            nodes.emplace_back(x, y);
        }
    }
    return nodes;
}

void
writeOutputToFile(
    const std::string& filename,
    std::vector<std::pair<float,float>>& nodes,
    std::vector<int>& path
)
{
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Failed to open file for writing." << '\n';
        return;
    }
    /* Dump nodes */
    for (auto& [x,y] : nodes) {
        outFile << x << ':' << y << '\n';
    }
    /* Dump path */
    std::string pathStr;
    for (int& ix: path) {
         pathStr += std::to_string(ix) + '-';
    }
    if (!pathStr.empty()) {
        pathStr.pop_back();
        pathStr += '\n';
    }
    outFile << pathStr;

    outFile.close();
}

bool
parseArgs(
    int argc,
    char *argv[],
    std::string& nodesInPath,
    std::string& outDir,
    bool& verbose
)
{
    std::string usage = std::string("Usage:\n")
         + "./vanilla-tsp <DATA_FILEPATH> <OUT_DIR>\n"
         + "  <DATA_FILEPATH> is the path to a file containing user-defined data\n"
         + "  in the form of node locations. One node definition per line in the format:\n"
         + "  'x,y', where x and y are floating-point numbers.\n"
         + "  <OUT_DIR> is an optional argument specifying the output directory\n"
         + "  where the results will be saved.\n";

    if (argc < 2) {
        std::cerr << "ERROR: Unexpected number of input arguments.\n";
        std::cout << usage << '\n';
        return false;
    }

    nodesInPath = argv[1];

    for (int i=2; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--verbose") {
            verbose = true;
        }
        else if (arg.substr(0, 2) == "--") {
            std::cerr << "ERROR: Unknown argument " << arg << "\n";
            return false;
        }
        else {
            outDir = arg;
        }
    }
    if (!isValidFilePath(nodesInPath)) {
        std::cerr << (boost::format(
                      "ERROR: Provided data filepath does not exist:\n"
                      "        %1%\n")
                      % nodesInPath).str();
        return false;
    }
    if (!isValidDirectoryPath(outDir)) {
        std::cerr << (boost::format(
                      "ERROR: Provided output directory does not exist:\n"
                      "        %1%\n")
                      % outDir).str();
        return false;
    }

    return true;
}

int main(int argc, char *argv[]) {
    auto start = std::chrono::high_resolution_clock::now();

    std::cout << "== Travel Salesman Problem ==\n"
              << "== The following problem resolution is based on \n"
              << "== Christofides' algorithm.\n";

    std::string nodesInPath   = "nodes_data";
    std::string outDir        = ".";
    bool verbose = false;

    if (!parseArgs(argc, argv, nodesInPath, outDir, verbose))
        return -1;

    std::string solverOutPath = (fs::path(outDir) / fs::path("solver_data")).string();

    auto nodes = readInputFromFile(nodesInPath);
    std::vector<int> path = resolveTSP(nodes, verbose);
    writeOutputToFile(solverOutPath, nodes, path);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> d = end - start;
    std::cout << "Program execution time: " << d.count() << " ms" << '\n';

    return 0;
}
