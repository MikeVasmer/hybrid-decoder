#include "djikstra.h"
#include "hexLattice.h"

// c is the color of vertices NOT to include
graph_t buildGraph(const vpint &edgeToVertices, color c, int L)
{
    int numEdges = 0;
    int numNodes = L * L;
    for (auto p : edgeToVertices) // Count number of edges
    {
        if (p.first == 0 && p.second == 0) continue; // Edge deleted during unencoding
        if (vertexColor(p.first, L) == c || vertexColor(p.second, L) == c) continue;
        ++numEdges;
    }
    Edge edgeArray[numEdges];
    int weights[numEdges];
    int i = 0;
    for (auto p : edgeToVertices)
    {
        int v1 = p.first;
        int v2 = p.second;
        if (v1 == 0 && v2 == 0) continue; 
        if (vertexColor(v1, L) == c || vertexColor(v2, L) == c) continue;
        edgeArray[i] = Edge(v1, v2);
        weights[i] = 1;
        ++i;
    }
    int numArcs = sizeof(edgeArray) / sizeof(Edge);
    graph_t g(edgeArray, edgeArray + numArcs, weights, numNodes);
    return g;
}

void shortestPaths(const graph_t& g, std::map<vertex_descriptor, std::vector<vertex_descriptor>>& excitationToPaths, std::map<vertex_descriptor, vint>& excitationToDistances, const vint& excitations)
{
    std::vector<vertex_descriptor> p(num_vertices(g));
    std::vector<int> d(num_vertices(g));
    for (auto const excitation : excitations)
    {
        vertex_descriptor s = vertex(excitation, g);
        boost::dijkstra_shortest_paths(g, s, 
            predecessor_map(boost::make_iterator_property_map(p.begin(), get(boost::vertex_index, g)))
            .distance_map(boost::make_iterator_property_map(d.begin(), get(boost::vertex_index, g)))
        );
        excitationToPaths.insert(std::pair<vertex_descriptor, std::vector<vertex_descriptor>>(s, p));
        excitationToDistances.insert(std::pair<vertex_descriptor, std::vector<int>>(s, d));
    }
}