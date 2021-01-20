#ifndef DJIKSTRA_H
#define DJIKSTRA_H

#include "hexLattice.h"

#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>

typedef boost::adjacency_list< boost::listS, boost::vecS, boost::undirectedS, boost::no_property,
        boost::property< boost::edge_weight_t, int > >
        graph_t;
typedef boost::graph_traits< graph_t >::vertex_descriptor vertex_descriptor;
typedef std::pair< int, int > Edge;
typedef boost::graph_traits<graph_t>::vertex_iterator vertex_iter;
typedef boost::property_map<graph_t, boost::vertex_index_t>::type IndexMap;

// graph_t buildGraph(const vpint &edgeToVertices, color c, int L);
graph_t buildGraph(vpint edgeToVertices, color c, int L, const vpint &edgeToFaces, const sint &qubitIndices);

void shortestPaths(const graph_t& g, std::map<vertex_descriptor, std::vector<vertex_descriptor>>& excitationToPaths, std::map<vertex_descriptor, vint>& excitationToDistances, const vint& excitations);

#endif