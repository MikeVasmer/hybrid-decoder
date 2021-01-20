#include "gtest/gtest.h"
#include "djikstra.h"
#include "hexLattice.h"

TEST(buildGraph, no_unencoding)
{
    vint Ls = {6, 9, 12, 15};
    std::vector<color> colors = {r, g, b};
    for (auto L : Ls)
    {
        for (auto c : colors)
        {
            auto edgeToVertices = buildEdgeToVertices(L);
            auto edgeToFaces = buildEdgeToFaces(L);
            sint qubitIndices;
            auto g = buildGraph(edgeToVertices, c, L, edgeToFaces, qubitIndices);
            EXPECT_EQ(num_vertices(g), L * L);
            EXPECT_EQ(num_edges(g), L * L);
            for (int v = 0; v < L * L; ++v)
            {
                if (vertexColor(v, L) == c)
                {
                    EXPECT_FALSE(edge(v, neigh(v, x, 1, L), g).second);
                    EXPECT_FALSE(edge(v, neigh(v, y, 1, L), g).second);
                    EXPECT_FALSE(edge(v, neigh(v, xy, 1, L), g).second);
                    EXPECT_FALSE(edge(v, neigh(v, x, -1, L), g).second);
                    EXPECT_FALSE(edge(v, neigh(v, y, -1, L), g).second);
                    EXPECT_FALSE(edge(v, neigh(v, xy, -1, L), g).second);
                }
                else
                {
                    vint dirs = {x, y, xy};
                    vint signs = {1, -1};
                    for (auto dir : dirs)
                    {
                        for (auto sign : signs)
                        {
                            if (vertexColor(neigh(v, dir, sign, L), L) == c)
                            {
                                EXPECT_FALSE(edge(v, neigh(v, dir, sign, L), g).second);
                            }
                            else
                            {
                                EXPECT_TRUE(edge(v, neigh(v, dir, sign, L), g).second);
                            }
                        }
                    }
                }
            }
            // IndexMap index = boost::get(boost::vertex_index, g);
            // std::cout << "edges(g) = ";
            // boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
            // for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
            // {
            //     std::cout << "(" << index[source(*ei, g)] 
            //               << "," << index[target(*ei, g)] << ") ";
            // }
            // std::cout << std::endl;
        }
    }
}

TEST(buildGraph, full_unencoding)
{
    vint Ls = {6, 9, 12, 15};
    std::vector<color> colors = {b, g};
    std::random_device rd{}; 
    std::mt19937 engine{rd()};
    std::uniform_real_distribution<double> dist{0.0, 1.0};
    for (auto L : Ls)
    {
        auto vertexToQubits = buildVertexToFaces(L);
        auto edgeToVertices = buildEdgeToVertices(L);
        auto vertexToEdges = buildVertexToEdges(L);
        auto logicals = buildLogicals(L);
        sint unencodedVertices, qubitIndices;
        auto edgeToFaces = buildEdgeToFaces(L);
        vint qubits;
        double p = 1.0;
        auto faceToEdges = buildFaceToEdges(L);
        auto lift = buildLift(L, vertexToQubits, vertexToEdges, faceToEdges);
        unencode(vertexToQubits, edgeToVertices, unencodedVertices, qubitIndices, qubits, logicals, edgeToFaces, vertexToEdges, L, p, engine, dist, false, true, lift);
        for (auto c : colors)
        {
            auto g = buildGraph(edgeToVertices, c, L, edgeToFaces, qubitIndices);
            EXPECT_EQ(num_vertices(g), L * L);
            EXPECT_EQ(num_edges(g), 2 * L * L / 3); // Half of new edges
            for (int v = 0; v < L * L; ++v) // Check no CC edges present
            {
                vint dirs = {x, y, xy};
                vint signs = {1, -1};
                for (auto dir : dirs)
                {
                    for (auto sign : signs)
                    {
                        EXPECT_FALSE(edge(v, neigh(v, dir, sign, L), g).second);
                    }
                }
            }
            for (int v = 0; v < L * L; ++v)
            {
                if (vertexColor(v, L) == c || vertexColor(v, L) == r) continue;
                int neighbor = neigh(neigh(v, xy, 1, L), y, 1, L);
                EXPECT_TRUE(edge(v, neighbor, g).second);
                neighbor = neigh(neigh(v, x, 1, L), xy, 1, L);
                EXPECT_TRUE(edge(v, neighbor, g).second);
            }
        }
    }
}

TEST(shortestPaths, no_unencoding)
{
    int L = 6;
    color c = g; // rb restricted lattice
    auto edgeToVertices = buildEdgeToVertices(L);
    auto edgeToFaces = buildEdgeToFaces(L);
    sint qubitIndices;
    auto g = buildGraph(edgeToVertices, c, L, edgeToFaces, qubitIndices);
    std::map<vertex_descriptor, std::vector<vertex_descriptor>> excitationToPaths;
    std::map<vertex_descriptor, vint> excitationToDistances;
    vint excitations = {0};
    shortestPaths(g, excitationToPaths, excitationToDistances, excitations);
    EXPECT_EQ(excitationToDistances[0][0], 0);
    EXPECT_EQ(excitationToDistances[0][1], INT_MAX);  
    EXPECT_EQ(excitationToDistances[0][2], 3);
    EXPECT_EQ(excitationToDistances[0][3], 4);
    EXPECT_EQ(excitationToDistances[0][4], INT_MAX); 
    EXPECT_EQ(excitationToDistances[0][5], 1);
    EXPECT_EQ(excitationToDistances[0][6], INT_MAX); 
    EXPECT_EQ(excitationToDistances[0][7], 1);
    EXPECT_EQ(excitationToDistances[0][8], 2);
    EXPECT_EQ(excitationToDistances[0][9], INT_MAX); 
    EXPECT_EQ(excitationToDistances[0][10], 3);
    EXPECT_EQ(excitationToDistances[0][11], 2);
    EXPECT_EQ(excitationToDistances[0][12], 3);
    EXPECT_EQ(excitationToDistances[0][13], 2);
    EXPECT_EQ(excitationToDistances[0][14], INT_MAX); 
    EXPECT_EQ(excitationToDistances[0][15], 3);
    EXPECT_EQ(excitationToDistances[0][16], 4);
    EXPECT_EQ(excitationToDistances[0][17], INT_MAX);
    EXPECT_EQ(excitationToDistances[0][18], 4);
    EXPECT_EQ(excitationToDistances[0][19], INT_MAX);
    EXPECT_EQ(excitationToDistances[0][20], 3);
    EXPECT_EQ(excitationToDistances[0][21], 4);
    EXPECT_EQ(excitationToDistances[0][22], INT_MAX);
    EXPECT_EQ(excitationToDistances[0][23], 3);
    EXPECT_EQ(excitationToDistances[0][24], INT_MAX);
    EXPECT_EQ(excitationToDistances[0][25], 3);
    EXPECT_EQ(excitationToDistances[0][26], 4);
    EXPECT_EQ(excitationToDistances[0][27], INT_MAX);
    EXPECT_EQ(excitationToDistances[0][28], 3);
    EXPECT_EQ(excitationToDistances[0][29], 2);
    EXPECT_EQ(excitationToDistances[0][30], 1);
    EXPECT_EQ(excitationToDistances[0][31], 2);
    EXPECT_EQ(excitationToDistances[0][32], INT_MAX);
    EXPECT_EQ(excitationToDistances[0][33], 3);
    EXPECT_EQ(excitationToDistances[0][34], 2);
    EXPECT_EQ(excitationToDistances[0][35], INT_MAX);
}

TEST(shortestPaths, full_unencoding)
{
    int L = 6;
    color c = b; // rg restricted lattice
    auto edgeToVertices = buildEdgeToVertices(L);
    auto vertexToEdges = buildVertexToEdges(L);
    std::random_device rd{}; 
    std::mt19937 engine{rd()};
    std::uniform_real_distribution<double> dist{0.0, 1.0};
    auto vertexToQubits = buildVertexToFaces(L);
    auto logicals = buildLogicals(L);
    sint unencodedVertices, qubitIndices;
    auto edgeToFaces = buildEdgeToFaces(L);
    vint qubits;
    double p = 1.0;
    auto faceToEdges = buildFaceToEdges(L);
    auto lift = buildLift(L, vertexToQubits, vertexToEdges, faceToEdges);
    unencode(vertexToQubits, edgeToVertices, unencodedVertices, qubitIndices, qubits, logicals, edgeToFaces, vertexToEdges, L, p, engine, dist, false, true, lift);
    auto g = buildGraph(edgeToVertices, c, L, edgeToFaces, qubitIndices);
    std::map<vertex_descriptor, std::vector<vertex_descriptor>> excitationToPaths;
    std::map<vertex_descriptor, vint> excitationToDistances;
    vint excitations = {14};
    shortestPaths(g, excitationToPaths, excitationToDistances, excitations);
    EXPECT_EQ(excitationToDistances[14][0], INT_MAX);
    EXPECT_EQ(excitationToDistances[14][1], 1);
    EXPECT_EQ(excitationToDistances[14][2], INT_MAX);
    EXPECT_EQ(excitationToDistances[14][3], INT_MAX);
    EXPECT_EQ(excitationToDistances[14][4], 2);
    EXPECT_EQ(excitationToDistances[14][5], INT_MAX);
    EXPECT_EQ(excitationToDistances[14][6], 1);
    EXPECT_EQ(excitationToDistances[14][7], INT_MAX);
    EXPECT_EQ(excitationToDistances[14][8], INT_MAX);
    EXPECT_EQ(excitationToDistances[14][9], 2);
    EXPECT_EQ(excitationToDistances[14][10], INT_MAX);
    EXPECT_EQ(excitationToDistances[14][11], INT_MAX);
    EXPECT_EQ(excitationToDistances[14][12], INT_MAX);
    EXPECT_EQ(excitationToDistances[14][13], INT_MAX);
    EXPECT_EQ(excitationToDistances[14][14], 0);
    EXPECT_EQ(excitationToDistances[14][15], INT_MAX);
    EXPECT_EQ(excitationToDistances[14][16], INT_MAX);
    EXPECT_EQ(excitationToDistances[14][17], 3);
    EXPECT_EQ(excitationToDistances[14][18], INT_MAX);
    EXPECT_EQ(excitationToDistances[14][19], 2);
    EXPECT_EQ(excitationToDistances[14][20], INT_MAX);
    EXPECT_EQ(excitationToDistances[14][21], INT_MAX);
    EXPECT_EQ(excitationToDistances[14][22], 1);
    EXPECT_EQ(excitationToDistances[14][23], INT_MAX);
    EXPECT_EQ(excitationToDistances[14][24], 2);
    EXPECT_EQ(excitationToDistances[14][25], INT_MAX);
    EXPECT_EQ(excitationToDistances[14][26], INT_MAX);
    EXPECT_EQ(excitationToDistances[14][27], 1);
    EXPECT_EQ(excitationToDistances[14][28], INT_MAX);
    EXPECT_EQ(excitationToDistances[14][29], INT_MAX);
    EXPECT_EQ(excitationToDistances[14][30], INT_MAX);
    EXPECT_EQ(excitationToDistances[14][31], INT_MAX);
    EXPECT_EQ(excitationToDistances[14][32], 3);
    EXPECT_EQ(excitationToDistances[14][33], INT_MAX);
    EXPECT_EQ(excitationToDistances[14][34], INT_MAX);
    EXPECT_EQ(excitationToDistances[14][35], 2);
}