#include "gtest/gtest.h"
#include "hexLattice.h"
#include <algorithm>

TEST(coordToIndex, correct_output)
{
    int L = 9;
    coord c = {0, 0};
    EXPECT_EQ(coordToIndex(c, L), 0);
    c = {4, 5};
    EXPECT_EQ(coordToIndex(c, L), 49);
}

TEST(indexToCoord, correct_output)
{
    int L = 12;
    int v = 0;
    coord expected = {0, 0};
    EXPECT_EQ(indexToCoord(v, L), expected);
    v = 69;
    expected = {9, 5};
    EXPECT_EQ(indexToCoord(v, L), expected);
}

TEST(neigh, correct_output)
{
    int L = 6;
    int v = 0;
    EXPECT_EQ(neigh(v, x, 1, L), 1);
    EXPECT_EQ(neigh(v, y, 1, L), 6);
    EXPECT_EQ(neigh(v, xy, 1, L), 7);
    EXPECT_EQ(neigh(v, x, -1, L), 5);
    EXPECT_EQ(neigh(v, y, -1, L), 30);
    EXPECT_EQ(neigh(v, xy, -1, L), 35);
}

TEST(edgeIndex, correct_output)
{
    int L = 3;
    int v = 8;
    EXPECT_EQ(edgeIndex(v, x, 1, L), 24);
    EXPECT_EQ(edgeIndex(v, y, 1, L), 25);
    EXPECT_EQ(edgeIndex(v, xy, 1, L), 26);
    EXPECT_EQ(edgeIndex(v, x, -1, L), 21);
    EXPECT_EQ(edgeIndex(v, y, -1, L), 16);
    EXPECT_EQ(edgeIndex(v, xy, -1, L), 14);
}

TEST(faceIndex, correct_output)
{
    int v = 77;
    EXPECT_EQ(faceIndex(v, x), 154);
    EXPECT_EQ(faceIndex(v, y), 155);
}

TEST(vertexColor, correct_output)
{
    int L = 15;
    EXPECT_EQ(vertexColor(221, L), g);
    EXPECT_EQ(vertexColor(17, L), r);
    EXPECT_EQ(vertexColor(70, L), b);
}

TEST(buildFaceToEdges, correct_L12)
{
    int L = 12;
    auto faceToEdges = buildFaceToEdges(L);
    EXPECT_EQ(faceToEdges[69][0], 103);
    EXPECT_EQ(faceToEdges[69][1], 104);
    EXPECT_EQ(faceToEdges[69][2], 138);
    EXPECT_EQ(faceToEdges[42][0], 63);
    EXPECT_EQ(faceToEdges[42][1], 65);
    EXPECT_EQ(faceToEdges[42][2], 67);
}

TEST(buildFaceToEdges, correct_size)
{
    vint Ls = {3, 6, 9, 12, 15};
    for (auto const L : Ls)
    {
        auto faceToEdges = buildFaceToEdges(L);
        for (auto const &e : faceToEdges)
        {
            EXPECT_EQ(e.size(), 3);
        }
        EXPECT_EQ(faceToEdges.size(), 2 * L * L);
    }
}

TEST(buildVertexToEdges, correct_L6)
{
    int L = 6;
    auto vertexToEdges = buildVertexToEdges(L);
    EXPECT_EQ(vertexToEdges[33][0], 80);
    EXPECT_EQ(vertexToEdges[33][1], 82);
    EXPECT_EQ(vertexToEdges[33][2], 96);
    EXPECT_EQ(vertexToEdges[33][3], 99);
    EXPECT_EQ(vertexToEdges[33][4], 100);
    EXPECT_EQ(vertexToEdges[33][5], 101);
}

TEST(buildVertexToEdges, correct_size)
{
    vint Ls = {3, 6, 9, 12, 15};
    for (auto const L : Ls)
    {
        auto vertexToEdges = buildVertexToEdges(L);
        for (auto const &e : vertexToEdges)
        {
            EXPECT_EQ(e.size(), 6);
        }
        EXPECT_EQ(vertexToEdges.size(), L * L);
    }
}

TEST(buildVertexToFaces, correct_L3)
{
    int L = 3;
    auto vertexToFaces = buildVertexToFaces(L);
    vint expectedFaces = {0, 1, 4, 13, 16, 17};
    int i = 0;
    for (auto f : vertexToFaces[0])
    {
        EXPECT_EQ(f, expectedFaces[i]);
        ++i;
    }
}

TEST(buildVertexToFaces, correct_sizes)
{
    vint Ls = {3, 6, 9, 12, 15};
    for (auto const L : Ls)
    {
        auto vertexToFaces = buildVertexToFaces(L);
        for (auto const &f : vertexToFaces)
        {
            EXPECT_EQ(f.size(), 6);
        }
        EXPECT_EQ(vertexToFaces.size(), L * L);
    }
}

TEST(buildEdgeToVertices, correct_L6)
{
    int L = 6;
    auto edgeToVertices = buildEdgeToVertices(L);
    std::pair<int, int> expected{18, 23};
    EXPECT_EQ(edgeToVertices[69], expected);
    expected = {14, 20};
    EXPECT_EQ(edgeToVertices[43], expected);
    expected = {9, 16};
    EXPECT_EQ(edgeToVertices[29], expected);
}

TEST(buildEdgeToVertices, correct_size)
{
    vint Ls = {3, 6, 9, 12, 15};
    for (auto const L : Ls)
    {
        auto edgeToVertices = buildEdgeToVertices(L);
        EXPECT_EQ(edgeToVertices.size(), 3 * L * L);
    }
}

TEST(buildLogicals, correct_L3)
{
    int L = 3;
    auto logicals = buildLogicals(L);
    sint expectedL1{0, 1, 10, 11, 14, 15};
    sint expectedL2{2, 5, 6, 9, 13, 16};
    sint expectedL3{2, 3, 6, 7, 16, 17};
    sint expectedL4{1, 4, 8, 11, 12, 15};
    EXPECT_EQ(logicals[0], expectedL1);
    EXPECT_EQ(logicals[1], expectedL2);
    EXPECT_EQ(logicals[2], expectedL3);
    EXPECT_EQ(logicals[3], expectedL4);
}

TEST(buildLogicals, correct_size)
{
    vint Ls = {3, 6, 9, 12, 15};
    for (auto const L : Ls)
    {
        auto logicals = buildLogicals(L);
        for (auto const &l : logicals)
        {
            EXPECT_EQ(l.size(), 2 * L);
        }
    }
}

TEST(modifyVertexToQubits, correct_output_L6)
{
    int L = 6;
    auto vertexToQubits = buildVertexToFaces(L);
    int v = 0;
    int e = 3 * L * L;
    vvint vertexToEdgesL;
    vint localVertexMap;
    vpint edgeToVerticesL;
    std::map<std::pair<int, int>, std::vector<int>> logicalOperatorMapX;
    std::map<std::pair<int, int>, std::vector<int>> logicalOperatorMapZ;
    std::mt19937 engine;
    localUnencoding(vertexToEdgesL, edgeToVerticesL, localVertexMap, logicalOperatorMapX, logicalOperatorMapZ, v, L, engine, false);
    modifyVertexToQubits(vertexToQubits, vertexToEdgesL, localVertexMap, v, e, L);
    sint expectedQ = {e, e + 1, e + 2, e + 3};
    EXPECT_EQ(vertexToQubits[v], expectedQ);
    expectedQ = {2, 3, 60, 63, e};
    EXPECT_EQ(vertexToQubits[neigh(v, x, 1, L)], expectedQ);
    expectedQ  = {11, 12, 13, 22, e + 1};
    EXPECT_EQ(vertexToQubits[neigh(v, y, 1, L)], expectedQ);
    expectedQ = {3, 12, 14, 15, e + 2, e + 3};
    EXPECT_EQ(vertexToQubits[neigh(v, xy, 1, L)], expectedQ);
    expectedQ = {8, 11, 68, 69, e + 3};
    EXPECT_EQ(vertexToQubits[neigh(v, x, -1, L)], expectedQ);
    expectedQ = {49, 58, 59, 60, e + 2};
    EXPECT_EQ(vertexToQubits[neigh(v, y, -1, L)], expectedQ);
    expectedQ = {56, 57, 59, 68, e, e + 1};
    EXPECT_EQ(vertexToQubits[neigh(v, xy, -1, L)], expectedQ);
}

TEST(modifyVertexToQubits, full_unencode)
{
    // All vertices should be 4-valent after full unencode (gives 2x toric codes)
    vint Ls = {3, 6, 9, 12, 15};
    for (auto const L : Ls)
    {
        auto vertexToQubits = buildVertexToFaces(L);
        int e = 3 * L * L;
        for (int v = 0; v < L * L; ++v)
        {
            if (vertexColor(v, L) != r) continue;
            vvint vertexToEdgesL;
            vint localVertexMap;
            vpint edgeToVerticesL;
            std::map<std::pair<int, int>, std::vector<int>> logicalOperatorMapX;
            std::map<std::pair<int, int>, std::vector<int>> logicalOperatorMapZ;
            std::mt19937 engine;
            localUnencoding(vertexToEdgesL, edgeToVerticesL, localVertexMap, logicalOperatorMapX, logicalOperatorMapZ, v, L, engine, false);
            modifyVertexToQubits(vertexToQubits, vertexToEdgesL, localVertexMap, v, e, L);
            e += 4;
        }
        for (auto const qs : vertexToQubits)
        {
            EXPECT_EQ(qs.size(), 4);
        }
    }
}

TEST(modifyEdgeToVertices, correct_L6)
{
    int L = 6;
    auto edgeToVertices = buildEdgeToVertices(L);
    // Unencode (0,0)
    int v = 0;
    int e = 3 * L * L;
    vvint vertexToEdgesL;
    vint localVertexMap;
    vpint edgeToVerticesL;
    std::map<std::pair<int, int>, std::vector<int>> logicalOperatorMapX;
    std::map<std::pair<int, int>, std::vector<int>> logicalOperatorMapZ;
    std::mt19937 engine;
    localUnencoding(vertexToEdgesL, edgeToVerticesL, localVertexMap, logicalOperatorMapX, logicalOperatorMapZ, v, L, engine, false);
    modifyEdgeToVertices(edgeToVertices, edgeToVerticesL, v, L, localVertexMap);
    std::pair<int, int> expected{1, 35};
    EXPECT_EQ(edgeToVertices[e], expected);
    expected = {6, 35};
    EXPECT_EQ(edgeToVertices[e + 1], expected);
    expected = {7, 30};
    EXPECT_EQ(edgeToVertices[e + 2], expected);
    expected = {5, 7};
    EXPECT_EQ(edgeToVertices[e + 3], expected);
    expected = {0, 0};
    EXPECT_EQ(edgeToVertices[0], expected);
    EXPECT_EQ(edgeToVertices[1], expected);
    EXPECT_EQ(edgeToVertices[2], expected);
    EXPECT_EQ(edgeToVertices[15], expected);
    EXPECT_EQ(edgeToVertices[91], expected);
    EXPECT_EQ(edgeToVertices[107], expected);
    e += 4;
    // Unencode (1, 5)
    v = 31;
    localUnencoding(vertexToEdgesL, edgeToVerticesL, localVertexMap, logicalOperatorMapX, logicalOperatorMapZ, v, L, engine, false);
    modifyEdgeToVertices(edgeToVertices, edgeToVerticesL, v, L, localVertexMap);
    expected = {24, 32};
    EXPECT_EQ(edgeToVertices[e], expected);
    expected = {1, 24};
    EXPECT_EQ(edgeToVertices[e + 1], expected);
    expected = {2, 25};
    EXPECT_EQ(edgeToVertices[e + 2], expected);
    expected = {2, 30};
    EXPECT_EQ(edgeToVertices[e + 3], expected);
    expected = {0, 0};
    EXPECT_EQ(edgeToVertices[74], expected);
    EXPECT_EQ(edgeToVertices[76], expected);
    EXPECT_EQ(edgeToVertices[90], expected);
    EXPECT_EQ(edgeToVertices[93], expected);
    EXPECT_EQ(edgeToVertices[94], expected);
    EXPECT_EQ(edgeToVertices[95], expected);
}

TEST(buildEdgeToFaces, correctL6)
{
    int L = 6;
    auto edgeToFaces = buildEdgeToFaces(L);
    std::pair<int, int> expected{0, 61};
    EXPECT_EQ(edgeToFaces[0], expected);
    expected = {1, 10};
    EXPECT_EQ(edgeToFaces[1], expected);
    expected = {0, 1};
    EXPECT_EQ(edgeToFaces[2], expected);
    expected = {10, 71};
    EXPECT_EQ(edgeToFaces[15], expected);
    expected = {61, 70};
    EXPECT_EQ(edgeToFaces[91], expected);
    expected = {70, 71};
    EXPECT_EQ(edgeToFaces[107], expected);
    expected = {48, 49};
    EXPECT_EQ(edgeToFaces[74], expected);
    expected = {48, 51};
    EXPECT_EQ(edgeToFaces[76], expected);
    expected = {49, 60};
    EXPECT_EQ(edgeToFaces[90], expected);
    expected = {51, 62};
    EXPECT_EQ(edgeToFaces[93], expected);
    expected = {60, 63};
    EXPECT_EQ(edgeToFaces[94], expected);
    expected = {62, 63};
    EXPECT_EQ(edgeToFaces[95], expected);
}

TEST(modifyLogicalOperators, correct_L6)
{
    int L = 6;
    auto logicals = buildLogicals(L);
    auto l3 = logicals[2];
    // Unencode (0,0)
    int v = 0;
    int e = 3 * L * L;
    vvint vertexToEdgesL;
    vint localVertexMap;
    vpint edgeToVerticesL;
    std::map<std::pair<int, int>, std::vector<int>> logicalOperatorMapX;
    std::map<std::pair<int, int>, std::vector<int>> logicalOperatorMapZ;
    std::mt19937 engine;
    localUnencoding(vertexToEdgesL, edgeToVerticesL, localVertexMap, logicalOperatorMapX, logicalOperatorMapZ, v, L, engine, false);
    modifyLogicalOperators(logicals, logicalOperatorMapX, e);
    e += 4;
    // Unencode (1, 5)
    v = 31;
    localUnencoding(vertexToEdgesL, edgeToVerticesL, localVertexMap, logicalOperatorMapX, logicalOperatorMapZ, v, L, engine, false);
    modifyLogicalOperators(logicals, logicalOperatorMapX, e);
    sint expectedL1 = {22, 23, 32, 33, 42, 43, 52, 53, e - 2, e - 1, e + 2, e + 3};
    sint expectedL2 = {2, 5, 18, 21, 25, 34, 38, 41, 54, 57, e - 2};
    sint expectedL4 = {4, 7, 20, 23, 24, 27, 40, 43, 56, 59, e + 1};
    EXPECT_EQ(logicals[0], expectedL1);
    EXPECT_EQ(logicals[1], expectedL2);
    EXPECT_EQ(logicals[2], l3);
    EXPECT_EQ(logicals[3], expectedL4);
}


// Test unencode for p=1 and p=0 
// e.g. all CC edges removed etc

TEST(unencode, no_unencode)
{
    vint Ls = {3, 6, 9, 12, 15};
    std::random_device rd{}; 
    std::mt19937 engine{rd()};
    std::uniform_real_distribution<double> dist{0.0, 1.0};
    for (auto const L : Ls)
    {
        auto vertexToQubits = buildVertexToFaces(L);
        auto edgeToVertices = buildEdgeToVertices(L);
        auto vertexToEdges = buildVertexToEdges(L);
        auto logicals = buildLogicals(L);
        sint unencodedVertices, qubitIndices;
        auto edgeToFaces = buildEdgeToFaces(L);
        auto faceToEdges = buildFaceToEdges(L);
        vint qubits;
        double p = 0.0;
        auto lift = buildLift(L, vertexToQubits, vertexToEdges, faceToEdges);
        unencode(vertexToQubits, edgeToVertices, unencodedVertices, qubitIndices, qubits, logicals, edgeToFaces, vertexToEdges, L, p, engine, dist, false, true, lift);
        EXPECT_EQ(qubitIndices.size(), 2 * L * L);
        for (int i = 0; i < 2 * L * L; ++i) EXPECT_TRUE(qubitIndices.find(i) != qubitIndices.end());
        EXPECT_EQ(qubits.size(), 3 * L * L);
        EXPECT_EQ(edgeToVertices, buildEdgeToVertices(L));
    }
}

TEST(unencode, full_unencode)
{
    // vint Ls = {3, 6, 9, 12, 15};
    vint Ls = {3};
    std::random_device rd{}; 
    std::mt19937 engine{rd()};
    std::uniform_real_distribution<double> dist{0.0, 1.0};
    for (auto const L : Ls)
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
        // unencodedVertices
        for (int v = 0; v < L * L; ++v)
        {
            if (vertexColor(v, L) == r) EXPECT_TRUE(unencodedVertices.find(v) != unencodedVertices.end());
        }
        // qubitIndices
        EXPECT_EQ(qubitIndices.size(), 4 * L * L / 3);
        for (int eL = 3 * L * L; eL < 3 * L * L + 4 * L * L / 3; ++eL)
        {
            EXPECT_TRUE(qubitIndices.find(eL) != qubitIndices.end());
        }
        // vertexToQubits
        for (auto const &qs : vertexToQubits) EXPECT_EQ(qs.size(), 4);
        // edgeToVertices 
        EXPECT_EQ(edgeToVertices.size(), 3 * L * L + 4 * L * L / 3);
        // WARNING: We only longer delete the 'redundant edges' for r-only unencoding
        std::pair<int, int> defaultPair{0, 0};
        for (int e = 0; e < 3 * L * L; ++e) EXPECT_EQ(edgeToVertices[e], defaultPair); 
        // logicalOperators 
        EXPECT_EQ(logicals[0].size(), 2 * L);
        EXPECT_EQ(logicals[1].size(), L);
        EXPECT_EQ(logicals[2].size(), 2 * L);
        EXPECT_EQ(logicals[3].size(), L);
        EXPECT_EQ(qubits.size(), 3 * L * L + 4 * L * L / 3);
    }
}

TEST(buildVertexToEdges, same_both_ways)
{
    vint Ls = {6, 9, 12, 15};
    for (auto const L : Ls)
    {
        auto vertexToEdgesA = buildVertexToEdges(L);
        auto edgeToVertices = buildEdgeToVertices(L);
        auto vertexToEdgesB = buildVertexToEdges(edgeToVertices, L);
        EXPECT_EQ(vertexToEdgesA.size(), vertexToEdgesB.size());
        for (int v = 0; v < vertexToEdgesA.size(); ++v)
        {
            auto edgesA = vertexToEdgesA[v];
            auto edgesB = vertexToEdgesB[v];
            std::sort(edgesA.begin(), edgesA.end());
            std::sort(edgesB.begin(), edgesB.end());
            EXPECT_EQ(edgesA, edgesB);
        }
    }
}

TEST(pairToEdge, full_unencode)
{
    int L = 6;
    std::random_device rd{}; 
    std::mt19937 engine{rd()};
    std::uniform_real_distribution<double> dist{0.0, 1.0};
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
    int offset = 3 * L * L;
    EXPECT_EQ(pairToEdge(20, 28, vertexToEdges, edgeToVertices), offset + 31);
    EXPECT_EQ(pairToEdge(19, 27, vertexToEdges, edgeToVertices), offset + 32);
}

TEST(ccNeighbors, expected_output)
{
    int L = 10;
    int v = 0;
    auto neighbors = ccNeighbors(v, L);
    EXPECT_TRUE(neighbors.find(1) != neighbors.end()); // +x 
    EXPECT_TRUE(neighbors.find(9) != neighbors.end()); // -x
    EXPECT_TRUE(neighbors.find(11) != neighbors.end()); // +xy
    EXPECT_TRUE(neighbors.find(99) != neighbors.end()); // -xy
    EXPECT_TRUE(neighbors.find(10) != neighbors.end()); // +y
    EXPECT_TRUE(neighbors.find(90) != neighbors.end()); // -y
}

TEST(combinationsUpToK, expected_output)
{
    int n = 6;
    int k = 3;
    auto combs = combinationsUpToK(n, k);
    vvint expectedOut = {
        {0},
        {1},
        {2},
        {3},
        {4},
        {5},
        {0, 1},
        {0, 2},
        {0, 3},
        {0, 4},
        {0, 5},
        {1, 2},
        {1, 3},
        {1, 4},
        {1, 5},
        {2, 3},
        {2, 4},
        {2, 5},
        {3, 4},
        {3, 5},
        {4, 5},
        {0, 1, 2},
        {0, 1, 3},
        {0, 1, 4},
        {0, 1, 5},
        {0, 2, 3},
        {0, 2, 4},
        {0, 2, 5},
        {0, 3, 4},
        {0, 3, 5},
        {0, 4, 5},
        {1, 2, 3},
        {1, 2, 4},
        {1, 2, 5},
        {1, 3, 4},
        {1, 3, 5},
        {1, 4, 5},
        {2, 3, 4},
        {2, 3, 5},
        {2, 4, 5},
        {3, 4, 5},
    };
    // for (auto comb : combs)
    // {
    //     for (auto c : comb)
    //     {
    //         std::cerr << c << " ";
    //     }
    //     std::cerr << std::endl;
    // }
    EXPECT_EQ(combs, expectedOut);
}

TEST(buildLift, expected_output_L6)
{
    int L = 6;
    auto vertexToQubits = buildVertexToFaces(L);
    sint unencodedVertices;
    auto vertexToEdges = buildVertexToEdges(L);
    sint qubitIndices;
    for (int i = 0; i < 2 * L * L; ++i) qubitIndices.insert(i);
    auto faceToEdges = buildFaceToEdges(L);
    auto lift = buildLift(L, vertexToQubits, vertexToEdges, faceToEdges);
    // std::cerr << lift.size() << std::endl;
    // for (auto it = lift.begin(); it != lift.end(); ++it)
    // {
    //     auto edges = it->first;
    //     auto error = it->second;
    //     std::cerr << "Edges : ";
    //     for (auto e : edges) std::cerr << e << " ";
    //     std::cerr << ", Error : ";
    //     for (auto e : error) std::cerr << e << " ";
    //     std::cerr << std::endl;
    // }

    // Vertex (5, 1) = 11
    // Single qubit errors
    vint expectedQs = {9};
    vint edges = {14, 30};
    EXPECT_EQ(lift[edges], expectedQs);
    expectedQs = {8};
    edges = {14, 16};
    EXPECT_EQ(lift[edges], expectedQs);
    expectedQs = {11};
    edges = {16, 33};
    EXPECT_EQ(lift[edges], expectedQs);
    expectedQs = {20};
    edges = {30, 34};
    EXPECT_EQ(lift[edges], expectedQs);
    expectedQs = {22};
    edges = {33, 35};
    EXPECT_EQ(lift[edges], expectedQs);
    expectedQs = {23};
    edges = {34, 35};
    EXPECT_EQ(lift[edges], expectedQs);
    // Two qubit errors
    expectedQs = {8, 9};
    edges = {16, 30};
    EXPECT_EQ(lift[edges], expectedQs);
    expectedQs = {9, 11};
    edges = {14, 16, 30, 33};
    EXPECT_EQ(lift[edges], expectedQs);
    expectedQs = {9, 20};
    edges = {14, 34};
    EXPECT_EQ(lift[edges], expectedQs);
    expectedQs = {9, 22};
    edges = {14, 30, 33, 35};
    EXPECT_EQ(lift[edges], expectedQs);
    expectedQs = {9, 23};
    edges = {14, 30, 34, 35};
    EXPECT_EQ(lift[edges], expectedQs);
    expectedQs = {8, 11};
    edges = {14, 33};
    EXPECT_EQ(lift[edges], expectedQs);
    expectedQs = {8, 20};
    edges = {14, 16, 30, 34};
    EXPECT_EQ(lift[edges], expectedQs);
    expectedQs = {8, 22};
    edges = {14, 16, 33, 35};
    EXPECT_EQ(lift[edges], expectedQs);
    expectedQs = {8, 23};
    edges = {14, 16, 34, 35};
    EXPECT_EQ(lift[edges], expectedQs);
    expectedQs = {11, 20};
    edges = {16, 30, 33, 34};
    EXPECT_EQ(lift[edges], expectedQs);
    expectedQs = {11, 22};
    edges = {16, 35};
    EXPECT_EQ(lift[edges], expectedQs);
    expectedQs = {11, 23};
    edges = {16, 33, 34, 35};
    EXPECT_EQ(lift[edges], expectedQs);
    expectedQs = {20, 22};
    edges = {30, 33, 34, 35};
    EXPECT_EQ(lift[edges], expectedQs);
    expectedQs = {20, 23};
    edges = {30, 35};
    EXPECT_EQ(lift[edges], expectedQs);
    expectedQs = {22, 23};
    edges = {33, 34};
    EXPECT_EQ(lift[edges], expectedQs);
    // Three qubit errors
    expectedQs = {8, 9, 11};
    vint expectedQsAlt = {20, 22, 23}; // The same up to stabilizer
    edges = {30, 33};
    EXPECT_TRUE(lift[edges] == expectedQs || lift[edges] == expectedQsAlt);
    expectedQs = {8, 9, 20};
    expectedQsAlt = {11, 22, 23};
    edges = {16, 34};
    EXPECT_TRUE(lift[edges] == expectedQs || lift[edges] == expectedQsAlt);
    expectedQs = {8, 9, 22};
    expectedQsAlt = {11, 20, 23};
    edges = {16, 30, 33, 35};
    EXPECT_TRUE(lift[edges] == expectedQs || lift[edges] == expectedQsAlt);
    expectedQs = {8, 9, 23};
    expectedQsAlt = {11, 20, 22};
    edges = {16, 30, 34, 35};
    EXPECT_TRUE(lift[edges] == expectedQs || lift[edges] == expectedQsAlt);
    expectedQs = {8, 11, 20};
    expectedQsAlt = {9, 22, 23};
    edges = {14, 30, 33, 34};
    EXPECT_TRUE(lift[edges] == expectedQs || lift[edges] == expectedQsAlt);
    expectedQs = {8, 11, 22};
    expectedQsAlt = {9, 20, 23};
    edges = {14, 35};
    EXPECT_TRUE(lift[edges] == expectedQs || lift[edges] == expectedQsAlt);
    expectedQs = {8, 11, 23};
    expectedQsAlt = {9, 20, 22};
    edges = {14, 33, 34, 35};
    EXPECT_TRUE(lift[edges] == expectedQs || lift[edges] == expectedQsAlt);
    expectedQs = {8, 20, 22};
    expectedQsAlt = {9, 11, 23};
    edges = {14, 16, 30, 33, 34, 35};
    EXPECT_TRUE(lift[edges] == expectedQs || lift[edges] == expectedQsAlt);
    expectedQs = {8, 20, 23};
    expectedQsAlt = {9, 11, 22};
    edges = {14, 16, 30, 35};
    EXPECT_TRUE(lift[edges] == expectedQs || lift[edges] == expectedQsAlt);
    expectedQs = {8, 22, 23};
    expectedQsAlt = {9, 11, 20};
    edges = {14, 16, 33, 34};
    EXPECT_TRUE(lift[edges] == expectedQs || lift[edges] == expectedQsAlt);
}

// TEST(buildLift, empty_lift_full_unencode)
// {
//     int L = 6;
//     auto vertexToQubits = buildVertexToFaces(L);
//     auto edgeToVertices = buildEdgeToVertices(L);
//     auto vertexToEdges = buildVertexToEdges(L);
//     auto logicals = buildLogicals(L);
//     sint unencodedVertices, qubitIndices;
//     auto edgeToFaces = buildEdgeToFaces(L);
//     vint qubits;
//     double p = 1.0;
//     auto faceToEdges = buildFaceToEdges(L);
//     auto lift = buildLift(L, vertexToQubits, vertexToEdges, faceToEdges);
//     unencode(vertexToQubits, edgeToVertices, unencodedVertices, qubitIndices, qubits, logicals, edgeToFaces, vertexToEdges, L, p, engine, dist, false, true, lift);
//     EXPECT_EQ(lift.size(), 0);
// }

TEST(buildLift, L18_problem)
{
    int L = 18;
    auto vertexToQubits = buildVertexToFaces(L);
    sint unencodedVertices;
    auto vertexToEdges = buildVertexToEdges(L);
    sint qubitIndices;
    for (int i = 0; i < 2 * L * L; ++i) qubitIndices.insert(i);
    auto faceToEdges = buildFaceToEdges(L);
    auto lift = buildLift(L, vertexToQubits, vertexToEdges, faceToEdges);

    vint problemEdges = {41, 97};
    EXPECT_NO_THROW(auto qs = lift.at(problemEdges));
    auto qs = lift.at(problemEdges);
    vint expectedQs = {27, 62};
    EXPECT_EQ(qs, expectedQs);
}