#include "gtest/gtest.h"
#include "djikstra.h"
#include "hexLattice.h"
#include "decoder.h"
#include <numeric>

std::random_device rd{}; 
std::mt19937 engine{rd()};
std::uniform_real_distribution<double> dist{0.0, 1.0};

TEST(calcSyndrome, expected_behaviour_no_unencoding)
{
    int L = 6;
    auto vertexToFaces = buildVertexToFaces(L);
    vint error(2 * L * L, 0);
    sint unencodedVertices;
    // vint syndrome(L * L, 0);
    vint syndrome;
    error[0] = 1;
    calcSyndrome(syndrome, vertexToFaces, error, unencodedVertices, L);
    // for (int i = 0; i < syndrome.size(); ++i)
    // {
    //     if (i == 0 || i == 1 || i == 7)
    //     {
    //         EXPECT_EQ(syndrome[i], 1);
    //     }
    //     else
    //     {
    //         EXPECT_EQ(syndrome[i], 0);
    //     }
    // }
    vint expected = {0, 1, 7};
    EXPECT_EQ(syndrome, expected);
    error[3] = 1;
    calcSyndrome(syndrome, vertexToFaces, error, unencodedVertices, L);
    // for (int i = 0; i < syndrome.size(); ++i)
    // {
    //     if (i == 0 || i == 8)
    //     {
    //         EXPECT_EQ(syndrome[i], 1);
    //     }
    //     else
    //     {
    //         EXPECT_EQ(syndrome[i], 0);
    //     }
    // }
    expected = {0, 8};
    EXPECT_EQ(syndrome, expected);
    error[1] = 1;
    error[12] = 1;
    error[14] = 1;
    error[15] = 1;
    calcSyndrome(syndrome, vertexToFaces, error, unencodedVertices, L);
    // for (int i = 0; i < syndrome.size(); ++i)
    // {
    //     EXPECT_EQ(syndrome[i], 0);
    // }
    EXPECT_EQ(syndrome.size(), 0);
}

TEST(calcSyndrome, expected_behaviour_full_unencoding)
{
    int L = 6;
    auto vertexToQubits = buildVertexToFaces(L);
    auto edgeToVertices = buildEdgeToVertices(L);
    auto vertexToEdges = buildVertexToEdges(L);
    auto logicals = buildLogicals(L);
    std::map<std::pair<int, int>, int> vertexPairToEdgeU;
    sint unencodedVertices, qubitIndices;
    auto edgeToFaces = buildEdgeToFaces(L);
    vint qubits;
    double p = 1.0;
    unencode(vertexPairToEdgeU, vertexToQubits, edgeToVertices, unencodedVertices, qubitIndices, qubits, logicals, edgeToFaces, vertexToEdges, L, p, engine, dist);
    // vint syndrome(L * L, 0);
    vint syndrome;
    // Errors on CC qubits should have no effect
    qubits[1] = 1;
    qubits[12] = 1;
    qubits[14] = 1;
    qubits[15] = 1;
    calcSyndrome(syndrome, vertexToQubits, qubits, unencodedVertices, L);
    // for (int i = 0; i < syndrome.size(); ++i)
    // {
    //     EXPECT_EQ(syndrome[i], 0);
    // }
    EXPECT_EQ(syndrome.size(), 0);
    // TC green edge error
    qubits[108] = 1;
    calcSyndrome(syndrome, vertexToQubits, qubits, unencodedVertices, L);
    // for (int i = 0; i < syndrome.size(); ++i)
    // {
    //     if (i == 1 || i == 35)
    //     {
    //         EXPECT_EQ(syndrome[i], 1);
    //     }
    //     else
    //     {
    //         EXPECT_EQ(syndrome[i], 0);
    //     }
    // }
    vint expected = {1, 35};
    EXPECT_EQ(syndrome, expected);
    // for (auto uv : unencodedVertices, L) std::cerr << uv << " ";
    // std::cerr << std::endl;
    // for (auto q : vertexToQubits[0]) std::cerr << q << " ";
    // std::cerr << std::endl;
    // TC blue edge error
    qubits[114] = 1;
    calcSyndrome(syndrome, vertexToQubits, qubits, unencodedVertices, L);
    // for (int i = 0; i < syndrome.size(); ++i)
    // {
    //     if (i == 1 || i == 35 || i == 10 || i == 33)
    //     {
    //         EXPECT_EQ(syndrome[i], 1);
    //     }
    //     else
    //     {
    //         EXPECT_EQ(syndrome[i], 0);
    //     }
    // }
    expected = {1, 10, 33, 35};
    EXPECT_EQ(syndrome, expected);
    // No non-trivial syndrome at red vertices
    p = 0.5;
    // for (auto q : qubitIndices) std::cerr << q << " ";
    // std::cerr << std::endl;
    generateError(qubits, qubitIndices, p, engine, dist);
    calcSyndrome(syndrome, vertexToQubits, qubits, unencodedVertices, L);
    for (int i = 0; i < L * L; ++i)
    {
        if (vertexColor(i, L) == r)
        {
            EXPECT_TRUE(std::find(syndrome.begin(), syndrome.end(), i) == syndrome.end());
        }
    }
}

// pairExcitations

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

TEST(binToList, expected_behaviour)
{
    vint bin(10, 0);
    bin[0] = 1;
    bin[3] = 1;
    bin[9] = 1;
    auto list = binToList(bin);
    vint expectedOut = {0, 3, 9};
    EXPECT_EQ(list, expectedOut);
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
    auto lift = buildLift(L, vertexToQubits, unencodedVertices, vertexToEdges, qubitIndices, faceToEdges);
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

TEST(buildLift, empty_lift_full_unencode)
{
    int L = 6;
    auto vertexToQubits = buildVertexToFaces(L);
    auto edgeToVertices = buildEdgeToVertices(L);
    auto vertexToEdges = buildVertexToEdges(L);
    auto logicals = buildLogicals(L);
    std::map<std::pair<int, int>, int> vertexPairToEdgeU;
    sint unencodedVertices, qubitIndices;
    auto edgeToFaces = buildEdgeToFaces(L);
    vint qubits;
    double p = 1.0;
    unencode(vertexPairToEdgeU, vertexToQubits, edgeToVertices, unencodedVertices, qubitIndices, qubits, logicals, edgeToFaces, vertexToEdges, L, p, engine, dist);
    auto faceToEdges = buildFaceToEdges(L);
    auto lift = buildLift(L, vertexToQubits, unencodedVertices, vertexToEdges, qubitIndices, faceToEdges);
    EXPECT_EQ(lift.size(), 0);
}

TEST(localLift, correct_L6)
{
    int L = 6;
    auto vertexToQubits = buildVertexToFaces(L);
    sint unencodedVertices;
    auto vertexToEdges = buildVertexToEdges(L);
    sint qubitIndices;
    for (int i = 0; i < 2 * L * L; ++i) qubitIndices.insert(i);
    auto faceToEdges = buildFaceToEdges(L);
    auto lift = buildLift(L, vertexToQubits, unencodedVertices, vertexToEdges, qubitIndices, faceToEdges);
    int v = 11;
    vint paths = {0, 2, 33, 34};
    auto out = localLift(v, L, paths, vertexToEdges, lift);
    vint expected = {22, 23};
    EXPECT_EQ(out, expected);
    v = 0;
    out = localLift(v, L, paths, vertexToEdges, lift);
    expected = {0};
    EXPECT_EQ(out, expected);
}

TEST(generateError, p0_no_errors)
{
    double p = 0;
    int L = 12;
    sint qubitIndices;
    for (int i = 0; i < 2 * L * L; ++i) qubitIndices.insert(qubitIndices.end(), i);
    vint qubits(2 * L * L, 0);
    generateError(qubits, qubitIndices, p, engine, dist);
    for (int i = 0; i < qubits.size(); ++i) EXPECT_EQ(qubits[i], 0);
}

TEST(generateError, p1_all_errors)
{
    double p = 1;
    int L = 15;
    sint qubitIndices;
    for (int i = 0; i < 2 * L * L; ++i) qubitIndices.insert(qubitIndices.end(), i);
    vint qubits(2 * L * L, 0);
    generateError(qubits, qubitIndices, p, engine, dist);
    for (int i = 0; i < qubits.size(); ++i) EXPECT_EQ(qubits[i], 1);
}

TEST(pairExcitations, no_unencoding)
{
    int L = 6;
    auto edgeToVertices = buildEdgeToVertices(L);
    auto vertexToEdges = buildVertexToEdges(L);
    vint paths, redVertices;
    // r-r
    vint excitations = {16, 26};
    color c = b; 
    pairExcitations(excitations, edgeToVertices, vertexToEdges, c, L, paths, redVertices);
    std::sort(paths.begin(), paths.end());
    vint expectedPath = {49, 63, 64, 78};
    EXPECT_EQ(paths, expectedPath);
    vint expectedReds = {16, 21, 26};
    EXPECT_EQ(redVertices, expectedReds);
    c = g;
    paths = {};
    redVertices = {};
    pairExcitations(excitations, edgeToVertices, vertexToEdges, c, L, paths, redVertices);
    std::sort(paths.begin(), paths.end());
    expectedPath = {50, 56, 69, 75};
    EXPECT_EQ(paths, expectedPath);
    expectedReds = {16, 18, 26};
    EXPECT_EQ(redVertices, expectedReds);
    // g-g & b-b
    excitations = {4, 9, 2, 7};
    c = b;
    paths = {};
    redVertices = {};
    pairExcitations(excitations, edgeToVertices, vertexToEdges, c, L, paths, redVertices);
    std::sort(paths.begin(), paths.end());
    expectedPath = {9, 10};
    EXPECT_EQ(paths, expectedPath);
    expectedReds = {3};
    EXPECT_EQ(redVertices, expectedReds);
    c = g;
    paths = {};
    redVertices = {};
    pairExcitations(excitations, edgeToVertices, vertexToEdges, c, L, paths, redVertices);
    std::sort(paths.begin(), paths.end());
    expectedPath = {7, 21};
    EXPECT_EQ(paths, expectedPath);
    expectedReds = {8};
    EXPECT_EQ(redVertices, expectedReds);
    // r-b-g
    excitations = {24, 30, 31};
    c = b;
    paths = {};
    redVertices = {};
    pairExcitations(excitations, edgeToVertices, vertexToEdges, c, L, paths, redVertices);
    std::sort(paths.begin(), paths.end());
    expectedPath = {74};
    EXPECT_EQ(paths, expectedPath);
    expectedReds = {31};
    EXPECT_EQ(redVertices, expectedReds);
    c = g;
    paths = {};
    redVertices = {};
    pairExcitations(excitations, edgeToVertices, vertexToEdges, c, L, paths, redVertices);
    std::sort(paths.begin(), paths.end());
    expectedPath = {90};
    EXPECT_EQ(paths, expectedPath);
    expectedReds = {31};
    EXPECT_EQ(redVertices, expectedReds);
}

TEST(pairExcitations, fully_unencoded)
{
    int L = 6;
    auto vertexToQubits = buildVertexToFaces(L);
    auto edgeToVertices = buildEdgeToVertices(L);
    auto vertexToEdges = buildVertexToEdges(L);
    auto logicals = buildLogicals(L);
    std::map<std::pair<int, int>, int> vertexPairToEdgeU;
    sint unencodedVertices, qubitIndices;
    auto edgeToFaces = buildEdgeToFaces(L);
    vint qubits;
    double p = 1.0;
    unencode(vertexPairToEdgeU, vertexToQubits, edgeToVertices, unencodedVertices, qubitIndices, qubits, logicals, edgeToFaces, vertexToEdges, L, p, engine, dist);
    // b-b & g-g
    vint excitations = {25, 33, 14, 22};
    vint paths, redVertices;
    color c = b; 
    pairExcitations(excitations, edgeToVertices, vertexToEdges, c, L, paths, redVertices);
    std::sort(paths.begin(), paths.end());
    vint expectedPath = {136};
    EXPECT_EQ(paths, expectedPath);
    EXPECT_EQ(redVertices.size(), 0);
    c = g;
    paths = {};
    redVertices = {};
    pairExcitations(excitations, edgeToVertices, vertexToEdges, c, L, paths, redVertices);
    std::sort(paths.begin(), paths.end());
    expectedPath = {143};
    EXPECT_EQ(paths, expectedPath);
    EXPECT_EQ(redVertices.size(), 0);
}

TEST(findCorrection, no_unencoding)
{
    int L = 6;
    auto edgeToVertices = buildEdgeToVertices(L);
    auto vertexToEdges = buildVertexToEdges(L);
    auto vertexToQubits = buildVertexToFaces(L);
    sint unencodedVertices, qubitIndices;
    for (int i = 0; i < 2 * L * L; ++i) qubitIndices.insert(i);
    auto faceToEdges = buildFaceToEdges(L);
    auto lift = buildLift(L, vertexToQubits, unencodedVertices, vertexToEdges, qubitIndices, faceToEdges);

    // one qubit error
    vint excitations = {24, 30, 31};
    auto outCorr = findCorrection(excitations, edgeToVertices, vertexToEdges, L, lift);
    vint expectedCorr = {49};
    EXPECT_EQ(outCorr, expectedCorr);

    // 2x two qubit errors 
    excitations = {2, 7, 4, 9};
    outCorr = findCorrection(excitations, edgeToVertices, vertexToEdges, L, lift);
    std::sort(outCorr.begin(), outCorr.end());
    expectedCorr = {2, 3, 6, 7};
    EXPECT_EQ(outCorr, expectedCorr);

    // multi-qubit error
    excitations = {16, 26};
    outCorr = findCorrection(excitations, edgeToVertices, vertexToEdges, L, lift);
    std::sort(outCorr.begin(), outCorr.end());
    expectedCorr = {33, 37, 38, 39, 41, 42, 43, 46};
    EXPECT_EQ(outCorr, expectedCorr);
}

TEST(findCorrection, full_unencode)
{
    int L = 6;
    auto vertexToQubits = buildVertexToFaces(L);
    auto edgeToVertices = buildEdgeToVertices(L);
    auto vertexToEdges = buildVertexToEdges(L);
    auto logicals = buildLogicals(L);
    std::map<std::pair<int, int>, int> vertexPairToEdgeU;
    sint unencodedVertices, qubitIndices;
    auto edgeToFaces = buildEdgeToFaces(L);
    vint qubits;
    double p = 1.0;
    unencode(vertexPairToEdgeU, vertexToQubits, edgeToVertices, unencodedVertices, qubitIndices, qubits, logicals, edgeToFaces, vertexToEdges, L, p, engine, dist);
    auto faceToEdges = buildFaceToEdges(L);
    auto lift = buildLift(L, vertexToQubits, unencodedVertices, vertexToEdges, qubitIndices, faceToEdges);

    // 2x one qubit errors (b and g)
    vint excitations = {15, 23, 19, 32};
    auto outCorr = findCorrection(excitations, edgeToVertices, vertexToEdges, L, lift);
    std::sort(outCorr.begin(), outCorr.end());
    vint expectedCorr = {131, 141};
    EXPECT_EQ(outCorr, expectedCorr);
}
