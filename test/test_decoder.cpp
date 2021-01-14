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
    sint unencodedVertices, qubitIndices;
    auto edgeToFaces = buildEdgeToFaces(L);
    vint qubits;
    double p = 1.0;
    auto faceToEdges = buildFaceToEdges(L);
    auto lift = buildLift(L, vertexToQubits, vertexToEdges, faceToEdges);
    unencode(vertexToQubits, edgeToVertices, unencodedVertices, qubitIndices, qubits, logicals, edgeToFaces, vertexToEdges, L, p, engine, dist, false, true, lift);
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

TEST(localLift, correct_L6)
{
    int L = 6;
    auto vertexToQubits = buildVertexToFaces(L);
    sint unencodedVertices;
    auto vertexToEdges = buildVertexToEdges(L);
    sint qubitIndices;
    for (int i = 0; i < 2 * L * L; ++i) qubitIndices.insert(i);
    auto faceToEdges = buildFaceToEdges(L);
    auto lift = buildLift(L, vertexToQubits, vertexToEdges, faceToEdges);
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
    graph_t gr = buildGraph(edgeToVertices, b, L);
    std::map<vertex_descriptor, std::vector<vertex_descriptor>> excitationToPathsRG;
    std::map<vertex_descriptor, vint> excitationToDistancesRG;
    vint vertices(L * L);
    std::iota (std::begin(vertices), std::end(vertices), 0); // Populate with 0, 1, ..., (L * L) - 1
    vint verticesRG;
    for (auto const v : vertices)
    {
        if (vertexColor(v, L) != b) verticesRG.push_back(v);
    }
    shortestPaths(gr, excitationToPathsRG, excitationToDistancesRG, verticesRG);
    gr = buildGraph(edgeToVertices, g, L);
    std::map<vertex_descriptor, std::vector<vertex_descriptor>> excitationToPathsRB;
    std::map<vertex_descriptor, vint> excitationToDistancesRB;
    vint verticesRB;
    for (auto const v : vertices)
    {
        if (vertexColor(v, L) != g) verticesRB.push_back(v);
    }
    shortestPaths(gr, excitationToPathsRB, excitationToDistancesRB, verticesRB);
    vint paths, redVertices;
    // r-r
    vint excitations = {16, 26};
    color c = b; 
    pairExcitations(excitations, edgeToVertices, vertexToEdges, c, L, paths, redVertices, excitationToPathsRG, excitationToDistancesRG);
    std::sort(paths.begin(), paths.end());
    vint expectedPath = {49, 63, 64, 78};
    EXPECT_EQ(paths, expectedPath);
    vint expectedReds = {16, 21, 26};
    EXPECT_EQ(redVertices, expectedReds);
    c = g;
    paths = {};
    redVertices = {};
    pairExcitations(excitations, edgeToVertices, vertexToEdges, c, L, paths, redVertices, excitationToPathsRB, excitationToDistancesRB);
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
    pairExcitations(excitations, edgeToVertices, vertexToEdges, c, L, paths, redVertices, excitationToPathsRG, excitationToDistancesRG);
    std::sort(paths.begin(), paths.end());
    expectedPath = {9, 10};
    EXPECT_EQ(paths, expectedPath);
    expectedReds = {3};
    EXPECT_EQ(redVertices, expectedReds);
    c = g;
    paths = {};
    redVertices = {};
    pairExcitations(excitations, edgeToVertices, vertexToEdges, c, L, paths, redVertices, excitationToPathsRB, excitationToDistancesRB);
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
    pairExcitations(excitations, edgeToVertices, vertexToEdges, c, L, paths, redVertices, excitationToPathsRG, excitationToDistancesRG);
    std::sort(paths.begin(), paths.end());
    expectedPath = {74};
    EXPECT_EQ(paths, expectedPath);
    expectedReds = {31};
    EXPECT_EQ(redVertices, expectedReds);
    c = g;
    paths = {};
    redVertices = {};
    pairExcitations(excitations, edgeToVertices, vertexToEdges, c, L, paths, redVertices, excitationToPathsRB, excitationToDistancesRB);
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
    sint unencodedVertices, qubitIndices;
    auto edgeToFaces = buildEdgeToFaces(L);
    vint qubits;
    double p = 1.0;
    auto faceToEdges = buildFaceToEdges(L);
    auto lift = buildLift(L, vertexToQubits, vertexToEdges, faceToEdges);
    unencode(vertexToQubits, edgeToVertices, unencodedVertices, qubitIndices, qubits, logicals, edgeToFaces, vertexToEdges, L, p, engine, dist, false, true, lift);
    graph_t gr = buildGraph(edgeToVertices, b, L);
    std::map<vertex_descriptor, std::vector<vertex_descriptor>> excitationToPathsRG;
    std::map<vertex_descriptor, vint> excitationToDistancesRG;
    vint vertices(L * L);
    std::iota (std::begin(vertices), std::end(vertices), 0); // Populate with 0, 1, ..., (L * L) - 1
    vint verticesRG;
    for (auto const v : vertices)
    {
        if (vertexColor(v, L) != b) verticesRG.push_back(v);
    }
    shortestPaths(gr, excitationToPathsRG, excitationToDistancesRG, verticesRG);
    gr = buildGraph(edgeToVertices, g, L);
    std::map<vertex_descriptor, std::vector<vertex_descriptor>> excitationToPathsRB;
    std::map<vertex_descriptor, vint> excitationToDistancesRB;
    vint verticesRB;
    for (auto const v : vertices)
    {
        if (vertexColor(v, L) != g) verticesRB.push_back(v);
    }
    shortestPaths(gr, excitationToPathsRB, excitationToDistancesRB, verticesRB);
    // b-b & g-g
    vint excitations = {25, 33, 14, 22};
    vint paths, redVertices;
    color c = b; 
    pairExcitations(excitations, edgeToVertices, vertexToEdges, c, L, paths, redVertices, excitationToPathsRG, excitationToDistancesRG);
    std::sort(paths.begin(), paths.end());
    vint expectedPath = {136};
    EXPECT_EQ(paths, expectedPath);
    EXPECT_EQ(redVertices.size(), 0);
    c = g;
    paths = {};
    redVertices = {};
    pairExcitations(excitations, edgeToVertices, vertexToEdges, c, L, paths, redVertices, excitationToPathsRB, excitationToDistancesRB);
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
    auto lift = buildLift(L, vertexToQubits, vertexToEdges, faceToEdges);
    graph_t gr = buildGraph(edgeToVertices, b, L);
    std::map<vertex_descriptor, std::vector<vertex_descriptor>> excitationToPathsRG;
    std::map<vertex_descriptor, vint> excitationToDistancesRG;
    vint vertices(L * L);
    std::iota (std::begin(vertices), std::end(vertices), 0); // Populate with 0, 1, ..., (L * L) - 1
    vint verticesRG;
    for (auto const v : vertices)
    {
        if (vertexColor(v, L) != b) verticesRG.push_back(v);
    }
    shortestPaths(gr, excitationToPathsRG, excitationToDistancesRG, verticesRG);
    gr = buildGraph(edgeToVertices, g, L);
    std::map<vertex_descriptor, std::vector<vertex_descriptor>> excitationToPathsRB;
    std::map<vertex_descriptor, vint> excitationToDistancesRB;
    vint verticesRB;
    for (auto const v : vertices)
    {
        if (vertexColor(v, L) != g) verticesRB.push_back(v);
    }
    shortestPaths(gr, excitationToPathsRB, excitationToDistancesRB, verticesRB);

    // one qubit error
    vint excitations = {24, 30, 31};
    auto outCorr = findCorrection(excitations, edgeToVertices, vertexToEdges, L, lift, unencodedVertices, excitationToPathsRG, excitationToDistancesRG, excitationToPathsRB, excitationToDistancesRB);
    sint expectedCorr = {49};
    EXPECT_EQ(outCorr, expectedCorr);

    // 2x two qubit errors 
    excitations = {2, 7, 4, 9};
    outCorr = findCorrection(excitations, edgeToVertices, vertexToEdges, L, lift, unencodedVertices, excitationToPathsRG, excitationToDistancesRG, excitationToPathsRB, excitationToDistancesRB);
    expectedCorr = {2, 3, 6, 7};
    EXPECT_EQ(outCorr, expectedCorr);

    // multi-qubit error
    excitations = {16, 26};
    outCorr = findCorrection(excitations, edgeToVertices, vertexToEdges, L, lift, unencodedVertices, excitationToPathsRG, excitationToDistancesRG, excitationToPathsRB, excitationToDistancesRB);
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
    sint unencodedVertices, qubitIndices;
    auto edgeToFaces = buildEdgeToFaces(L);
    vint qubits;
    double p = 1.0;
    auto faceToEdges = buildFaceToEdges(L);
    auto lift = buildLift(L, vertexToQubits, vertexToEdges, faceToEdges);
    unencode(vertexToQubits, edgeToVertices, unencodedVertices, qubitIndices, qubits, logicals, edgeToFaces, vertexToEdges, L, p, engine, dist, false, true, lift);
    graph_t gr = buildGraph(edgeToVertices, b, L);
    std::map<vertex_descriptor, std::vector<vertex_descriptor>> excitationToPathsRG;
    std::map<vertex_descriptor, vint> excitationToDistancesRG;
    vint vertices(L * L);
    std::iota (std::begin(vertices), std::end(vertices), 0); // Populate with 0, 1, ..., (L * L) - 1
    vint verticesRG;
    for (auto const v : vertices)
    {
        if (vertexColor(v, L) != b) verticesRG.push_back(v);
    }
    shortestPaths(gr, excitationToPathsRG, excitationToDistancesRG, verticesRG);
    gr = buildGraph(edgeToVertices, g, L);
    std::map<vertex_descriptor, std::vector<vertex_descriptor>> excitationToPathsRB;
    std::map<vertex_descriptor, vint> excitationToDistancesRB;
    vint verticesRB;
    for (auto const v : vertices)
    {
        if (vertexColor(v, L) != g) verticesRB.push_back(v);
    }
    shortestPaths(gr, excitationToPathsRB, excitationToDistancesRB, verticesRB);

    // 2x one qubit errors (b and g)
    vint excitations = {15, 23, 19, 32};
    auto outCorr = findCorrection(excitations, edgeToVertices, vertexToEdges, L, lift, unencodedVertices, excitationToPathsRG, excitationToDistancesRG, excitationToPathsRB, excitationToDistancesRB);
    sint expectedCorr = {131, 141};
    EXPECT_EQ(outCorr, expectedCorr);
}

TEST(checkCorrection, no_unencode)
{
    int L = 6;
    auto logicals = buildLogicals(L);
    vint qubits(2 * L * L, 0);
    // Stabilizer error
    // Stab with index 15
    qubits[30] = 1;
    qubits[31] = 1;
    qubits[16] = 1;
    qubits[17] = 1;
    qubits[19] = 1;
    qubits[28] = 1;
    EXPECT_TRUE(checkCommutation(qubits, logicals));
    // L1 g 
    vint indices = {0, 1, 22, 23, 32, 33, 42, 43, 52, 54, 62, 63};
    for (auto i : indices)
    {
        qubits[i] = (qubits[i] + 1) % 2;
    }
    EXPECT_FALSE(checkCommutation(qubits, logicals));
    // L1 b
    indices = {2, 3, 12, 13, 34, 35, 44, 45, 54, 55, 64, 65};
    for (auto i : indices)
    {
        qubits[i] = (qubits[i] + 1) % 2;
    }
    EXPECT_FALSE(checkCommutation(qubits, logicals));
    // L2 g
    indices = {2, 5, 18, 21, 34, 35, 38, 41, 54, 57, 70, 61};
    for (auto i : indices)
    {
        qubits[i] = (qubits[i] + 1) % 2;
    }
    EXPECT_FALSE(checkCommutation(qubits, logicals));
    // L2 b
    indices = {1, 14, 17, 30, 33, 46, 37, 50, 53, 66, 69};
    for (auto i : indices)
    {
        qubits[i] = (qubits[i] + 1) % 2;
    }
    EXPECT_FALSE(checkCommutation(qubits, logicals));
    // L red diagonal
    std::fill(qubits.begin(), qubits.end(), 0);
    indices = {10, 11, 20, 21, 30, 31, 40, 41, 50, 51, 60, 61};
    for (auto i : indices)
    {
        qubits[i] = (qubits[i] + 1) % 2;
    }
    EXPECT_FALSE(checkCommutation(qubits, logicals));
}

TEST(checkCorrection, full_unencode)
{
    int L = 6;
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
    // Stabilizer error b
    int offset = 3 * L * L;
    qubits[offset + 30] = 1;
    qubits[offset + 31] = 1;
    qubits[offset + 11] = 1;
    qubits[offset + 18] = 1;
    EXPECT_TRUE(checkCommutation(qubits, logicals));
    // Stabilizer error g
    qubits[offset + 17] = 1;
    qubits[offset + 16] = 1;
    qubits[offset + 32] = 1;
    qubits[offset + 29] = 1;
    EXPECT_TRUE(checkCommutation(qubits, logicals));
    // L1 b
    vint indices = {35, 6, 15, 26};
    for (auto i : indices)
    {
        qubits[offset + i] = (qubits[offset + i] + 1) % 2;
    }
    EXPECT_FALSE(checkCommutation(qubits, logicals));
    // L2 b
    indices = {2, 18, 34, 6, 22, 38};
    for (auto i : indices)
    {
        qubits[offset + i] = (qubits[offset + i] + 1) % 2;
    }
    EXPECT_FALSE(checkCommutation(qubits, logicals));
    // L1 g
    indices = {9, 28, 37, 0};
    for (auto i : indices)
    {
        qubits[offset + i] = (qubits[offset + i] + 1) % 2;
    }
    EXPECT_FALSE(checkCommutation(qubits, logicals));
    // L2 g
    indices = {16, 28, 36, 40, 4, 12};
    for (auto i : indices)
    {
        qubits[offset + i] = (qubits[offset + i] + 1) % 2;
    }
    EXPECT_FALSE(checkCommutation(qubits, logicals));
}

TEST(reRoute, correct_L6)
{
    // Initialization
    int L = 6;
    auto vertexToQubits = buildVertexToFaces(L);
    auto edgeToVertices = buildEdgeToVertices(L);
    auto vertexToEdges = buildVertexToEdges(L);
    auto logicals = buildLogicals(L);
    sint unencodedVertices, qubitIndices;
    auto edgeToFaces = buildEdgeToFaces(L);
    vint qubits;

    // This is basically the unencoding function
    for (int i = 0; i < 2 * L * L; ++i) qubitIndices.insert(i);
    int e = 3 * L * L;
    vint vertices = {20, 22}; // vertices to unencode in the test
    for (auto const v : vertices)
    {
        auto neighV = ccNeighbors(v, L);
        for (auto u : neighV)
        {
            // If a neighboring vertex has already been unencoded we can't unencode
            if (unencodedVertices.find(u) != unencodedVertices.end()) continue; 
        }
        unencodedVertices.insert(v);
        // qubitIndices
        auto faces = vertexToQubits[v];
        for (auto f : faces) qubitIndices.erase(f);
        qubitIndices.insert(e);
        qubitIndices.insert(e + 1);
        qubitIndices.insert(e + 2);
        qubitIndices.insert(e + 3);
        // Specify unencoding locally 
        vvint vertexToEdgesL;
        vpint edgeToVerticesL;
        vint localVertexMap;
        std::map<std::pair<int, int>, std::vector<int>> logicalOperatorMapX;
        std::map<std::pair<int, int>, std::vector<int>> logicalOperatorMapZ;
        localUnencoding(vertexToEdgesL, edgeToVerticesL, localVertexMap, logicalOperatorMapX, logicalOperatorMapZ, v, L, engine, false); // Use default unencoding
        // Modify data structures to account for unencoding
        modifyVertexToQubits(vertexToQubits, vertexToEdgesL, localVertexMap, v, e, L);
        modifyEdgeToVertices(edgeToVertices, edgeToVerticesL, v, L, localVertexMap);
        modifyLogicalOperators(logicals, logicalOperatorMapX, e);
        e += 4;
    }
    qubits.assign(e, 0); // Init qubits with length equal to final edge index + 1
    vertexToEdges = buildVertexToEdges(edgeToVertices, L); // Use overloaded function to invert edgeToVertices

    auto newPath = reRoute(13, 26, 20, b, L, vertexToEdges, edgeToVertices);
    std::sort(newPath.begin(), newPath.end());
    vvint expPaths = {{39, 78, 110}, {40, 78, 111}};
    EXPECT_TRUE(newPath == expPaths[0] || newPath == expPaths[1]);
    newPath = reRoute(26, 13, 20, b, L, vertexToEdges, edgeToVertices);
    std::sort(newPath.begin(), newPath.end());
    EXPECT_TRUE(newPath == expPaths[0] || newPath == expPaths[1]);

    newPath = reRoute(21, 13, 20, b, L, vertexToEdges, edgeToVertices);
    std::sort(newPath.begin(), newPath.end());
    expPaths = {{39, 64, 110}, {40, 64, 111}};
    EXPECT_TRUE(newPath == expPaths[0] || newPath == expPaths[1]);
    newPath = reRoute(13, 21, 20, b, L, vertexToEdges, edgeToVertices);
    std::sort(newPath.begin(), newPath.end());
    EXPECT_TRUE(newPath == expPaths[0] || newPath == expPaths[1]);
    
    newPath = reRoute(16, 29, 22, g, L, vertexToEdges, edgeToVertices);
    std::sort(newPath.begin(), newPath.end());
    expPaths = {{45, 70, 112}, {45, 84, 113}};
    EXPECT_TRUE(newPath == expPaths[0] || newPath == expPaths[1]);
    newPath = reRoute(29, 16, 22, g, L, vertexToEdges, edgeToVertices);
    std::sort(newPath.begin(), newPath.end());
    EXPECT_TRUE(newPath == expPaths[0] || newPath == expPaths[1]);

    newPath = reRoute(21, 29, 22, g, L, vertexToEdges, edgeToVertices);
    std::sort(newPath.begin(), newPath.end());
    expPaths = {{46, 70, 112}, {46, 84, 113}};
    EXPECT_TRUE(newPath == expPaths[0] || newPath == expPaths[1]);
    newPath = reRoute(29, 21, 22, g, L, vertexToEdges, edgeToVertices);
    std::sort(newPath.begin(), newPath.end());
    EXPECT_TRUE(newPath == expPaths[0] || newPath == expPaths[1]);
}