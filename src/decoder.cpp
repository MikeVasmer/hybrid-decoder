#include "decoder.h"

// void calcSyndrome(vint &syndrome, const vsint &vertexToQubits, const vint &qubits, const sint &undencodedVertices)
// {
//     std::fill(syndrome.begin(), syndrome.end(), 0);
//     for (int i = 0; i < syndrome.size(); ++i)
//     {
//         if (undencodedVertices.find(i) != undencodedVertices.end()) continue;
//         auto support = vertexToQubits[i];
//         int parity = 0;
//         for (auto q : support)
//         {
//             parity = (parity + qubits[q]) % 2;
//         }
//         syndrome[i] = parity;
//     }
// }

// We use excitations later as a list not a bin so better calculate it that way!
void calcSyndrome(vint &syndrome, const vsint &vertexToQubits, const vint &qubits, const sint &undencodedVertices, int L)
{
    syndrome = {};
    for (int i = 0; i < L * L; ++i)
    {
        if (undencodedVertices.find(i) != undencodedVertices.end()) continue;
        auto support = vertexToQubits[i];
        int parity = 0;
        for (auto q : support)
        {
            parity = (parity + qubits[q]) % 2;
        }
        if (parity == 1) syndrome.push_back(i);
    }
}

// Color c to exclude!
void pairExcitations(const vint &excitations, const vpint &edgeToVertices, const vvint &vertexToEdges, color c, int L, vint &paths, vint &redVertices)
{
    graph_t g = buildGraph(edgeToVertices, c, L);
    vint restrictedExcitations;
    for (auto const e : excitations)
    {
        if (vertexColor(e, L) != c) restrictedExcitations.push_back(e);
    }
    std::map<vertex_descriptor, std::vector<vertex_descriptor>> excitationToPaths;
    std::map<vertex_descriptor, vint> excitationToDistances;
    shortestPaths(g, excitationToPaths, excitationToDistances, restrictedExcitations);
    // Blossom V
    vint edges;
    vint weights;
    int nodeNum = restrictedExcitations.size();
    for (int i = 0; i < nodeNum; ++i)
    {
        for (int j = i + 1; j < nodeNum; ++j)
        {
            edges.push_back(i);
            edges.push_back(j);
            weights.push_back(excitationToDistances[restrictedExcitations[i]][restrictedExcitations[j]]);
        }
    }
    int edgeNum = edges.size() / 2;
    struct PerfectMatching::Options options;
    options.verbose = false;
    PerfectMatching *pm;
    pm = new PerfectMatching(nodeNum, edgeNum);
    for (int e = 0; e < edgeNum; ++e)
    {
        pm->AddEdge(edges[2*e], edges[2*e+1], weights[e]);
    }
    pm->options = options;
    pm->Solve();
    for (int i = 0; i < nodeNum; ++i)
    {
        int j = pm->GetMatch(i);
        if (i < j)
        {
            int v_tmp = restrictedExcitations[i];
            auto pathsToJ = excitationToPaths[restrictedExcitations[j]];
            while (v_tmp != restrictedExcitations[j])
            {
                paths.push_back(pairToEdge(v_tmp, pathsToJ[v_tmp], vertexToEdges, edgeToVertices));
                if (vertexColor(v_tmp, L) == r) redVertices.push_back(v_tmp);
                v_tmp = pathsToJ[v_tmp];
            }
            if (vertexColor(restrictedExcitations[j], L) == r) redVertices.push_back(restrictedExcitations[j]);
        }
    }
    delete pm;
    // There should be no duplicate edges in paths
    sint testSet(paths.begin(), paths.end());
    if (testSet.size() != paths.size()) throw std::length_error("Matching paths have duplicate edges!");
}

vvint combinationsUpToK(int n, int k)
{
    vvint combinations;
    for (int j = 1; j <= k; ++j)
    {
        std::vector<bool> selector(n);
        std::fill(selector.begin(), selector.begin() + j, true);
        do {
            vint comb;    
            for (int i = 0; i < n; ++i) {
                if (selector[i]) comb.push_back(i);
            }
            combinations.push_back(comb);
        } while (std::prev_permutation(selector.begin(), selector.end()));
    }
    return combinations;
}

vint binToList(vint &bin)
{
    vint list;
    for (int i = 0; i < bin.size(); ++i)
    {
        if (bin[i] == 1) list.push_back(i);
    }
    return list;
}

// Build a map from local syndrome edges to local errors (r vertices only)
std::map<vint, vint> buildLift(int L, const vsint &vertexToQubits, const sint &unencodedVertices, const vvint &vertexToEdges, const sint &qubitIndices, const vvint &faceToEdges)
{
    vvint indexCombinations = combinationsUpToK(6, 3);
    std::map<vint, vint> lift;
    for (int v = 0; v < L * L; ++v)
    {
        if (vertexColor(v, L) != r) continue;
        if (unencodedVertices.find(v) != unencodedVertices.end()) continue;
        vint qs(vertexToQubits[v].begin(), vertexToQubits[v].end());
        auto edges = vertexToEdges[v];
        // std::sort(edges.begin(), edges.end());
        for (auto const ic : indexCombinations)
        {
            vint error;
            for (auto i : ic) error.push_back(qs[i]);
            // Find boundary of the error
            vint boundary;
            for (auto f : error)
            {
                boundary.insert(boundary.end(), faceToEdges[f].begin(), faceToEdges[f].end());
            }
            std::sort(boundary.begin(), boundary.end());
            vint boundaryR; // restricted to v
            for (int i = 0; i < boundary.size(); ++i)
            {
                if (std::find(edges.begin(), edges.end(), boundary[i]) != edges.end())
                {
                    if (boundary[i + 1] == boundary[i]) 
                    {
                        // Edge is not in boundary mod 2 so skip over it 
                        ++i;
                    }
                    else
                    {
                        boundaryR.push_back(boundary[i]);
                    }
                }
            }
            // We only add if the boundary is not there already, as we don't care about errors that are the same up to stabilizers
            if (lift.find(boundaryR) == lift.end())
            {
                lift[boundaryR] = error;
            } 
        }
    }
    return lift;
}

vint localLift(int v, int L, const vint &paths, const vvint &vertexToEdges, const std::map<vint, vint> &lift)
{
    // Calculate lift at a given red vertex
    vint localEdges;
    auto edges = vertexToEdges[v];
    std::sort(edges.begin(), edges.end());
    for (auto const e : edges)
    {
        if (std::find(paths.begin(), paths.end(), e) != paths.end())
        {
            localEdges.push_back(e);
        }
    }
    return lift.at(localEdges);
}

void generateError(vint &qubits, sint &qubitIndices, double p, std::mt19937 &engine, std::uniform_real_distribution<double> &dist)
{
    std::fill(qubits.begin(), qubits.end(), 0);
    for (auto const q : qubitIndices)
    {
        if (dist(engine) < p)
        {
            qubits[q] = (qubits[q] + 1) % 2;
        }
    }
}

vint findCorrection(const vint &excitations, const vpint &edgeToVertices, const vvint &vertexToEdges, int L, const std::map<vint, vint> &lift)
{
    vint rgPaths, rbPaths, rgVertices, rbVertices;
    vint correction;
    // Compute matchings using pairExcitations 
    pairExcitations(excitations, edgeToVertices, vertexToEdges, b, L, rgPaths, rgVertices);
    pairExcitations(excitations, edgeToVertices, vertexToEdges, g, L, rbPaths, rbVertices);
    // Combine paths and red vertex vectors
    vint paths = rgPaths;
    paths.insert(paths.end(), rbPaths.begin(), rbPaths.end());
    vint rVerticesTmp = rgVertices;
    rVerticesTmp.insert(rVerticesTmp.end(), rbVertices.begin(), rbVertices.end());
    sint rVertices(rVerticesTmp.begin(), rVerticesTmp.end());
    // Add each g-g and b-b edge to correction
    for (auto const e : paths)
    {
        if (e >= 3 * L * L) correction.push_back(e); // g-g or b-b edges have indices larger than original max edge index of 3*L*L-1 
    }
    // For each red vertex use localLift to find qubits to be added to correction
    for (auto const rv : rVertices)
    {
        auto localCorr = localLift(rv, L, paths, vertexToEdges, lift);
        for (auto corr : localCorr) correction.push_back(corr);
    }
    return correction;
}