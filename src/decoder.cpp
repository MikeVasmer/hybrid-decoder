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

vint binToList(vint &bin)
{
    vint list;
    for (int i = 0; i < bin.size(); ++i)
    {
        if (bin[i] == 1) list.push_back(i);
    }
    return list;
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

void generateError(vint &qubits, const sint &qubitIndices, double p, std::mt19937 &engine, std::uniform_real_distribution<double> &dist)
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

// c is color of unencoded face
vint reRoute(const int v1, const int v2, const int centralV, const int c, const int L, const vvint &vertexToEdges, const vpint &edgeToVertices)
{
    // (v1, v2) is an r-edge
    auto v1Neighbors = ccNeighbors(v1, L);
    // Find a an alternate path v1 -> v2 via a {g, b} \ c edge
    vint newPath;
    int w1, w2;
    for (auto const v : v1Neighbors)
    {
        if (vertexColor(v, L) == c) continue;
        auto vNeighbors = ccNeighbors(v, L);
        if (vNeighbors.find(centralV) != vNeighbors.end())
        {
            w1 = v;
            newPath.push_back(pairToEdge(v1, w1, vertexToEdges, edgeToVertices));
            break;
        } 
    }
    auto w1Edges = vertexToEdges[w1];
    for (auto const e : w1Edges)
    {
        auto vs = edgeToVertices[e];
        if (vertexColor(vs.first, L) != vertexColor(vs.second, L)) continue;
        if (vs.first == w1) w2 = vs.second;
        else w2 = vs.first;
        auto w2Neighbors = ccNeighbors(w2, L);
        if (w2Neighbors.find(centralV) != w2Neighbors.end() && w2Neighbors.find(v2) != w2Neighbors.end()) 
        {
            newPath.push_back(e);
            newPath.push_back(pairToEdge(v2, w2, vertexToEdges, edgeToVertices));
            break; 
        }
    }
    if (newPath.size() == 1) // Sometimes it goes the wrong way
    {
        newPath.clear();
        for (auto const v : v1Neighbors)
        {
            if (vertexColor(v, L) == c) continue;
            auto vNeighbors = ccNeighbors(v, L);
            // Enforce other choice for w1
            if (vNeighbors.find(centralV) != vNeighbors.end() && v != w1)
            {
                w1 = v;
                newPath.push_back(pairToEdge(v1, w1, vertexToEdges, edgeToVertices));
                break;
            } 
        }
        w1Edges = vertexToEdges[w1];
        for (auto const e : w1Edges)
        {
            auto vs = edgeToVertices[e];
            if (vertexColor(vs.first, L) != vertexColor(vs.second, L)) continue;
            if (vs.first == w1) w2 = vs.second;
            else w2 = vs.first;
            auto w2Neighbors = ccNeighbors(w2, L);
            if (w2Neighbors.find(centralV) != w2Neighbors.end() && w2Neighbors.find(v2) != w2Neighbors.end()) 
            {
                newPath.push_back(e);
                newPath.push_back(pairToEdge(v2, w2, vertexToEdges, edgeToVertices));
                break; 
            }
        }
    }
    return newPath;
}

sint findCorrection(const vint &excitations, const vpint &edgeToVertices, const vvint &vertexToEdges, int L, const std::map<vint, vint> &lift, const sint &unencodedVertices)
{
    vint rgPaths, rbPaths, rgVertices, rbVertices;
    sint correction;
    // Compute matchings using pairExcitations 
    pairExcitations(excitations, edgeToVertices, vertexToEdges, b, L, rgPaths, rgVertices);
    pairExcitations(excitations, edgeToVertices, vertexToEdges, g, L, rbPaths, rbVertices);
    // Process matchings
    vvint pathTup = {rgPaths, rbPaths};
    vint cols = {g, b};
    vint changeEdges;
    for (int i = 0; i < 2; ++i)
    {
        for (auto const e : pathTup[i])
        {
            if (e >= 3 * L * L) // toric code edges have indices larger than original max edge index of 3*L*L-1
            {
                if (vertexColor(edgeToVertices[e].first, L) != r)
                {
                    // Add each g-g and b-b edge to correction
                    if (correction.find(e) != correction.end()) correction.erase(e);
                    else correction.insert(e); 
                } 
                else {
                    int v1 = edgeToVertices[e].first;
                    int v2 = edgeToVertices[e].second;
                    auto v1Neighbors = ccNeighbors(v1, L); // Color code neighbors
                    auto v2Neighbors = ccNeighbors(v2, L);
                    // Find the vertex of the unencoded face
                    int cV, c;
                    for (auto const v : v1Neighbors)
                    {
                        if (v2Neighbors.find(v) != v2Neighbors.end())
                        {
                            if (unencodedVertices.find(v) != unencodedVertices.end())
                            {
                                cV = v;
                                break;
                            }
                        }
                    }
                    c = vertexColor(cV, L); // color 
                    // If r-edge in rgParing is from a g-face then add to corr, analogous for b case
                    if (c == cols[i]) correction.insert(e); 
                    // Get the new path connecting v1 and v2 and add to the edges to be changed
                    auto newEdges = reRoute(v1, v2, cV, cols[i], L, vertexToEdges, edgeToVertices);
                    for (auto const ne : newEdges)
                    {
                        changeEdges.push_back(ne);
                        // Add any re-routed g-g or b-b edges to the correction
                        if (ne >= 3 * L * L)
                        {
                            // We add modulo 2 (these edges could have been added earlier)
                            if (correction.find(ne) != correction.end()) correction.erase(ne);
                            else correction.insert(ne);
                        } 
                    } 
                    changeEdges.push_back(e); // We want to remove e
                }
            }    
        }
    }
    // Combine paths and red vertex vectors
    vint paths = rgPaths;
    paths.insert(paths.end(), rbPaths.begin(), rbPaths.end());
    vint rVerticesTmp = rgVertices;
    rVerticesTmp.insert(rVerticesTmp.end(), rbVertices.begin(), rbVertices.end());
    sint rVertices(rVerticesTmp.begin(), rVerticesTmp.end());
    // Change the pairing according to edges computed earlier
    for (auto const e : changeEdges)
    {
        if (std::find(paths.begin(), paths.end(), e) == paths.end())
        {
            paths.push_back(e);
        }
        else
        {
            paths.erase(std::find(paths.begin(), paths.end(), e));
        }
    }
    // For each red vertex use localLift to find qubits to be added to correction
    for (auto const rv : rVertices)
    {
        auto localCorr = localLift(rv, L, paths, vertexToEdges, lift);
        for (auto corr : localCorr)
        {
            // We add modulo 2
            if (correction.find(corr) != correction.end()) correction.erase(corr);
            else correction.insert(corr);
        }
    }
    return correction;
}

bool checkCommutation(vint &qubits, vsint &logicals)
{
    for (auto &logical : logicals)
    {
        int parity = 0;
        for (auto q : logical)
        {
            if (qubits[q] == 1) ++parity;
        }
        if (parity % 2 == 1) return false;
    }
    return true;
}