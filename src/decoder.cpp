#include "decoder.h"

void calcSyndrome(vint &syndrome, const vvint &vertexToQubits, const vint &qubits, const sint &undencodedVertices)
{
    for (int i = 0; i < syndrome.size(); ++i)
    {
        if (undencodedVertices.find(i) != undencodedVertices.end()) continue;
        auto qubits = vertexToQubits[i];
        int parity = 0;
        for (auto q : qubits)
        {
            parity = (parity + q) % 2;
        }
        syndrome[i] = parity;
    }
}

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
    pm = new PerfectMatching(2 * nodeNum, edgeNum);
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
        }
        // We need to know the vertices in the matching
        // We need to know the paths
        // What happens to edges in 2 paths?
        // NB Edges is rb are not in rg and vice versa
    }
    delete pm;
}

vint findCorrection()
{
    // Compute matchings using pairExcitations 
    // Add each g-g and b-b edge to correction
    // For each red vertex use localLift to find qubits to be added to correction
    // Fin
}

vint localLift()
{
    // Calculate lift at a given red vertex
}