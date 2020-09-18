#include "hexLattice.h"
#include <math.h>
#include <algorithm>
#include <iostream>

int coordToIndex(const coord &c, int L)
{
    return c.xi[0] + c.xi[1] * L;
}

coord indexToCoord(int i, int L)
{
    coord c;
    c.xi[0] = i % L;
    c.xi[1] = (int)floor(i / L) % L;
    return c;
}

int neigh(int v, int dir, int sign, int L)
{
    coord c = indexToCoord(v, L);
    if (dir == x)
    {
        c.xi[0] = (L + ((c.xi[0] + sign) % L)) % L;
    }
    else if (dir == y)
    {
        c.xi[1] = (L + ((c.xi[1] + sign) % L)) % L;
    }
    else if (dir == xy)
    {
        c.xi[0] = (L + ((c.xi[0] + sign) % L)) % L;
        c.xi[1] = (L + ((c.xi[1] + sign) % L)) % L;
    }
    return coordToIndex(c, L);
}

int edgeIndex(int v, int dir, int sign, int L)
{
    if (sign < 0)
    {
        v = neigh(v, dir, sign, L);
    }
    return 3 * v + dir;
}

int faceIndex(int v, int dir)
{
    return 2 * v + dir;
}

vvint buildVertexToEdges(int L)
{
    vvint vertexToEdges;
    for (int v = 0; v < L * L; ++v)
    {
        vint edges;
        edges.push_back(edgeIndex(v, x, 1, L));
        edges.push_back(edgeIndex(v, y, 1, L));
        edges.push_back(edgeIndex(v, xy, 1, L));
        edges.push_back(edgeIndex(v, x, -1, L));
        edges.push_back(edgeIndex(v, y, -1, L));
        edges.push_back(edgeIndex(v, xy, -1, L));
        std::sort(edges.begin(), edges.end());
        vertexToEdges.push_back(edges);
    }
    return vertexToEdges;
}

vsint buildVertexToFaces(int L)
{
    vsint vertexToFaces;
    for (int v = 0; v < L * L; ++v)
    {
        std::set<int> faces;
        faces.insert(faceIndex(v, x));
        faces.insert(faceIndex(v, y));
        faces.insert(faceIndex(neigh(v, x, -1, L), x));
        faces.insert(faceIndex(neigh(v, y, -1, L), y));
        faces.insert(faceIndex(neigh(v, xy, -1, L), x));
        faces.insert(faceIndex(neigh(v, xy, -1, L), y));
        vertexToFaces.push_back(faces);
    }
    return vertexToFaces;
}

vvint buildFaceToEdges(int L)
{
    vvint faceToEdges;
    for (int f = 0; f < 2 * L * L; ++f)
    {
        vint edges;
        int v = f / 2;
        int dir = f % 2;
        edges.push_back(edgeIndex(v, xy, 1, L));
        if (dir == x)
        {
            edges.push_back(edgeIndex(v, x, 1, L));
            edges.push_back(edgeIndex(neigh(v, x, 1, L), y, 1, L));
        }
        else if (dir == y)
        {
            edges.push_back(edgeIndex(v, y, 1, L));
            edges.push_back(edgeIndex(neigh(v, y, 1, L), x, 1, L));
        }
        std::sort(edges.begin(), edges.end());
        faceToEdges.push_back(edges);
    }
    return faceToEdges;
} 

vpint buildEdgeToFaces(int L)
{
    vpint edgeToFaces;
    for (int e = 0; e < 3 * L * L; ++e)
    {
        int v = e / 3;
        int dir = e % 3;
        int f1, f2;
        if (dir == x)
        {
            f1 = faceIndex(v, x);
            f2 = faceIndex(neigh(v, y, -1, L), y);
        }
        else if (dir == y)
        {
            f1 = faceIndex(v, y);
            f2 = faceIndex(neigh(v, x, -1, L), x);
        }
        else if (dir == xy)
        {
            f1 = faceIndex(v, x);
            f2 = faceIndex(v, y);
        }
        if (f1 < f2)
        {
            edgeToFaces.push_back({f1, f2});
        }
        else if (f2 < f1)
        {
            edgeToFaces.push_back({f2, f1});
        }
    }
    return edgeToFaces;
}

int vertexColor(int v, int L)
{
    coord c = indexToCoord(v, L);
    return (c.xi[0] + c.xi[1]) % 3;
}

vpint buildEdgeToVertices(int L)
{
    vpint edgeToVertices;
    for (int e = 0; e < 3 * L * L; ++e)
    {
        int v1 = e / 3;
        int dir = e % 3;
        int v2 = neigh(v1, dir, 1, L);
        if (v1 < v2)
        {
            edgeToVertices.push_back({v1, v2});
        }
        else if (v2 < v1)
        {
            edgeToVertices.push_back({v2, v1});
        }
    }
    return edgeToVertices;
}

vsint buildLogicals(int L)
{
    vsint logicals;
    vint startingVs = {1, 2}; // (1, 0) green (2, 0) blue
    for (auto const v : startingVs)
    {
        sint logical;
        int u = v;
        for (int i = 0; i < L; ++i)
        {
            logical.insert(faceIndex(neigh(u, x, -1, L), x));
            logical.insert(faceIndex(neigh(u, x, -1, L), y));
            u = neigh(neigh(u, x, -1, L), y, 1, L);
        }
        logicals.push_back(logical);
        logical.clear();
        u = v;
        for (int i = 0; i < L; ++i)
        {
            logical.insert(faceIndex(u, x));
            logical.insert(faceIndex(neigh(u, x, 1, L), y));
            u = neigh(neigh(neigh(u, x, 1, L), x, 1, L), y, 1, L);            
        }
        logicals.push_back(logical);
    }
    return logicals;
}

void localUnencoding(vvint &vertexToEdgesLocal, vpint &edgeToVerticesLocal, vint &localVertexMap, std::map<std::pair<int, int>, std::vector<int>> &logicalOperatorMap, int v, int L)
{
    // ToDo: Randomize the unencoding
    // Make below two are sorted
    vertexToEdgesLocal = {{0, 1}, {2}, {3}, {2, 3}, {0}, {1}};
    edgeToVerticesLocal = {{0, 4}, {0, 5}, {1, 3}, {2, 3}};
    localVertexMap = {
        neigh(v, xy, -1, L), 
        neigh(v, y, -1, L), 
        neigh(v, x, -1, L), 
        neigh(v, xy, 1, L),
        neigh(v, x, 1, L),
        neigh(v, y, 1, L)
    };
    logicalOperatorMap = {
        {{faceIndex(localVertexMap[0], x), faceIndex(localVertexMap[1], y)}, {0}},
        {{faceIndex(v, x), faceIndex(v, y)}, {0, 1}},
        {{faceIndex(localVertexMap[2], x), faceIndex(v, y)}, {3}},
        {{faceIndex(localVertexMap[0], x), faceIndex(localVertexMap[0], y)}, {2, 3}}
    };
}

void modifyVertexToQubits(vsint &vertexToQubits, const vvint &vertexToEdgesL, const vint &localVertexMap, int v, int e, int L)
{
    vertexToQubits[v] = {e, e + 1, e + 2, e + 3}; // Useful information for decoder
    int u = localVertexMap[0];
    vertexToQubits[u].erase(faceIndex(u, x));
    vertexToQubits[u].erase(faceIndex(u, y));
    for (auto eL : vertexToEdgesL[0]) vertexToQubits[u].insert(e + eL);
    u = localVertexMap[1];
    vertexToQubits[u].erase(faceIndex(u, y));
    vertexToQubits[u].erase(faceIndex(neigh(u, x, -1, L), x));
    for (auto eL : vertexToEdgesL[1]) vertexToQubits[u].insert(e + eL);
    u = localVertexMap[2];
    vertexToQubits[u].erase(faceIndex(u, x));
    vertexToQubits[u].erase(faceIndex(neigh(u, y, -1, L), y));
    for (auto eL : vertexToEdgesL[2]) vertexToQubits[u].insert(e + eL);
    u = localVertexMap[3];
    vertexToQubits[u].erase(faceIndex(v, x));
    vertexToQubits[u].erase(faceIndex(v, y));
    for (auto eL : vertexToEdgesL[3]) vertexToQubits[u].insert(e + eL);
    u = localVertexMap[4];
    vertexToQubits[u].erase(faceIndex(v, x));
    vertexToQubits[u].erase(faceIndex(neigh(v, y, -1, L), y));
    for (auto eL : vertexToEdgesL[4]) vertexToQubits[u].insert(e + eL);
    u = localVertexMap[5];
    vertexToQubits[u].erase(faceIndex(v, y));
    vertexToQubits[u].erase(faceIndex(neigh(v, x, -1, L), x));
    for (auto eL : vertexToEdgesL[5]) vertexToQubits[u].insert(e + eL);
}

void modifyEdgeToVertices(vpint &edgeToVerticesU, const vpint &edgeToVerticesL, int v, int L, const vint &localVerticesMap)
{
    // U for unencoded, L for local
    for (int i = 0; i < 4; ++i) 
    {
        int v1 = localVerticesMap[edgeToVerticesL[i].first];
        int v2 = localVerticesMap[edgeToVerticesL[i].second];
        if (v1 < v2)
        {
            edgeToVerticesU.push_back({v1, v2});
        }
        else if (v2 < v1)
        {
            edgeToVerticesU.push_back({v2, v1});
        }
    }
    // Delete edges connected to unencoded vertex NB {} -> {0, 0} (default pair)
    edgeToVerticesU[edgeIndex(v, x, 1, L)] = {};
    edgeToVerticesU[edgeIndex(v, y, 1, L)] = {};
    edgeToVerticesU[edgeIndex(v, xy, 1, L)] = {};
    edgeToVerticesU[edgeIndex(v, x, -1, L)] = {};
    edgeToVerticesU[edgeIndex(v, y, -1, L)] = {};
    edgeToVerticesU[edgeIndex(v, xy, -1, L)] = {};
}

vvint buildVertexToEdges(const vpint &edgeToVertices, int L)
{
    vvint vertexToEdges;
    vertexToEdges.assign(L * L, {});
    for (int e = 0; e < edgeToVertices.size(); ++e)
    {
        int v1 = edgeToVertices[e].first;
        int v2 = edgeToVertices[e].second;
        if (v1 == 0 && v2 == 0) continue;
        vertexToEdges[v1].push_back(e);
        vertexToEdges[v2].push_back(e);
    }
    return vertexToEdges;
}

// Get rid of edges belonging to two qubits that have been unencoded
void removeRedundantEdges(vpint &edgeToVertices, const sint &qubitIndices, const vpint &edgeToFaces, int L)
{
    for (int e = 0; e < 3 * L * L; ++e)
    {
        auto faces = edgeToFaces[e];
        if ((qubitIndices.find(faces.first) != qubitIndices.end()) || (qubitIndices.find(faces.second) != qubitIndices.end())) continue;
        edgeToVertices[e] = {};
    }
}

void modifyLogicalOperators(vsint &logicals, const std::map<std::pair<int, int>, std::vector<int>> &logicalOperatorMap, int e)
{
    for (auto &logical : logicals)
    {
        for (auto const &item : logicalOperatorMap)
        {
            auto key = item.first;
            auto value = item.second;
            if ((logical.find(key.first) != logical.end()) && (logical.find(key.second) != logical.end()))
            {
                logical.erase(key.first);
                logical.erase(key.second);
                for (auto const eL : value) logical.insert(e + eL);
            }
        }
    }
}

// vertexToQubits should be a copy of vertexToFaces
void unencode(std::map<std::pair<int, int>, int> &vertexPairToEdgeU, vsint &vertexToQubits, vpint &edgeToVertices, sint &unencodedVertices, sint &qubitIndices, vint &qubits, vsint& logicals, const vpint &edgeToFaces, vvint &vertexToEdges, int L, double p, std::mt19937& engine, std::uniform_real_distribution<double>& dist)
{
    for (int i = 0; i < 2 * L * L; ++i) qubitIndices.insert(i);
    int e = 3 * L * L;
    for (int v = 0; v < L * L; ++v)
    {
        if (vertexColor(v, L) != r) continue;
        if (dist(engine) <= p)
        {
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
            std::map<std::pair<int, int>, std::vector<int>> logicalOperatorMap;
            localUnencoding(vertexToEdgesL, edgeToVerticesL, localVertexMap, logicalOperatorMap, v, L);
            // Modify data structures to account for unencoding
            modifyVertexToQubits(vertexToQubits, vertexToEdgesL, localVertexMap, v, e, L);
            modifyEdgeToVertices(edgeToVertices, edgeToVerticesL, v, L, localVertexMap);
            modifyLogicalOperators(logicals, logicalOperatorMap, e);
            // vertexPairToEdgeU
            for (int i = 0; i < 4; ++i)
            {
                int v1 = localVertexMap[edgeToVerticesL[i].first];
                int v2 = localVertexMap[edgeToVerticesL[i].second];
                if (v1 < v2)
                {
                    vertexPairToEdgeU.insert({{v1, v2}, e + i});
                }
                else if (v2 < v1)
                {
                    vertexPairToEdgeU.insert({{v2, v1}, e + i});
                }
                else
                {
                    throw std::invalid_argument("v1 should not equal v2");
                }
                
            }
            e += 4;
        }
    }
    removeRedundantEdges(edgeToVertices, qubitIndices, edgeToFaces, L);
    qubits.assign(e, 0); // Init qubits with length equal to final edge index + 1
    vertexToEdges = buildVertexToEdges(edgeToVertices, L); // Use overloaded function to invert edgeToVertices
}

int pairToEdge(int v, int u, const vvint &vertexToEdges, const vpint &edgeToVertices)
{
    for (auto e : vertexToEdges[v])
    {
        if (edgeToVertices[e].first == u || edgeToVertices[e].second == u) return e;
    }
    throw std::invalid_argument("no edge exists between u and v");
    return -1;
}