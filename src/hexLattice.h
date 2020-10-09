#ifndef HEXLATTICE_H
#define HEXLATTICE_H

#include <vector>
#include <set>
#include <random>
#include <map>

typedef std::vector<int> vint;
typedef std::vector<std::vector<int>> vvint;
typedef std::set<int> sint;
typedef std::vector<std::set<int>> vsint;
typedef std::vector<std::pair<int, int>> vpint;

struct coord
{
    int xi[2];

    bool operator==(const coord& c) const
    {
        return (xi[0] == c.xi[0] && xi[1] == c.xi[1]);
    }
};

enum direction {x, y, xy};

enum color {r, g, b};

int coordToIndex(const coord &c, int L);

coord indexToCoord(int i, int L);

int neigh(int v, int dir, int sign, int L);

int edgeIndex(int v, int dir, int sign, int L);

int faceIndex(int v, int dir);

vvint buildVertexToEdges(int L);

vvint buildVertexToEdges(const vpint &edgeToVertices, int L);

vsint buildVertexToFaces(int L);

vvint buildFaceToEdges(int L);

int vertexColor(int v, int L);

vpint buildEdgeToFaces(int L);

vpint buildEdgeToVertices(int L);

vsint buildLogicals(int L);

void localUnencoding(vvint &vertexToEdgesLocal, vpint &edgeToVerticesLocal, vint &localVertexMap, std::map<std::pair<int, int>, std::vector<int>> &logicalOperatorMap, int v, int L);

void modifyVertexToQubits(vsint &vertexToQubits, const vvint &vertexToEdgesL, const vint &localVertexMap, int v, int e, int L);

void modifyEdgeToVertices(vpint &edgeToVerticesU, const vpint &edgeToVerticesL, int v, int L, const vint &localVerticesMap);

void modifyLogicalOperators(vsint &logicals, const std::map<std::pair<int, int>, std::vector<int>> &logicalOperatorMap, int e);

void removeRedundantEdges(vpint &edgeToVertices, const sint &qubitIndices, const vpint &edgeToFaces, int L);

void unencode(vsint &vertexToQubits, vpint &edgeToVertices, sint &unencodedVertices, sint &qubitIndices, vint &qubits, vsint& logicals, const vpint &edgeToFaces, vvint &vertexToEdges, int L, double p, std::mt19937& engine, std::uniform_real_distribution<double>& dist);

int pairToEdge(int v, int u, const vvint &vertexToEdges, const vpint &edgeToVertices);

#endif