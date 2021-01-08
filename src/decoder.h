#ifndef DECODER_H
#define DECODER_H

#include "djikstra.h"
#include "hexLattice.h"
#include "PerfectMatching.h"

void calcSyndrome(vint &syndrome, const vsint &vertexToQubits, const vint &qubits, const sint &undencodedVertices, int L);

void generateError(vint &qubits, const sint &qubitIndices, double p, std::mt19937 &engine, std::uniform_real_distribution<double> &dist);

vint binToList(vint &bin);

vint localLift(int v, int L, const vint &paths, const vvint &vertexToEdges, const std::map<vint, vint> &lift);

void pairExcitations(const vint &excitations, const vpint &edgeToVertices, const vvint &vertexToEdges, color c, int L, vint &paths, vint &redVertices);

sint findCorrection(const vint &excitations, const vpint &edgeToVertices, const vvint &vertexToEdges, int L, const std::map<vint, vint> &lift, const sint &unencodedVertices);

bool checkCommutation(vint &qubits, vsint &logicals);

vint reRoute(const int v1, const int v2, const int centralV, const int c, const int L, const vvint &vertexToEdges, const vpint &edgeToVertices);

#endif