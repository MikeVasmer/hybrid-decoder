#ifndef DECODER_H
#define DECODER_H

#include "djikstra.h"
#include "hexLattice.h"
#include "PerfectMatching.h"

void calcSyndrome(vint &syndrome, const vsint &vertexToQubits, const vint &qubits, const sint &undencodedVertices, int L);

void generateError(vint &qubits, sint &qubitIndices, double p, std::mt19937 &engine, std::uniform_real_distribution<double> &dist);

vvint combinationsUpToK(int n, int k);

std::map<vint, vint> buildLift(int L, const vsint &vertexToQubits, const sint &unencodedVertices, const vvint &vertexToEdges, const sint &qubitIndices, const vvint &faceToEdges);

vint binToList(vint &bin);

vint localLift(int v, int L, const vint &paths, const vvint &vertexToEdges, const std::map<vint, vint> &lift);

void pairExcitations(const vint &excitations, const vpint &edgeToVertices, const vvint &vertexToEdges, color c, int L, vint &paths, vint &redVertices);

vint findCorrection(const vint &excitations, const vpint &edgeToVertices, const vvint &vertexToEdges, int L, const std::map<vint, vint> &lift);

#endif