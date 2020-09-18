#ifndef DECODER_H
#define DECODER_H

#include "djikstra.h"
#include "hexLattice.h"
#include "PerfectMatching.h"

void calcSyndrome(vint &syndrome, const vvint &vertexToQubits, const vint &qubits, const sint &undencodedVertices);

#endif