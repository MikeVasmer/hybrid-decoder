#include <iostream>
#include "djikstra.h"
#include "hexLattice.h"
#include "decoder.h"
#include <fstream>
#include <ctime>
#include <chrono>


// Returns true (1) if EC successful, false (0) if not
bool oneRun(vint &qubits, const sint &qubitIndices, double errorP, std::mt19937 &engine, std::uniform_real_distribution<double> &dist, const vsint &vertexToQubits, const vpint &edgeToVertices, const vvint &vertexToEdges, const std::map<vint, vint> &lift, vsint &logicals, const sint &unencodedVertices, int L, std::map<vertex_descriptor, std::vector<vertex_descriptor>>& excitationToPathsRG, std::map<vertex_descriptor, vint>& excitationToDistancesRG, std::map<vertex_descriptor, std::vector<vertex_descriptor>>& excitationToPathsRB, std::map<vertex_descriptor, vint>& excitationToDistancesRB)
{
    generateError(qubits, qubitIndices, errorP, engine, dist);
    vint excitations;
    calcSyndrome(excitations, vertexToQubits, qubits, unencodedVertices, L);
    auto correction = findCorrection(excitations, edgeToVertices, vertexToEdges, L, lift, unencodedVertices, excitationToPathsRG, excitationToDistancesRG, excitationToPathsRB, excitationToDistancesRB);
    for (auto q : correction)
    {
        qubits[q] = (qubits[q] + 1) % 2;
    }
    calcSyndrome(excitations, vertexToQubits, qubits, unencodedVertices, L);
    if (excitations.size() != 0) throw std::length_error("Syndrome should be trivial after correction");
    return checkCommutation(qubits, logicals);
}

void saveLattice(vpint edgeToVertices, vvint vertexToEdges, std::string fileName)
{
    std::string path = "data/lattices/" + fileName;
    // std::cout << path << std::endl;
    std::ofstream file(path);
    // Save adjacency matrix in sparse form
    for (int v = 0; v < vertexToEdges.size(); ++ v)
    {
        vint neighbours;
        for (auto edge : vertexToEdges[v])
        {
            auto vertices = edgeToVertices[edge];
            if (vertices.first == v) neighbours.push_back(vertices.second);
            else neighbours.push_back(vertices.first);
        }
        for (auto nv : neighbours)
        {
            file << nv << " ";
            // std::cout << nv << " ";
        }
        file << std::endl;
        // std::cout << std::endl;
    }
    file.close();
}

int main(int argc, char* argv[])
{
    // Get run parameters
    int L = atoi(argv[1]);
    double unencodingP = atof(argv[2]);
    double errorP = atof(argv[3]);
    bool randomizeUnencoding;
    std::stringstream ss(argv[4]);
    if (!(ss >> std::boolalpha >> randomizeUnencoding)) throw std::invalid_argument("Problem with randomizeUnencoding.");
    int trials = atoi(argv[5]);
    int job = atoi(argv[6]);
    bool rOnly = false;
    // bool rOnly = true;

    // Parameter checks
    if (L % 3 != 0 || L <= 0) throw std::invalid_argument("L must be a positive multiple of three.");
    if (unencodingP > 1 || unencodingP < 0) throw std::invalid_argument("Invalid unencoding probability.");
    if (errorP > 1 || errorP < 0) throw std::invalid_argument("Invalid error probability.");
    if (trials < 0) throw std::invalid_argument("Number of trials must be positive.");

    auto start = std::chrono::high_resolution_clock::now(); 

    // Init rng
    std::random_device rd{}; 
    std::mt19937 engine{rd()};
    // std::mt19937 engine{0}; // Deterministic seed for testing
    std::uniform_real_distribution<double> dist{0.0, 1.0};
    std::uniform_int_distribution<long long> bigDice{0, LLONG_MAX}; 

    // Build lattice
    auto vertexToQubits = buildVertexToFaces(L);
    auto edgeToVertices = buildEdgeToVertices(L);
    auto vertexToEdges = buildVertexToEdges(L);
    auto edgeToFaces = buildEdgeToFaces(L);
    auto logicals = buildLogicals(L);
    // Build Lift
    auto faceToEdges = buildFaceToEdges(L);
    auto lift = buildLift(L, vertexToQubits, vertexToEdges, faceToEdges);
    // Unencode
    sint unencodedVertices, qubitIndices;
    vint qubits;
    unencode(vertexToQubits, edgeToVertices, unencodedVertices, qubitIndices, qubits, logicals, edgeToFaces, vertexToEdges, L, unencodingP, engine, dist, randomizeUnencoding, rOnly, lift);
    // Compute paths and distances
    graph_t gr = buildGraph(edgeToVertices, b, L, edgeToFaces, qubitIndices);
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
    gr = buildGraph(edgeToVertices, g, L, edgeToFaces, qubitIndices);
    std::map<vertex_descriptor, std::vector<vertex_descriptor>> excitationToPathsRB;
    std::map<vertex_descriptor, vint> excitationToDistancesRB;
    vint verticesRB;
    for (auto const v : vertices)
    {
        if (vertexColor(v, L) != g) verticesRB.push_back(v);
    }
    shortestPaths(gr, excitationToPathsRB, excitationToDistancesRB, verticesRB);
    
    // Save lattice
    // Random number as a file name
    auto rand = bigDice(engine);
    std::time_t now = time(0);
    std::tm *ltm = std::localtime(&now);
    std::string date = std::to_string(1 + ltm->tm_mday) + "_" 
        + std::to_string(1 + ltm->tm_mon) 
        + "_" + std::to_string(1900 + ltm->tm_year);
    std::string fileName = date + "_" 
        + "L=" + std::to_string(L) + "_"
        + std::to_string(rand) + ".txt";
    saveLattice(edgeToVertices, vertexToEdges, fileName);

    // Monte Carlo
    int fails = 0;
    int t = 0;
    while (t < trials)
    {
        fails += !oneRun(qubits, qubitIndices, errorP, engine, dist, vertexToQubits, edgeToVertices, vertexToEdges, lift, logicals, unencodedVertices, L, excitationToPathsRG, excitationToDistancesRG, excitationToPathsRB, excitationToDistancesRB);
        ++t;
        double pfail = fails / t;
        double err = sqrt(pfail * (1-pfail) / t);
        if (err != 0 && err / pfail < 0.01) break;
    }

    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start); 

    // Save results
    std::string path = "data/L=" 
        + std::to_string(L) + "_"
        + "uP=" + std::to_string(unencodingP) + "_"
        + "eP=" + std::to_string(errorP) + "_"
        + "job=" + std::to_string(job) + "_"
        + "date=" + date + "_"
        + "lat=" + std::to_string(rand) + ".csv";
    std::ofstream file(path);
    // L, error p, unencoding p, trials, fails, job
    file << "L,uP,eP,randomU,trials,fails,runtime(s),job" << std::endl;
    file << std::to_string(L) << ","
        << std::to_string(unencodingP) << ","
        << std::to_string(errorP) << ","
        << std::to_string(randomizeUnencoding) << ","
        << std::to_string(t) << ","
        << std::to_string(fails) << ","
        << std::to_string(duration.count()) << ","
        << std::to_string(job) << std::endl;

	return 0;
}
