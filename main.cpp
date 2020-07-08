#include <iostream>
#include "djikstra.h"

int main(int argc, char* argv[])
{
    std::cout << "Helo Byd!" << std::endl;

    // create the graph given in above figure
	int V = 9;
	Graph g(V);

	// making above shown graph
	g.addEdge(0, 1, 4);
	g.addEdge(0, 7, 8);
	g.addEdge(1, 2, 8);
	g.addEdge(1, 7, 11);
	g.addEdge(2, 3, 7);
	g.addEdge(2, 8, 2);
	g.addEdge(2, 5, 4);
	g.addEdge(3, 4, 9);
	g.addEdge(3, 5, 14);
	g.addEdge(4, 5, 10);
	g.addEdge(5, 6, 2);
	g.addEdge(6, 7, 1);
	g.addEdge(6, 8, 6);
	g.addEdge(7, 8, 7);

    int source = 0;
    int destination = 7;
    std::cout << "Shortest path from " << source << " to " << destination << ":\n";
    std::vector<int> path = g.shortestPath(source, destination);
    std::cout << "[";
    for (auto v : path)
    {
        std::cout << v << ",";
    }
    std::cout << "]" << std::endl;
	return 0;
}