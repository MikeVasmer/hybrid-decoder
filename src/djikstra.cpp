#include "djikstra.h"

// Allocates memory for adjacency list
Graph::Graph(int V)
{
	this->V = V;
	adj = new list<iPair> [V];
}

void Graph::addEdge(int u, int v, int w)
{
	adj[u].push_back(make_pair(v, w));
	adj[v].push_back(make_pair(u, w));
}

// Returns shortest path from src to dest
vector<int> Graph::shortestPath(int src, int dest)
{
    if (src == dest) return {};
	// Create a priority queue to store vertices that
	// are being preprocessed. This is weird syntax in C++.
	// Refer below link for details of this syntax
	// https://www.geeksforgeeks.org/implement-min-heap-using-stl/
	priority_queue< iPair, vector <iPair> , greater<iPair> > pq;

	// Create a vector for distances and initialize all
	// distances as infinite (INF)
	vector<int> dist(V, INF);

    // Create vector for previous vertices and initialize all to -1
    vector<int> prev(V, -1);

	// Insert source itself in priority queue and initialize
	// its distance as 0.
	pq.push(make_pair(0, src));
	dist[src] = 0;
	
	vector<bool> flags(V, false);

	/* Looping till priority queue becomes empty (or all
	distances are not finalized) */
	while (!pq.empty())
	{
		// The first vertex in pair is the minimum distance
		// vertex, extract it from priority queue.
		// vertex label is stored in second of pair (it
		// has to be done this way to keep the vertices
		// sorted distance (distance must be first item
		// in pair)
		int u = pq.top().second;
        if (u == dest) break;
		pq.pop();
		flags[u] = true;

		// 'i' is used to get all adjacent vertices of a vertex
		list< pair<int, int> >::iterator i;
		for (i = adj[u].begin(); i != adj[u].end(); ++i)
		{
			// Get vertex label and weight of current adjacent
			// of u.
			int v = (*i).first;
			int weight = (*i).second;

			// If there is shorter path to v through u.
			if (flags[v] == false && dist[v] > dist[u] + weight)
			{
				// Updating distance of v
				dist[v] = dist[u] + weight;
                prev[v] = u;
				pq.push(make_pair(dist[v], v));
			}
		}
	}

	// Reconstruct path
    vector<int> path;
    int u = dest;
    if (prev[dest] == -1)
    {
        throw std::invalid_argument("Source and destination nodes are not connected!");
    } 
    else
    {
        while (u != src)
        {
            path.push_back(u);
            u = prev[u];
        }
    }
    path.push_back(src);
    return path;
}