#ifndef DINIC_H
#define DINIC_H

#include <algorithm>
#include <vector>
#include <queue>

class Edge
{
public:
    Edge(int v_, int cap_, int rev_): v(v_), rev(rev_), cap(cap_), init_cap(cap_) {}
    ~Edge(){};
    int v;
    int rev;
	int cap;
	int init_cap;
};

class Dinic
{
public:
    Dinic(int n_): g(n_), level(n_), n(n_) {}
    ~Dinic(){}
    std::vector<std::vector<Edge>> g;
	std::vector<int> level;
	std::queue<int> q;
	int flow;
	int n;
	void reset();
    void addEdge(int u, int v, int cap);
	bool buildLevelGraph(int src, int sink);
	int blockingFlow(int u, int sink, int f);
	int maxFlow(int src, int sink);
};
#endif
