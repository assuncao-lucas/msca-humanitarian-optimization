#include <limits>
#include <iostream>
#include "dinic.h"
#include "general.h"

void Dinic::addEdge(int u, int v, int cap)
{
    if(u == v) return;
    Edge e(v, cap, int(g[v].size()));
    Edge r(u, 0, int(g[u].size()));
    g[u].push_back(e);
    g[v].push_back(r);
}

void Dinic::reset()
{
    for(size_t u = 0; u < g.size(); u++)
    {
        for(auto& e : g[u])
        {
            e.cap = e.init_cap;
        }
    }
}

bool Dinic::buildLevelGraph(int src, int sink)
{
    fill(level.begin(), level.end(), -1);
    while(not q.empty()) q.pop();
    level[src] = 0;
    q.push(src);
    while(not q.empty())
    {
        int u = q.front();
        q.pop();
        for(auto& e : g[u])
        {
            if(not e.cap or level[e.v] != -1) continue;
            level[e.v] = level[u] + 1;
            if(e.v == sink) return true;
            q.push(e.v);
        }
    }
    return false;
}

int Dinic::blockingFlow(int u, int sink, int f)
{
    if(u == sink or not f) return f;
    int fu = f;
    for(auto& e : g[u])
    {
        if(not e.cap or level[e.v] != level[u] + 1) continue;
        int mincap = blockingFlow(e.v, sink, std::min(fu, e.cap));
        if(mincap)
        {
            g[e.v][e.rev].cap += mincap;
            e.cap -= mincap;
            fu -= mincap;
        }
    }
    if(f == fu) level[u] = -1;
    return f - fu;
}

int Dinic::maxFlow(int src, int sink)
{
    flow = 0;
    while(buildLevelGraph(src, sink))
    {
        flow += blockingFlow(src, sink, std::numeric_limits<int>::max());
        //std::cout << "* " << flow << std::endl; getchar();getchar();
    }
    return flow;
}
