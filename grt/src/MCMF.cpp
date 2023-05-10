#include <limits>
#include <algorithm>
#include <iostream>

#include "MCMF.h"

/* if define reweight_edge plus SPFA, this algorithm can also use for negative edge */
#undef reweight_edge

template <typename weight_t>
MCMF<weight_t>::MCMF(const int node_num) { this->init(node_num); }

template <typename weight_t>
void MCMF<weight_t>::init(const int node_num)
{
    gnode_num = node_num;
    
    capacity.resize(node_num);
    cost.resize(node_num);
    flow.resize(node_num);
    for (auto& element : capacity) { element.resize(node_num, 0); }
    for (auto& element : cost) { element.resize(node_num, 0); }
    for (auto& element : flow) { element.resize(node_num, 0); }
    
    d.resize(node_num, 0);
    path.resize(node_num, 0);
    reweight.resize(node_num, 0);
}

template <typename weight_t>
void MCMF<weight_t>::clear()
{
    capacity.clear();
    cost.clear();
    flow.clear();
    d.clear();
    path.clear();
    reweight.clear();
}

/* add a new flow path: from u to v	*
    * flow capacity of this path: cap	*
    * flow cost of this path: c		*/
template <typename weight_t>
void MCMF<weight_t>::addPath(const int u, const int v, const weight_t c, const weight_t cap)
{
    capacity[u][v] = cap;
    cost[u][v] = c;
}

/* relax v through u, cost[u][v] = w */
template <typename weight_t>
bool MCMF<weight_t>::relax(const int u, const int v, const weight_t w)
{
    if(d[u]+w < d[v]){
        d[v] = d[u] + w;
        return true;
    }
    return false;
}

/* reweight all edge for non-negative */
template <typename weight_t>
void MCMF<weight_t>::reweighting()
{
    for(int i=0; i<gnode_num; ++i){
        if(reweight[i] < std::numeric_limits<int>::max())
            reweight[i] += d[i];
    }
}

/* Use Dijkstra to find min-cost path */
template <typename weight_t>
bool MCMF<weight_t>::Dijkstra(const int source)
{
    auto compare = [&](const int x, const int y) -> bool { return d[x] > d[y]; } ;
    std::vector<int> pque; // priority queue(min-heap)

    for (auto& element : d) { element = std::numeric_limits<int>::max(); }
    for (auto& element : path) { element = -1; }
    d[source] = 0; // important! distance to source is 0
    path[source] = source;

    // push all unvisited node into heap
    for(int i=0; i<gnode_num; ++i) pque.push_back(i);
    // adjust to mean heap
    std::make_heap(pque.begin(), pque.end(), compare);
    while(!pque.empty()){
        // pop an element from heap which d[i] is minimum
        std::pop_heap(pque.begin(), pque.end(), compare);
        const int node = pque.back();
        pque.pop_back();
        // all node is infinity, break out
        if(d[node] >= std::numeric_limits<int>::max()) return false;
        for(auto iter=pque.begin(); iter!=pque.end(); ++iter){
            // relax all vertex that is unvisited(in heap) and neighbour to node
            // try residue flow
            if((flow[*iter][node] > 0) && relax(node, *iter, potential(node, *iter) - cost[*iter][node]))
                path[*iter] = node;
            // try forward flow
            if((flow[node][*iter] < capacity[node][*iter]) && relax(node, *iter, potential(node, *iter) + cost[node][*iter]))
                path[*iter] = node;
        }
        std::make_heap(pque.begin(), pque.end(), compare);
    }
    
    #ifdef reweight_edge
        // reweighting
        reweighting();
    #endif
    
    return true;
}

template <typename weight_t>
int MCMF<weight_t>::findMCMF(const int source, const int target, weight_t &minc)
{
    for (auto& element : reweight) { element = 0; }
    for (auto& element : flow) {
        for (auto& node : element) { node = 0; }
    }

    // max flow, min cost
    int maxf = 0;
    minc = 0;
    while(Dijkstra(source) && (path[target] > 0)){
        int df = std::numeric_limits<int>::max();
        // find bottle neck
        for(int j=target, i=path[j]; i!=j; i=path[(j=i)])
            df = std::min(df, (flow[j][i]? flow[j][i] : (capacity[i][j] - flow[i][j])) );
        // update
        for(int j=target, i=path[j]; i!=j; i=path[j=i]){
            if(flow[j][i]) { flow[j][i] -= df; minc -= df * cost[j][i]; }
            else { flow[i][j] += df; minc += df * cost[i][j]; }
        }
        maxf += df;
    }
    return maxf;
}

template class MCMF<int>;