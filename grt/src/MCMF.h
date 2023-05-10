#ifndef MCMF_H 
#define MCMF_H

#include <vector>

template <typename weight_t>
class MCMF
{
public:
    MCMF(const int node_num);
    
    void init(const int node_num);

    /* clear all resources */
    void clear();
    
    /* add a new flow path: from u to v	*/
    void addPath(const int u, const int v, const weight_t c, const weight_t cap);
    
    /************************ Minimum Cost Maximum Flow *************************
        * To find a max-flow with min cost											*
        * allow only non-negative weight!(Dijkstra version)						*
        * Parameter:																*
        * 1. source: flow source													*
        * 2. target: flow sink														*
        * 3. minc: reference paramenter, will be the min-cost after executing MCMF	*
        * return value: max-flow value												*
        * Note: can use reweighting edge and SPFA to fit in negative-edge			*
        ****************************************************************************/
    int findMCMF(const int source, const int target, weight_t &minc);
    
    /* retrieve the flow value of the edge start -> end */
    inline int getFlow(const int start, const int end) const { return flow[start][end]; }

private:
    int gnode_num;	// number of vertex in this graph

    typedef std::vector< std::vector<weight_t> > Graph;
    
    /* flow capacity, cost per unit and flow */		
    std::vector< std::vector<int> > capacity, flow;
    std::vector< std::vector<weight_t> > cost;
    /* cost from source, and will store the min-cost augmenting path in path */
    /* reweight is used for negative cost re-weighting function */
    std::vector<weight_t> d, reweight;
    std::vector<int> path;
    
    inline weight_t potential(const int u, const int v) const
    { return (reweight[(u)] - reweight[(v)]); }
    
    /* relax v through u, cost[u][v] = w */
    bool relax(const int u, const int v, const weight_t w);
    /* reweight all edge for non-negative */
    void reweighting();
    /* Use Dijkstra to find min-cost path */
    bool Dijkstra(const int source);
};

#endif