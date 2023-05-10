#include <vector>
#include <cstdint>
#include <vector>
#include <fstream>
#include <ctime>

#include "ortools/graph/min_cost_flow.h"
#include "grt/GlobalRouter.h"
#include "MCMF.h"
#include "Net.h"
#include "odb/db.h"
#include "utl/Logger.h"

namespace grt {
// The main quadtree class
class Quad {
    // Hold details of the boundary of this node
    odb::Rect cur_area;
 
    // Contains details of node
    std::vector<int> cur_nets;

    // MCMF
    // std::vector<std::vector<int>> costs;
    // std::set<int> valid_bts;
 
    // Children of this tree
    Quad* topLeftTree;
    Quad* topRightTree;
    Quad* botLeftTree;
    Quad* botRightTree;
 
public:
    inline static std::vector<Net*> nets;
    inline static std::vector<odb::Rect> nets_bbox;
    inline static std::vector<bool> assigned; // record whether the bterm is assigned.
    inline static std::vector<int> btnum4nets;
    inline static std::vector<int> nets_length;
    inline static int min_wh;
    inline static odb::dbDatabase* db;
    inline static utl::Logger* logger;
    inline static time_t cost_time;
    inline static time_t mcf_time;
    
    Quad()
    {
        cur_area = odb::Rect(0, 0, 0, 0);
        topLeftTree = NULL;
        topRightTree = NULL;
        botLeftTree = NULL;
        botRightTree = NULL;
    }
    Quad(odb::Rect area)
    {
        cur_area = area;
        topLeftTree = NULL;
        topRightTree = NULL;
        botLeftTree = NULL;
        botRightTree = NULL;
    }
    Quad(odb::Rect area, std::vector<Net*> nets_3d, std::vector<odb::Rect> nets_3d_bbox, odb::dbDatabase* db_, utl::Logger* logger_, int bterm_num, int min_wh_)
    {
        cur_area = area;
        topLeftTree = NULL;
        topRightTree = NULL;
        botLeftTree = NULL;
        botRightTree = NULL;
        nets.assign(nets_3d.begin(), nets_3d.end());
        nets_bbox.assign(nets_3d_bbox.begin(), nets_3d_bbox.end());
        db = db_;
        logger = logger;
        assigned.resize(bterm_num, false);
        btnum4nets.resize(nets_3d.size(), 0);
        nets_length.resize(nets_3d.size(), 0);
        min_wh = min_wh_;
        cost_time = 0;
        mcf_time = 0;
    }
    void init(std::vector<int> nets);
    void update();
    void clear();
    void solveMatchMCFOrtool();
    void solveMatchMCMF();
    void XWL(Net* net, std::vector<size_t> valid_bts, std::vector<size_t>& assign_cost);
    void prim(Net* net, std::vector<size_t> valid_bts, std::vector<size_t>& assign_cost);
    int primHelper(std::vector<std::vector<size_t>>& graph, int root_idx, std::vector<size_t>& bt_cost);
};
}