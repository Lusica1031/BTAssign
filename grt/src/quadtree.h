#include <vector>
#include <cstdint>
#include <vector>
#include <fstream>
#include <ctime>

#include "ortools/graph/min_cost_flow.h"
#include "ortools/base/logging.h"
#include "ortools/linear_solver/linear_solver.h"
#include "grt/GlobalRouter.h"
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
 
    // Children of this tree
    Quad* topLeftTree;
    Quad* topRightTree;
    Quad* botLeftTree;
    Quad* botRightTree;
 
public:
    inline static int width;

    inline static std::vector<Net*> nets;
    inline static std::vector<odb::Rect> nets_bbox;
    inline static std::vector<bool> assigned; // record whether the bterm is assigned.
    inline static std::vector<int> btnum4nets; // how many bts have been assigned
    inline static std::vector<int> nets_length; // current length
    inline static int min_wh;
    inline static odb::dbDatabase* db;
    inline static utl::Logger* logger;
    inline static time_t cost_time;
    inline static time_t mcf_time;
    inline static bool update_flag; // whether it has been updated for the next solve iteration

    // ispd
    
    inline static std::vector<odb::dbBTerm *> bterms;
    inline static int min_disp;
    inline static int max_disp;
    inline static int total_disp;

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
    // Initialization for the bterms assignment
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
        update_flag = false;
        width = area.xMax();
    }

    // Initialization for via legalization
    Quad(odb::Rect area, std::vector<Net*> nets_3d, std::vector<odb::Rect> nets_3d_bbox, std::vector<odb::dbBTerm *>& pins, odb::dbDatabase* db_, int min_wh_) {
        cur_area = area;
        bterms.assign(pins.begin(), pins.end());
        db = db_;

        nets.assign(nets_3d.begin(), nets_3d.end());
        nets_bbox.assign(nets_3d_bbox.begin(), nets_3d_bbox.end());
        min_wh = min_wh_;

        cost_time = 0;
        mcf_time = 0;
        min_disp = std::numeric_limits<int>::max();
        max_disp = 0;
        total_disp = 0;

        assigned.resize(std::ceil(area.xMax() * 1.0 / BT_PITCH) * std::ceil(area.yMax() * 1.0 / BT_PITCH), false);
        width = area.xMax();
    }

    void init(std::vector<int> nets);
    void update();
    void clear();
    void solveMatchMCFOrtool(bool bot_up = true);
    void legalize();
    void AssignMip();
    std::set<int> findDuplicates(std::vector<int> vec);
    void XWL(Net* net, std::vector<int> valid_bts, std::vector<int>& assign_cost);
    void prim(Net* net, std::vector<int> valid_bts, std::vector<int>& assign_cost);
    void primMip(Net* net, std::vector<int> valid_bts, std::vector<std::vector<int>>& assign_cost);
    int primHelper(std::vector<std::vector<int>>& graph, int root_idx, std::vector<int>& bt_cost);
    int primHelper(std::vector<std::vector<int>>& graph, int root_idx, std::vector<int>& bt_cost_i, std::vector<int>& bt_cost_j, int bt_dist);
};
}