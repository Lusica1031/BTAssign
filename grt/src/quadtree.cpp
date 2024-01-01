#include "quadtree.h"

namespace grt {

using utl::GRT;
using namespace operations_research;

// Init the quad tree
void Quad::init(std::vector<int> prev_nets_idxs) {
    // partitioning boundary
    int xl = cur_area.xMin();
    int xh = cur_area.xMax();
    int yl = cur_area.yMin();
    int yh = cur_area.yMax();

    if(prev_nets_idxs.size() == 0) {
        return;
    }

    if(((xh - xl) < min_wh) || ((yh - yl) < min_wh)) {
        cur_nets.assign(prev_nets_idxs.begin(), prev_nets_idxs.end());
        return;
    }

    int xmid = (xl + xh) / 2;
    int ymid = (yl + yh) / 2;

    
    std::vector<int> lb_nets;
    std::vector<int> rb_nets;
    std::vector<int> lt_nets;
    std::vector<int> rt_nets;

    odb::Rect lb_bbox(xl, yl, xmid, ymid);
    odb::Rect lt_bbox(xl, ymid, xmid, yh);
    odb::Rect rb_bbox(xmid, yl, xh, ymid);
    odb::Rect rt_bbox(xmid, ymid, xh, yh);

    for (auto net_id : prev_nets_idxs) {

        odb::Rect net_bb = nets_bbox[net_id];

        if(lb_bbox.contains(nets_bbox[net_id])) {
            // left bot
            lb_nets.push_back(net_id);
            continue;
        }
        else if(rb_bbox.contains(nets_bbox[net_id])) {
            // right bot
            rb_nets.push_back(net_id);
            continue;
        }
        else if(lt_bbox.contains(nets_bbox[net_id])) {
            // left top
            lt_nets.push_back(net_id);
            continue;
        } else if (rt_bbox.contains(nets_bbox[net_id])) {
            // right top
            rt_nets.push_back(net_id);
            continue;
        } else {
            // cur
            cur_nets.push_back(net_id);
            continue;
        }
    }

    // init leaf nodes
    topLeftTree = new Quad(lt_bbox);
    topRightTree = new Quad(rt_bbox);
    botLeftTree = new Quad(lb_bbox);
    botRightTree = new Quad(rb_bbox);

    topLeftTree->init(lt_nets);
    topRightTree->init(rt_nets);
    botLeftTree->init(lb_nets);
    botRightTree->init(rb_nets);
}
// update the quad tree
void Quad::update() {

    update_flag = true;

    // partitioning boundary
    for (auto net = cur_nets.begin(); net != cur_nets.end();) {
        if(nets[(*net)]->getNumPins() <= 10 || nets[(*net)]->getNumPins() >= 100)
            net = cur_nets.erase(net);
        else
            net++;
    }

    if(topLeftTree)
        topLeftTree->update();
    if(topRightTree)
        topRightTree->update();
    if(botLeftTree)
        botLeftTree->update();
    if(botRightTree)
        botRightTree->update();
}

void Quad::clear() {
    cur_nets.clear();
    if(topLeftTree)
        topLeftTree->clear();
    if(topRightTree)
        topRightTree->clear();
    if(botLeftTree)
        botLeftTree->clear();
    if(botRightTree)
        botRightTree->clear();
    
    if(!nets.empty()) {
        nets.clear();
        assigned.clear();
        btnum4nets.clear();
    }
}

void Quad::solveMatchMCFOrtool(bool bot_up) {

    // solve the sub problems
    if(bot_up) {
        if(topLeftTree) {
            topLeftTree->solveMatchMCFOrtool(bot_up);
        }
        if(botLeftTree) {
            botLeftTree->solveMatchMCFOrtool(bot_up);
        }
        if(topRightTree) {
            topRightTree->solveMatchMCFOrtool(bot_up);
        }
        if(botRightTree) {
            botRightTree->solveMatchMCFOrtool(bot_up);
        }
    }

    if(cur_nets.size() == 0) {
        if (!bot_up) {
            if(topLeftTree) {
                topLeftTree->solveMatchMCFOrtool(bot_up);
            }
            if(botLeftTree) {
                botLeftTree->solveMatchMCFOrtool(bot_up);
            }
            if(topRightTree) {
                topRightTree->solveMatchMCFOrtool(bot_up);
            }
            if(botRightTree) {
                botRightTree->solveMatchMCFOrtool(bot_up);
            }
        }
        return;
    }
    
    time_t cost_s = time(NULL);
    // get all valid bterms
    std::vector<int> valid_bts;

    int start_x = std::ceil(cur_area.xMin() * 1.0 / BT_PITCH) * BT_PITCH;
    int start_y = std::ceil(cur_area.yMin() * 1.0 / BT_PITCH) * BT_PITCH;

    int bt_col_num = std::ceil((cur_area.xMax() - start_x) * 1.0 / BT_PITCH);
    int bt_row_num = std::ceil((cur_area.yMax() - start_y) * 1.0 / BT_PITCH);
    

    int skip_col  = 1;
    if((bt_col_num * bt_row_num > 1)) {
        // std::cout << " [DEBUG] nets number / bots number: " << (bt_col_num * bt_row_num) / cur_nets.size() << std::endl;
        if(update_flag)
            // skip_col = 1;
            // skip_col = std::max(std::round((bt_col_num * bt_row_num) / (cur_nets.size() * 4.0)), 1.0);
            skip_col = std::max(std::round((bt_col_num * bt_row_num) / (cur_nets.size() * 4.0)), 1.0);
        else
            // skip_col = 2;
            // skip_col = std::max(std::round((bt_col_num * bt_row_num) / (cur_nets.size() * 8.0)), 1.0);
            skip_col = std::max(std::round((bt_col_num * bt_row_num) / (cur_nets.size() * 8.0)), 1.0);
    }

    for(int row = 0; row < bt_row_num; row++) {
        int cur_y = start_y + row * BT_PITCH;

        for(int col = 0; col < bt_col_num; col += skip_col) {
            int cur_x = start_x + col * BT_PITCH + (row % skip_col) * BT_PITCH;
            if(cur_x > cur_area.xMax())
                break;
            int cur_idx = ((cur_y / BT_PITCH)) * std::ceil(width * 1.0 / BT_PITCH) + (cur_x / BT_PITCH) ;
            if(!assigned[cur_idx]) {
                valid_bts.push_back(cur_idx);
            }
        }
    }

    if(valid_bts.size() < cur_nets.size()) {
        std::cout << "the " << (valid_bts.size()) << " bt candidates are not enough for " << cur_nets.size() << " nets." << std::endl;
        skip_col = skip_col / 2;

        valid_bts.clear();

        for(int row = 0; row < bt_row_num; row++) {
            int cur_y = start_y + row * BT_PITCH;

            for(int col = 0; col < bt_col_num; col += skip_col) {
                int cur_x = start_x + col * BT_PITCH + (row % skip_col) * BT_PITCH;
                if(cur_x > cur_area.xMax())
                    break;
                int cur_idx = ((cur_y / BT_PITCH)) * std::ceil(width * 1.0 / BT_PITCH) + (cur_x / BT_PITCH) ;
                if(!assigned[cur_idx]) {
                    valid_bts.push_back(cur_idx);
                }
            }
        }
        
        if(valid_bts.size() < cur_nets.size()) {
            std::cout << "the " << (valid_bts.size()) << " bt candidates are not enough for " << cur_nets.size() << " nets." << std::endl;
            if(update_flag)
                return;
            exit(1);
        }
    }

    bool a = std::set<int>(valid_bts.begin(),valid_bts.end()).size()!=valid_bts.size();
    if(a) {
        std::cout << " [ERROR] repeat valid bonding terminals. " << std::endl;
        std::cout << " [ERROR] valid_bts.size(): " << valid_bts.size() << " with start_x,y: " << start_x << ", " << start_y << " and xMax, yMax: " << cur_area.xMax() << ", " << cur_area.yMax() << " bt_row_num: " << bt_row_num << " bt_col_num: " << bt_col_num << " with unique: " << std::set<int>(valid_bts.begin(),valid_bts.end()).size() << std::endl;
        auto dup_bts = findDuplicates(valid_bts);
        for (auto dup_bt : dup_bts) {
            std::cout << dup_bt << std::endl;
        }
        exit(1);
    }

    // calculate costs
    std::vector<std::vector<int>> costs(cur_nets.size());

    #pragma omp parallel for
    for(int i=0; i< cur_nets.size(); i++) {
        prim(nets[cur_nets[i]], valid_bts, costs[i]);
        for(int j=0; j<valid_bts.size(); j++) {
            costs[i][j] = costs[i][j] - nets_length[cur_nets[i]];
        }
    }

    // Instantiate a SimpleMinCostFlow solver.
    SimpleMinCostFlow min_cost_flow;

    // init MCF flow
    int vertex_num = cur_nets.size() + valid_bts.size() + 1;
    int min_cost;
    
    for (int from = 0; from < cur_nets.size(); from++) {
        min_cost_flow.SetNodeSupply(from, 1);
        for (int to = cur_nets.size(); to < vertex_num - 1; to++) {
            int arc = min_cost_flow.AddArcWithCapacityAndUnitCost(from, to, 1, costs[from][to - cur_nets.size()]);
            if (arc != from * valid_bts.size() + (to-cur_nets.size())) LOG(FATAL) << "Internal error";
        }
    }
    
    for (int to = cur_nets.size(); to < vertex_num - 1; to++) {
        int arc = min_cost_flow.AddArcWithCapacityAndUnitCost(to, vertex_num-1, 1, 1);
        if (arc != valid_bts.size()*cur_nets.size() + (to - cur_nets.size())) LOG(FATAL) << "Internal error";
    }

    min_cost_flow.SetNodeSupply(vertex_num - 1, -cur_nets.size());
    
    time_t cost_e = time(NULL);

    cost_time += (cost_e - cost_s);

    // Find the min cost flow.
    int status = min_cost_flow.Solve();
    
    time_t solve_e = time(NULL);
    mcf_time += (solve_e - cost_e);

    if(status == MinCostFlow::OPTIMAL) {  
        for (int u = 0; u < cur_nets.size(); ++u) {
            auto net = nets[cur_nets[u]];
            auto dbnet = net->getDbNet();
            for (int v = 0; v < valid_bts.size(); ++v) {
                if((btnum4nets[cur_nets[u]] > 0) && (costs[u][v] * 1.0 / (nets_length[cur_nets[u]] + 1e-8) > 0.5) )
                    continue;
                // if((btnum4nets[cur_nets[u]] > 0) && (costs[u][v] > 0) )
                //     continue;
                auto flow = min_cost_flow.Flow(u*valid_bts.size() + v);
                
                if(assigned[valid_bts[v]] && (flow > 0.5)) {
                    std::cout << "[ERROR] Bonding terminal " << valid_bts[v] << " has already been assigned with flow " << flow << " to net " << cur_nets[u] << std::endl;
                    exit(1);
                }

                if (flow > 0.5) {
                    // add pins to nets
                    // updateCost
                    assigned[valid_bts[v]] = true;
                    int bt_id = valid_bts[v];

                    std::string term_str = "BT"+std::to_string(bt_id);
                    
                    odb::dbBTerm* bterm = odb::dbBTerm::create(dbnet, term_str.c_str());
                    odb::Point bt = odb::Point((bt_id % int(std::ceil(width * 1.0 / BT_PITCH))) * BT_PITCH, std::floor(bt_id / std::ceil(width * 1.0 / BT_PITCH)) * BT_PITCH);
                    auto tech_layer = db->getTech()->findLayer("BT");

                    if (!bterm) {
                        std::cout << "[ORtool ERROR] Bonding terminal " << term_str << " create fail for net " << dbnet->getName() << std::endl;
                        exit(1);
                    }
                    // Added bpins
                    odb::dbBPin * btpin = odb::dbBPin::create(bterm);

                    odb::dbBox::create(btpin, tech_layer, bt.x() - (BT_PITCH / 2), bt.y() - (BT_PITCH / 2), bt.x() + (BT_PITCH / 2), bt.y() + (BT_PITCH / 2));
                    bterm->setSigType(odb::dbSigType::SIGNAL);
                    btpin->setPlacementStatus(odb::dbPlacementStatus::PLACED);
                    bterm->connect(dbnet);

                    btnum4nets[cur_nets[u]] += 1;
                    nets_length[cur_nets[u]] = nets_length[cur_nets[u]] + costs[u][v];
                    
                    break;
                }
            }
        }
    } else {
        std::cout << "[ORtool Warn] cannot obtain feasible solution." << std::endl;
        if(!update_flag)
            exit(1);
    }

    if(!bot_up) {
        if(topLeftTree) {
            topLeftTree->solveMatchMCFOrtool(bot_up);
        }
        if(botLeftTree) {
            botLeftTree->solveMatchMCFOrtool(bot_up);
        }
        if(topRightTree) {
            topRightTree->solveMatchMCFOrtool(bot_up);
        }
        if(botRightTree) {
            botRightTree->solveMatchMCFOrtool(bot_up);
        }
    }

    return;
}

void Quad::legalize() {
    // solve the sub problems
    if(topLeftTree) {
        topLeftTree->legalize();
    }
    if(botLeftTree) {
        botLeftTree->legalize();
    }
    if(topRightTree) {
        topRightTree->legalize();
    }
    if(botRightTree) {
        botRightTree->legalize();
    }

    if(cur_nets.size() == 0) {
        return;
    }
    
    int start_x = std::ceil(cur_area.xMin() * 1.0 / BT_PITCH) * BT_PITCH;
    int start_y = std::ceil(cur_area.yMin() * 1.0 / BT_PITCH) * BT_PITCH;

    int bt_col_num = std::ceil((cur_area.xMax() - start_x) * 1.0 / BT_PITCH);
    int bt_row_num = std::ceil((cur_area.yMax() - start_y) * 1.0 / BT_PITCH);

    std::vector<std::vector<int>> costs;
    costs.resize(cur_nets.size());

    time_t cost_s = time(NULL);

    std::vector<int> valid_bts;

    for(int row = 0; row < bt_row_num; row++) {
        int cur_y = start_y + row * BT_PITCH;

        for(int col = 0; col < bt_col_num; col++) {
            int cur_x = start_x + col * BT_PITCH;
            if(cur_x > cur_area.xMax())
                break;
            int cur_idx = ((cur_y / BT_PITCH)) * std::ceil(width * 1.0 / BT_PITCH) + (cur_x / BT_PITCH) ;
            if(!assigned[cur_idx]) {

                for(int i=0; i<cur_nets.size(); i++) {
                    if(assigned[cur_idx]) {
                        costs[i].push_back(std::numeric_limits<int>::max());
                        continue;
                    }
                    int x,y;
                    bterms[cur_nets[i]]->getFirstPinLocation(x, y);
                    costs[i].push_back(abs(x - cur_x) + abs(y - cur_y));
                }
                valid_bts.push_back(cur_idx);
            }
        }
    }

    if(valid_bts.size() < cur_nets.size()) {
        std::cout << "the " << (valid_bts.size()) << " bt candidates are not enough for " << cur_nets.size() << " nets." << std::endl;
        if(update_flag)
            return;
        exit(1);
    }

    bool a = std::set<int>(valid_bts.begin(),valid_bts.end()).size()!=valid_bts.size();
    if(a) {
        std::cout << " [ERROR] repeat valid bonding terminals. " << std::endl;
        std::cout << " [ERROR] valid_bts.size(): " << valid_bts.size() << " with start_x,y: " << start_x << ", " << start_y << " and xMax, yMax: " << cur_area.xMax() << ", " << cur_area.yMax() << " bt_row_num: " << bt_row_num << " bt_col_num: " << bt_col_num << " with unique: " << std::set<int>(valid_bts.begin(),valid_bts.end()).size() << std::endl;
        auto dup_bts = findDuplicates(valid_bts);
        for (auto dup_bt : dup_bts) {
            std::cout << dup_bt << std::endl;
        }
        exit(1);
    }

    // Instantiate a SimpleMinCostFlow solver.
    SimpleMinCostFlow min_cost_flow;

    // init MCF flow
    int vertex_num = cur_nets.size() + valid_bts.size() + 1;
    int min_cost;
    
    for (int from = 0; from < cur_nets.size(); from++) {
        min_cost_flow.SetNodeSupply(from, 1);
        for (int to = cur_nets.size(); to < vertex_num - 1; to++) {
            int arc = min_cost_flow.AddArcWithCapacityAndUnitCost(from, to, 1, costs[from][to - cur_nets.size()]);
            if (arc != from * costs[0].size() + (to-cur_nets.size())) LOG(FATAL) << "Internal error";
        }
    }
    
    for (int to = cur_nets.size(); to < vertex_num - 1; to++) {
        int arc = min_cost_flow.AddArcWithCapacityAndUnitCost(to, vertex_num-1, 1, 1);
        if (arc != costs[0].size()*cur_nets.size() + (to-cur_nets.size())) LOG(FATAL) << "Internal error";
    }

    min_cost_flow.SetNodeSupply(vertex_num - 1, -cur_nets.size());
    
    time_t cost_e = time(NULL);

    cost_time += (cost_e - cost_s);

    // Find the min cost flow.
    int status = min_cost_flow.Solve();
    
    time_t solve_e = time(NULL);
    mcf_time += (solve_e - cost_e);

    if(status == MinCostFlow::OPTIMAL) {

        auto block = db->getChip()->getBlock();

        for (int u = 0; u < cur_nets.size(); ++u) {
            auto net = nets[cur_nets[u]];
            auto dbnet = net->getDbNet();

            int x,y;
            bterms[cur_nets[u]]->getFirstPinLocation(x, y);
            
            for (int v = 0; v < valid_bts.size(); ++v) {
                
                auto flow = min_cost_flow.Flow(u*costs[0].size() + v);

                if(flow > 0.5 && assigned[valid_bts[v]]){
                    std::cout << "[ERROR] bt " << v << "has been assigned." << std::endl;
                    exit(1);
                }

                if(flow > 0.5) {

                    // add pins to nets
                    // updateCost

                    bterms[cur_nets[u]]->disconnect();
                    odb::dbBTerm::destroy(bterms[cur_nets[u]]);

                    assigned[valid_bts[v]] = true;
                    int bt_id = valid_bts[v];

                    std::string term_str = "BT"+std::to_string(bt_id);
                    
                    odb::dbBTerm* bterm = odb::dbBTerm::create(dbnet, term_str.c_str());
                    odb::Point bt = odb::Point((bt_id % int(std::ceil(width * 1.0 / BT_PITCH))) * BT_PITCH, std::floor(bt_id / std::ceil(width * 1.0 / BT_PITCH)) * BT_PITCH);
                    auto tech_layer = db->getTech()->findLayer("BT");

                    if (!bterm) {
                        std::cout << "[ORtool ERROR] Bonding terminal " << term_str << " create fail for net " << dbnet->getName() << std::endl;
                        exit(1);
                    }
                    // Added bpins
                    odb::dbBPin * btpin = odb::dbBPin::create(bterm);

                    odb::dbBox::create(btpin, tech_layer, bt.x() - (BT_PITCH / 2), bt.y() - (BT_PITCH / 2), bt.x() + (BT_PITCH / 2), bt.y() + (BT_PITCH / 2));
                    bterm->setSigType(odb::dbSigType::SIGNAL);
                    btpin->setPlacementStatus(odb::dbPlacementStatus::PLACED);
                    bterm->connect(dbnet);
                    
                    int disp = abs(x-bt.x()) + abs(y-bt.y());
                    if(disp < min_disp)
                        min_disp = disp;
                    if(disp > max_disp)
                        max_disp = disp;
                    total_disp += disp;
                }
            }
        }
    } else {
        std::cout << "[ORtool Warn] cannot obtain feasible solution." << std::endl;
        exit(1);
    }
}

void Quad::AssignMip() {
    // solve the sub problems
    if(topLeftTree) {
        topLeftTree->AssignMip();
    }
    if(botLeftTree) {
        botLeftTree->AssignMip();
    }
    if(topRightTree) {
        topRightTree->AssignMip();
    }
    if(botRightTree) {
        botRightTree->AssignMip();
    }

    if(cur_nets.size() == 0) {
        return;
    }

    time_t cost_s = time(NULL);
    // get all valid bterms
    std::vector<int> valid_bts;

    int start_x = std::ceil(cur_area.xMin() * 1.0 / BT_PITCH) * BT_PITCH;
    int start_y = std::ceil(cur_area.yMin() * 1.0 / BT_PITCH) * BT_PITCH;

    int bt_col_num = std::ceil((cur_area.xMax() - start_x) * 1.0 / BT_PITCH);
    int bt_row_num = std::ceil((cur_area.yMax() - start_y) * 1.0 / BT_PITCH);
    

    int skip_col  = 2;
    if((bt_col_num * bt_row_num > 1)) {
        skip_col = std::max(std::round((bt_col_num * bt_row_num) / (cur_nets.size() * 4.0)), 1.0);
    }

    for(int row = 0; row < bt_row_num; row++) {
        int cur_y = start_y + row * BT_PITCH;

        for(int col = 0; col < bt_col_num; col += skip_col) {
            int cur_x = start_x + col * BT_PITCH + (row % skip_col) * BT_PITCH;
            if(cur_x > cur_area.xMax())
                break;
            int cur_idx = ((cur_y / BT_PITCH)) * std::ceil(width * 1.0 / BT_PITCH) + (cur_x / BT_PITCH) ;
            if(!assigned[cur_idx]) {
                valid_bts.push_back(cur_idx);
            }
        }
    }

    if(valid_bts.size() < cur_nets.size()) {
        std::cout << "the " << (valid_bts.size()) << " bt candidates are not enough for " << cur_nets.size() << " nets." << std::endl;
        
        skip_col = skip_col / 2;

        valid_bts.clear();

        for(int row = 0; row < bt_row_num; row++) {
            int cur_y = start_y + row * BT_PITCH;

            for(int col = 0; col < bt_col_num; col += skip_col) {
                int cur_x = start_x + col * BT_PITCH + (row % skip_col) * BT_PITCH;
                if(cur_x > cur_area.xMax())
                    break;
                int cur_idx = ((cur_y / BT_PITCH)) * std::ceil(width * 1.0 / BT_PITCH) + (cur_x / BT_PITCH) ;
                if(!assigned[cur_idx]) {
                    valid_bts.push_back(cur_idx);
                }
            }
        }
        
        if(valid_bts.size() < cur_nets.size()) {
            std::cout << "the " << (valid_bts.size()) << " bt candidates are not enough for " << cur_nets.size() << " nets." << std::endl;
            if(update_flag)
                return;
            exit(1);
        }
            
    }

    bool a = std::set<int>(valid_bts.begin(),valid_bts.end()).size()!=valid_bts.size();
    if(a) {
        std::cout << " [ERROR] repeat valid bonding terminals. " << std::endl;
        std::cout << " [ERROR] valid_bts.size(): " << valid_bts.size() << " with start_x,y: " << start_x << ", " << start_y << " and xMax, yMax: " << cur_area.xMax() << ", " << cur_area.yMax() << " bt_row_num: " << bt_row_num << " bt_col_num: " << bt_col_num << " with unique: " << std::set<int>(valid_bts.begin(),valid_bts.end()).size() << std::endl;
        auto dup_bts = findDuplicates(valid_bts);
        for (auto dup_bt : dup_bts) {
            std::cout << dup_bt << std::endl;
        }
        exit(1);
    }

    // initial costs matrix
    std::vector<std::vector<std::vector<int>>> costs;

    costs.resize(cur_nets.size());
    for(int n = 0; n < cur_nets.size(); n++) {
        costs[n].resize(valid_bts.size());
        for(int i = 0; i < valid_bts.size(); i++) {
            costs[n][i].resize(i+1);
        }
        primMip(nets[cur_nets[n]], valid_bts, costs[n]);
    }

    const int num_nets = cur_nets.size();
    const int num_bts = valid_bts.size();

    // Solver
    // Create the mip solver with the Gurobi backend.
    std::unique_ptr<MPSolver> solver(MPSolver::CreateSolver("Gurobi"));
    if (!solver) {
        LOG(WARNING) << "Gurobi solver unavailable.";
        return;
    }

    std::cout << "[INFO] Solver is created." << std::endl; 
    // initial MIP solver for assignment
    // Variables
    // x[n][i][j] is an array of 0-1 variables, which will be 1
    // if net i is assigned to bonding terminal i and j.
    std::vector<std::vector<std::vector<const MPVariable*>>> x(
        num_nets, std::vector<std::vector<const MPVariable*>>(num_bts));
    for (int n = 0; n < num_nets; ++n) {
        for (int i = 0; i < num_bts; ++i) {
            x[n][i].resize(i+1);
            for(int j = 0; j < (i+1); ++j) {
                x[n][i][j] = solver->MakeIntVar(0, 1, "");
            }
        }
    }

    std::cout << "[INFO] Variables are created." << std::endl;

    // Constr
    // Each net should be assigned one candidate
    for(int n=0; n<num_nets; ++n){
        LinearExpr sum_candidate;
        for (int i = 0; i < num_bts; ++i) {
            for(int j = 0; j < (i+1); ++j) {
                sum_candidate += x[n][i][j];
            }
        }
        solver->MakeRowConstraint(sum_candidate == 1.0);
    }

    std::cout << "[INFO] Constr-1 is created." << std::endl;

    // Each bonding terminal can only assigned to one net
    for(int i=0; i < num_bts; ++i) {
        LinearExpr sum_assign;
        for(int n=0; n<num_nets; n++) {
            for(int idx=0; idx < num_bts; ++idx) {
                if(idx <= i)
                    sum_assign += x[n][i][idx];
                else
                    sum_assign += x[n][idx][i];
            }
        }
        solver->MakeRowConstraint(sum_assign <= 1.0);
    }

    std::cout << "[INFO] Constr-2 is created." << std::endl;

    // Objective
    MPObjective* const objective = solver->MutableObjective();
    for(int n=0; n<num_nets; ++n){
        for (int i = 0; i < num_bts; ++i) {
            for(int j = 0; j < (i+1); ++j) {
                objective->SetCoefficient(x[n][i][j], costs[n][i][j]);
            }
        }
    }
    objective->SetMinimization();

    std::cout << "[INFO] Obj is set." << std::endl;

    // Solve
    const MPSolver::ResultStatus result_status = solver->Solve();

    std::cout << "[INFO] MIP has been solved." << std::endl;

    // Print solution.
    // Check that the problem has a feasible solution.
    if (result_status != MPSolver::OPTIMAL &&
        result_status != MPSolver::FEASIBLE) {
        LOG(FATAL) << "No solution found in Mip.";
        exit(1);
    }

    // get the assignment results
    for(int n=0; n<num_nets; n++) {
        auto net = nets[cur_nets[n]];
        auto dbnet = net->getDbNet();
        
        for (int i = 0; i < num_bts; ++i) {
            for(int j = 0; j < (i+1); ++j) {
                if(x[n][i][j]->solution_value() > 0.5) {
                    if(assigned[valid_bts[i]] || assigned[valid_bts[j]]){
                        // error
                        std::cout << "[ERROR] Bonding terminal " << valid_bts[i] << " or " << valid_bts[j] << " has already been assigned to net " << cur_nets[n] << std::endl;
                        exit(1);
                    }

                    assigned[valid_bts[i]] = true;
                    assigned[valid_bts[j]] = true;

                    int bt_id_i = valid_bts[i];
                    int bt_id_j = valid_bts[j];

                    std::string term_str_i = "BT"+std::to_string(bt_id_i);
                    std::string term_str_j = "BT"+std::to_string(bt_id_j);
                    
                    odb::dbBTerm* bterm_i = odb::dbBTerm::create(dbnet, term_str_i.c_str());
                    odb::Point bt_i = odb::Point((bt_id_i % int(std::ceil(width * 1.0 / BT_PITCH))) * BT_PITCH, std::floor(bt_id_i / std::ceil(width * 1.0 / BT_PITCH)) * BT_PITCH);
                    auto tech_layer = db->getTech()->findLayer("BT");

                    if (!bterm_i) {
                        std::cout << "[ORtool ERROR] Bonding terminal " << term_str_i << " create fail for net " << dbnet->getName() << std::endl;
                        exit(1);
                    }
                    // Added bpins
                    odb::dbBPin * btpin_i = odb::dbBPin::create(bterm_i);

                    odb::dbBox::create(btpin_i, tech_layer, bt_i.x() - (BT_PITCH / 2), bt_i.y() - (BT_PITCH / 2), bt_i.x() + (BT_PITCH / 2), bt_i.y() + (BT_PITCH / 2));
                    bterm_i->setSigType(odb::dbSigType::SIGNAL);
                    btpin_i->setPlacementStatus(odb::dbPlacementStatus::PLACED);
                    bterm_i->connect(dbnet);

                    if(i != j) {
                        odb::dbBTerm* bterm_j = odb::dbBTerm::create(dbnet, term_str_j.c_str());
                        odb::Point bt_j= odb::Point((bt_id_j % int(std::ceil(width * 1.0 / BT_PITCH))) * BT_PITCH, std::floor(bt_id_j / std::ceil(width * 1.0 / BT_PITCH)) * BT_PITCH);
                        auto tech_layer = db->getTech()->findLayer("BT");

                        if (!bterm_j) {
                            std::cout << "[ORtool ERROR] Bonding terminal " << term_str_j << " create fail for net " << dbnet->getName() << std::endl;
                            exit(1);
                        }
                        // Added bpins
                        odb::dbBPin * btpin_j= odb::dbBPin::create(bterm_j);

                        odb::dbBox::create(btpin_j, tech_layer, bt_j.x() - (BT_PITCH / 2), bt_j.y() - (BT_PITCH / 2), bt_j.x() + (BT_PITCH / 2), bt_j.y() + (BT_PITCH / 2));
                        bterm_j->setSigType(odb::dbSigType::SIGNAL);
                        btpin_j->setPlacementStatus(odb::dbPlacementStatus::PLACED);
                        bterm_j->connect(dbnet);
                    }

                    break;
                }
            }
        }
    }
}

//===========================
// Helper functions
//===========================
void Quad::prim(Net* net, std::vector<int> valid_bts, std::vector<int>& assign_cost) {
    int node_num = net->getNumPins();
    auto pins = net->getPins();
    int root_idx = 0;

    assign_cost.resize(valid_bts.size());

    // init cost matrix
    std::vector<std::vector<int>> graph(node_num);

    #pragma omp parallel for
    for(int nid_1 = 0; nid_1 < node_num; ++nid_1) {
        graph[nid_1].resize(node_num);
        if(nid_1 != node_num - 1) {
        if(pins[nid_1].isDriver())
            root_idx = nid_1;
        }
    }

    #pragma omp parallel for collapse(2)
    for(int nid_1 = 0; nid_1 < node_num; ++nid_1) {
        for(int nid_2 = 0; nid_2 < node_num; ++nid_2) {
            if(nid_1 > nid_2)
                continue;

            if(nid_1 == nid_2) {
                graph[nid_1][nid_2] = 0;
            } else if((nid_1 != node_num - 1) && (nid_2 != node_num - 1)){
                auto p1 = pins[nid_1].getPosition();
                auto p2 = pins[nid_2].getPosition();
                
                int layer_1 = pins[nid_1].getConnectionLayer();
                int layer_2 = pins[nid_2].getConnectionLayer();

                if(layer_1 != layer_2) {
                    graph[nid_1][nid_2] = std::numeric_limits<int>::max();
                    graph[nid_2][nid_1] = std::numeric_limits<int>::max();
                } else {
                    graph[nid_1][nid_2] = (abs(p1.x()-p2.x())+abs(p1.y()-p2.y()));
                    graph[nid_2][nid_1] = (abs(p1.x()-p2.x())+abs(p1.y()-p2.y()));
                }
            }
        }
    }

    #pragma omp parallel for
    for(int i = 0; i < valid_bts.size(); ++i) {
        int bt_id = valid_bts[i];
        std::vector<int> bt_cost(node_num);
        // #pragma omp parallel for
        for(int nid_1 = 0; nid_1 < node_num; ++nid_1) {
            // if(pins[nid_1].getConnectionLayer() > 3 && pins[nid_1].getConnectionLayer() < 10) {
            //     bt_cost[nid_1] = std::numeric_limits<int>::max();
            //     continue;
            // }
            odb::Point bt = odb::Point((bt_id % int(std::ceil(width * 1.0 / BT_PITCH))) * BT_PITCH, std::floor(bt_id / std::ceil(width * 1.0 / BT_PITCH)) * BT_PITCH);
            auto p1 = pins[nid_1].getPosition();
            bt_cost[nid_1] = (int)(abs(p1.x()-bt.x()) + abs(p1.y() - bt.y()));
        }
        assign_cost[i]= primHelper(graph, root_idx, bt_cost);
    }

    // return assign_cost;
}

void Quad::primMip(Net* net, std::vector<int> valid_bts, std::vector<std::vector<int>>& assign_cost) {
    int node_num = net->getNumPins();
    auto pins = net->getPins();
    int root_idx = 0;

    // init cost matrix
    std::vector<std::vector<int>> graph(node_num);

    #pragma omp parallel for
    for(int nid_1 = 0; nid_1 < node_num; ++nid_1) {
        graph[nid_1].resize(node_num);
        if(nid_1 != node_num - 1) {
        if(pins[nid_1].isDriver())
            root_idx = nid_1;
        }
    }

    #pragma omp parallel for collapse(2)
    for(int nid_1 = 0; nid_1 < node_num; ++nid_1) {
        for(int nid_2 = 0; nid_2 < node_num; ++nid_2) {
            if(nid_1 > nid_2)
                continue;

            if(nid_1 == nid_2) {
                graph[nid_1][nid_2] = 0;
            } else if((nid_1 != node_num - 1) && (nid_2 != node_num - 1)){
                auto p1 = pins[nid_1].getPosition();
                auto p2 = pins[nid_2].getPosition();
                
                int layer_1 = pins[nid_1].getConnectionLayer();
                int layer_2 = pins[nid_2].getConnectionLayer();

                if(layer_1 != layer_2) {
                    graph[nid_1][nid_2] = std::numeric_limits<int>::max();
                    graph[nid_2][nid_1] = std::numeric_limits<int>::max();
                } else {
                    graph[nid_1][nid_2] = (abs(p1.x()-p2.x())+abs(p1.y()-p2.y()));
                    graph[nid_2][nid_1] = (abs(p1.x()-p2.x())+abs(p1.y()-p2.y()));
                }
            }
        }
    }

    #pragma omp parallel for
    for(int i = 0; i < valid_bts.size(); ++i) {
        int bt_id_i = valid_bts[i];
        std::vector<int> bt_cost_i(node_num);
        odb::Point bt_i = odb::Point((bt_id_i % int(std::ceil(width * 1.0 / BT_PITCH))) * BT_PITCH, std::floor(bt_id_i / std::ceil(width * 1.0 / BT_PITCH)) * BT_PITCH);
        for(int nid_1 = 0; nid_1 < node_num; ++nid_1) {
            auto p1 = pins[nid_1].getPosition();
            bt_cost_i[nid_1] = (int)(abs(p1.x()-bt_i.x()) + abs(p1.y() - bt_i.y()));
        }
        for(int j = 0; j < (i+1); ++j) {
            int bt_id_j = valid_bts[j];
            std::vector<int> bt_cost_j(node_num);
            odb::Point bt_j = odb::Point((bt_id_j % int(std::ceil(width * 1.0 / BT_PITCH))) * BT_PITCH, std::floor(bt_id_j / std::ceil(width * 1.0 / BT_PITCH)) * BT_PITCH);
            // #pragma omp parallel for
            for(int nid_1 = 0; nid_1 < node_num; ++nid_1) {
                auto p1 = pins[nid_1].getPosition();
                bt_cost_j[nid_1] = (int)(abs(p1.x()-bt_j.x()) + abs(p1.y() - bt_j.y()));
            }
            int bt_dist = abs(bt_i.x() - bt_j.x()) + abs(bt_i.y() - bt_j.y());
            assign_cost[i][j]= primHelper(graph, root_idx, bt_cost_i, bt_cost_j, bt_dist);
        }
        // process i == j
        assign_cost[i][i]= primHelper(graph, root_idx, bt_cost_i);
    }

}

int Quad::primHelper(std::vector<std::vector<int>>& graph, int root_idx, std::vector<int>& bt_cost) {
    int node_num = graph[0].size() + 1;
    int total_weight = 0;
    std::vector<bool> selected(node_num, false);
    std::vector<int> weight(node_num, std::numeric_limits<int>::max());
    std::vector<int> to(node_num, -1);

    weight[root_idx] = 0;

    for (int i=0; i<node_num; ++i) {
        int v = -1;
        for (int j = 0; j < node_num; ++j) {
            if (!selected[j] && (v == -1 || weight[j] < weight[v]))
                v = j;
        }

        if (weight[v] == std::numeric_limits<int>::max()) {
            std::cout << "No MST!" << std::endl;
            exit(0);
        }

        selected[v] = true;
        total_weight += weight[v];
        
        // #pragma omp parallel for
        for (int to_id = 0; to_id < node_num; ++to_id) {
            auto cost = 0;
            if((v == node_num - 1) && (to_id == (node_num - 1))) {
                cost = 0;
            }
            else if(v == (node_num - 1)) {
                cost = bt_cost[to_id];
            } else if (to_id == (node_num - 1)) {
                cost = bt_cost[v];
            } else {
                cost = graph[v][to_id];
            }

            if (cost < weight[to_id]) {
                weight[to_id] = cost;
                to[to_id] = v;
            }
        }
    }

    return total_weight;
}

int Quad::primHelper(std::vector<std::vector<int>>& graph, int root_idx, std::vector<int>& bt_cost_i, std::vector<int>& bt_cost_j, int bt_dist) {
    int node_num = graph[0].size() + 2;
    int total_weight = 0;
    std::vector<bool> selected(node_num, false);
    std::vector<int> weight(node_num, std::numeric_limits<int>::max());
    std::vector<int> to(node_num, -1);

    weight[root_idx] = 0;

    for (int i=0; i<node_num; ++i) {
        int v = -1;
        for (int j = 0; j < node_num; ++j) {
            if (!selected[j] && (v == -1 || weight[j] < weight[v]))
                v = j;
        }

        if (weight[v] == std::numeric_limits<int>::max()) {
            std::cout << "No MST!" << std::endl;
            exit(0);
        }

        selected[v] = true;
        total_weight += weight[v];
        
        // #pragma omp parallel for
        for (int to_id = 0; to_id < node_num; ++to_id) {
            auto cost = 0;
            if((v >= node_num - 2) && (to_id >= (node_num - 2))) {
                if(v == to_id)
                    cost = 0;
                else
                    cost = bt_dist;
            }
            else if(v == (node_num - 2)) {
                cost = bt_cost_i[to_id];
            } 
            else if(v == (node_num - 1)) {
                cost = bt_cost_j[to_id];
            }
            else if (to_id == (node_num - 2)) {
                cost = bt_cost_i[v];
            } 
            else if (to_id == (node_num - 1)) {
                cost = bt_cost_j[v];
            }
            else {
                cost = graph[v][to_id];
            }

            if (cost < weight[to_id]) {
                weight[to_id] = cost;
                to[to_id] = v;
            }
        }
    }
    return total_weight;
}

std::set<int> Quad::findDuplicates(std::vector<int> vec)        // 无引用，无常量
{
    std::set<int> duplicates;
    std::sort(vec.begin(), vec.end());
    std::set<int> distinct(vec.begin(), vec.end());
    std::set_difference(vec.begin(), vec.end(), distinct.begin(), distinct.end(),
        std::inserter(duplicates, duplicates.end()));
    return duplicates;
}

}
