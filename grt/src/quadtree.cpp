#include "quadtree.h"

namespace grt {

using utl::GRT;

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

// Init the quad tree
void Quad::update() {
    // partitioning boundary
    for (auto net = cur_nets.begin(); net != cur_nets.end();) {
        if(int(nets[(*net)]->getNumPins() / 2) <= btnum4nets[(*net)])
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

void Quad::solveMatchMCFOrtool() {
    if(cur_nets.size() == 0)
        return;

    // solve the sub problems
    if(topLeftTree) {
        topLeftTree->solveMatchMCFOrtool();
    }
    if(botLeftTree) {
        botLeftTree->solveMatchMCFOrtool();
    }
    if(topRightTree) {
        topRightTree->solveMatchMCFOrtool();
    }
    if(botRightTree) {
        botRightTree->solveMatchMCFOrtool();
    }
    
    time_t cost_s = time(NULL);
    // get all valid bterms
    std::vector<size_t> valid_bts;
    int width = db->getChip()->getBlock()->getDieArea().xMax();

    int start_x = std::ceil(cur_area.xMin() / BT_PITCH) * BT_PITCH;
    int start_y = std::ceil(cur_area.yMin() / BT_PITCH) * BT_PITCH;

    int bt_col_num = std::floor((cur_area.xMax() - start_x) / BT_PITCH);
    int bt_row_num = std::floor((cur_area.yMax() - start_y) / BT_PITCH);
    

    int skip_col  = 2;
    if((bt_col_num * bt_row_num > 1)) {
        skip_col = std::max(std::round((bt_col_num * bt_row_num) / (cur_nets.size() * 16)), 1.0);
    }

    for(int row = 0; row < bt_row_num; row++) {
        int cur_y = start_y + row * BT_PITCH;

        for(int col = 0; col < bt_col_num; col += skip_col) {
            int cur_x = start_x + col * BT_PITCH + (row % skip_col) * BT_PITCH;
            if(cur_x > width)
                break;
            size_t cur_idx = (int(cur_y / BT_PITCH)) * std::floor(width / BT_PITCH) + std::floor(cur_x / BT_PITCH);
            if(!assigned[cur_idx]) {
                valid_bts.push_back(cur_idx);
            }
        }
    }

    // calculate costs
    std::vector<std::vector<size_t>> costs(cur_nets.size());

    #pragma omp parallel for
    for(int i=0; i< cur_nets.size(); i++) {
        prim(nets[cur_nets[i]], valid_bts, costs[i]);
        for(int j=0; j<valid_bts.size(); j++) {
            costs[i][j] = costs[i][j] - nets_length[cur_nets[i]];
        }
    }

    // Instantiate a SimpleMinCostFlow solver.
    operations_research::SimpleMinCostFlow min_cost_flow;

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

    if(status == operations_research::MinCostFlow::OPTIMAL) {   
        for (int u = 0; u < cur_nets.size(); ++u) {
            auto net = nets[cur_nets[u]];
            auto dbnet = net->getDbNet();
            for (int v = 0; v < valid_bts.size(); ++v) {
                if(costs[u][v] > 0 && (btnum4nets[cur_nets[u]] > 0))
                    continue;
                auto flow = min_cost_flow.Flow(u*valid_bts.size() + v);
                // assigned[u][v] = flow;
                if(assigned[valid_bts[v]] && (flow > 0.5)) {
                    logger->error(GRT, 505, "Bonding terminal {} has already been assigned.", std::to_string(valid_bts[v]));
                    exit(1);
                }

                if (flow > 0.5) {
                    // add pins to nets
                    // updateCost
                    assigned[valid_bts[v]] = true;
                    int bt_id = valid_bts[v];
                    std::string term_str = "BT"+std::to_string(bt_id);
                    
                    odb::dbBTerm* bterm = odb::dbBTerm::create(dbnet, term_str.c_str());
                    odb::Point bt = odb::Point((bt_id % int(std::floor(width / BT_PITCH))) * BT_PITCH, std::floor(bt_id / std::floor(width / BT_PITCH)) * BT_PITCH);
                    auto tech_layer = db->getTech()->findRoutingLayer(10);

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
        std::cout << "[ORtool ERROR] cannot obtain feasible solution. " << std::endl;
        exit(1);
    }
}

// Init MCMF solver
// solve four sub-problems and then init current costs
// init current costs and flow
void Quad::solveMatchMCMF() {

    if(cur_nets.size() == 0)
        return;

    // solve the sub problems
    if(topLeftTree) {
        topLeftTree->solveMatchMCMF();
    }
    if(botLeftTree) {
        botLeftTree->solveMatchMCMF();
    }
    if(topRightTree) {
        topRightTree->solveMatchMCMF();
    }
    if(botRightTree) {
        botRightTree->solveMatchMCMF();
    }
    
    // get all valid bterms
    std::vector<size_t> valid_bts;
    int width = db->getChip()->getBlock()->getDieArea().xMax();

    int start_x = std::ceil(cur_area.xMin() / BT_PITCH) * BT_PITCH;
    int start_y = std::ceil(cur_area.yMin() / BT_PITCH) * BT_PITCH;

    int bt_col_num = std::floor((cur_area.xMax() - start_x) / BT_PITCH);
    int bt_row_num = std::floor((cur_area.yMax() - start_y) / BT_PITCH);
    

    int skip_col  = 2;
    if((bt_col_num * bt_row_num > 1)) {
        skip_col = std::max(std::round((bt_col_num * bt_row_num) / (cur_nets.size() * 8.0)), 1.0);
    }

    std::cout << "test-2 with skip_col " << skip_col << " " << (bt_col_num * bt_row_num) << " " << (cur_nets.size() * 10) << std::endl;

    for(int row = 0; row < bt_row_num; row++) {
        int cur_y = start_y + row * BT_PITCH;

        for(int col = 0; col < bt_col_num; col += skip_col) {
            int cur_x = start_x + col * BT_PITCH + (row % skip_col) * BT_PITCH;
            if(cur_x > width)
                break;
            size_t cur_idx = (int(cur_y / BT_PITCH)) * std::floor(width / BT_PITCH) + std::floor(cur_x / BT_PITCH);
            if(!assigned[cur_idx]) {
                valid_bts.push_back(cur_idx);
            }
        }
    }
    std::cout << "prim-0 with nets size " << cur_nets.size() << " and valid bts size: " << valid_bts.size() << std::endl;

    // calculate costs
    std::vector<std::vector<size_t>> costs(cur_nets.size());
    for(int i=0; i< cur_nets.size(); i++) {
        prim(nets[cur_nets[i]], valid_bts, costs[i]);
        // XWL(nets[cur_nets[i]], valid_bts, costs[i]);
    }

    std::cout << "prim-1" << std::endl;

    // #pragma omp parallel for
    // for(int i=0; i<cur_nets.size(); i++) {
    //     auto net = nets[cur_nets[i]];
    //     auto pins = net->getPins();
    //     costs[i].resize(valid_bts.size(), 0);
    //     #pragma omp parallel for
    //     for(int j=0; j<valid_bts.size(); j++) {
    //         int bt_id = valid_bts[j];
    //         odb::Point bt = odb::Point((bt_id % int(std::floor(width / BT_PITCH))) * BT_PITCH, std::ceil(bt_id / std::floor(width / BT_PITCH)) * BT_PITCH);
    //         for(auto pin : pins) {
    //             auto pin_pos = pin.getPosition();
    //             costs[i][j] += abs(pin_pos.x() - bt.x()) + abs(pin_pos.y() + bt.y());
    //         }
    //     }
    // }

    std::cout << "prim-2" << std::endl;

    // init MCMF flow
    int vertex_num = cur_nets.size() + valid_bts.size() + 2;
    MCMF<int> flow_network(vertex_num);

    int source = 0;
    int sink = vertex_num - 1;
    int min_cost;

    for (int from = 0; from < vertex_num; from++) {
        int to_s, to_e;
        bool isN2B = false;
        if (from == 0) {
            // source to nets
            isN2B = false;
            to_s = 1;
            to_e = costs.size() + 1;
        } else if (from < costs.size() + 1) {
            // nets to bterms
            isN2B = true;
            to_s = costs.size() + 1;
            to_e = vertex_num - 1;
        } else if (from < vertex_num - 1) {
            // bterms to sink
            isN2B = false;
            to_s = vertex_num - 1;
            to_e = vertex_num;
        } else
            continue;

        for (int to = to_s; to < to_e; to++) {
            if(isN2B)
                flow_network.addPath(from, to, costs[from-1][to-to_s], 1);
            else
                flow_network.addPath(from, to, 1, 1);
        }
    }

    std::cout << "finish mcmf initialization." << std::endl;

    const int max_flow = flow_network.findMCMF(source, sink, min_cost);

    std::cout << "finish mcmf." << std::endl;

    for (int u = 0; u < cur_nets.size(); ++u) {

        auto net = nets[cur_nets[u]];
        auto dbnet = net->getDbNet();
        for (int v = 0; v < valid_bts.size(); ++v) {
            auto flow = flow_network.getFlow(u+1, v+cur_nets.size()+1);
            // assigned[u][v] = flow;
            if(assigned[valid_bts[v]] && (flow > 0.5)) {
                logger->error(GRT, 502, "Bonding terminal {} has already been assigned.", std::to_string(valid_bts[v]));
                exit(1);
            }
            assigned[valid_bts[v]] = ((assigned[valid_bts[v]] + flow) > 0.5);

            if (flow > 0.5) {
                // add pins to nets
                // updateCost
                int bt_id = valid_bts[v];
                std::string term_str = "BT"+std::to_string(bt_id);
                
                odb::dbBTerm* bterm = odb::dbBTerm::create(dbnet, term_str.c_str());
                odb::Point bt = odb::Point((bt_id % int(std::floor(width / BT_PITCH))) * BT_PITCH, std::floor(bt_id / std::floor(width / BT_PITCH)) * BT_PITCH);
                // std::cout << "check-6 " << v << " with bt (" << bt.x() << ", " << bt.y() << ")" << std::endl;
                auto tech_layer = db->getTech()->findRoutingLayer(10);

                if (!bterm) {
                    std::cout << "[ERROR] Bonding terminal " << term_str << " create fail for net " << dbnet->getName() << std::endl;
                    exit(1);
                }
                // Added bpins
                odb::dbBPin * btpin = odb::dbBPin::create(bterm);

                odb::dbBox::create(btpin, tech_layer, bt.x() - (BT_PITCH / 2), bt.y() - (BT_PITCH / 2), bt.x() + (BT_PITCH / 2), bt.y() + (BT_PITCH / 2));
                bterm->setSigType(odb::dbSigType::SIGNAL);
                btpin->setPlacementStatus(odb::dbPlacementStatus::PLACED);
                bterm->connect(dbnet);
                // printf("flow value from %d -> %d: %d\n", u, v, flow);
            }
        }
    }
    flow_network.clear();
    std::cout << "finish bterm addition." << std::endl;
}

void Quad::XWL(Net* net, std::vector<size_t> valid_bts, std::vector<size_t>& assign_cost) {
    int node_num = net->getNumPins();
    auto pins = net->getPins();
    int root_idx = 0;
    int width = db->getChip()->getBlock()->getDieArea().xMax();

    assign_cost.resize(valid_bts.size(), 0);

    for(int i=0; i<valid_bts.size(); i++) {
        int bt_id = valid_bts[i];
        odb::Point bt = odb::Point((bt_id % int(std::floor(width / BT_PITCH))) * BT_PITCH, std::floor(bt_id / std::floor(width / BT_PITCH)) * BT_PITCH);
        for(auto pin : pins) {
            auto pin_pos = pin.getPosition();
            assign_cost[i] += abs(pin_pos.x() - bt.x()) + abs(pin_pos.y() + bt.y());
        }
    }

    // return assign_cost;
}

void Quad::prim(Net* net, std::vector<size_t> valid_bts, std::vector<size_t>& assign_cost) {
    int node_num = net->getNumPins();
    auto pins = net->getPins();
    int root_idx = 0;
    int width = db->getChip()->getBlock()->getDieArea().xMax();

    assign_cost.resize(valid_bts.size());

    // init cost matrix
    std::vector<std::vector<size_t>> graph(node_num);

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

                if(pins[nid_1].getConnectionLayer() != pins[nid_2].getConnectionLayer()) {
                    graph[nid_1][nid_2] = std::numeric_limits<size_t>::max();
                    graph[nid_2][nid_1] = std::numeric_limits<size_t>::max();
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
        std::vector<size_t> bt_cost(node_num);
        // #pragma omp parallel for
        for(int nid_1 = 0; nid_1 < node_num; ++nid_1) {
            odb::Point bt = odb::Point((bt_id % int(std::floor(width / BT_PITCH))) * BT_PITCH, std::floor(bt_id / std::floor(width / BT_PITCH)) * BT_PITCH);
            auto p1 = pins[nid_1].getPosition();
            bt_cost[nid_1] = (size_t)(abs(p1.x()-bt.x()) + abs(p1.y() - bt.y()));
        }
        assign_cost[i]= primHelper(graph, root_idx, bt_cost);
    }

    // return assign_cost;
}

int Quad::primHelper(std::vector<std::vector<size_t>>& graph, int root_idx, std::vector<size_t>& bt_cost) {
    int node_num = graph[0].size() + 1;
    int total_weight = 0;
    std::vector<bool> selected(node_num, false);
    std::vector<size_t> weight(node_num, std::numeric_limits<size_t>::max());
    std::vector<size_t> to(node_num, -1);

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

}