/////////////////////////////////////////////////////////////////////////////
//
// BSD 3-Clause License
//
// Copyright (c) 2019, The Regents of the Chinese University of Hong Kong
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <boost/icl/interval.hpp>
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

#include "RoutingTracks.h"
#include "odb/db.h"

using boost::icl::interval;

namespace grt {

class F2FGrid
{
 public:
  F2FGrid() = default;
  ~F2FGrid() = default;

  void init(const odb::Rect& die_area,
            const int tile_size_top,
            const int tile_size_bot,
            const int x_grids_top,
            const int x_grids_bot,
            const int y_grids_top,
            const int y_grids_bot,
            const bool perfect_regular_x_top,
            const bool perfect_regular_x_bot,
            const bool perfect_regular_y_top,
            const bool perfect_regular_y_bot,
            const int num_layers);

  void clear();

  // different tile size, grids number
  int getTileSize(bool tier) const { return tile_size_[tier]; }

  int getXGrids(bool tier) const { return x_grids_[tier]; }
  int getYGrids(bool tier) const { return y_grids_[tier]; }

  bool isPerfectRegularX(bool tier) const { return perfect_regular_x_[tier]; }
  bool isPerfectRegularY(bool tier) const { return perfect_regular_y_[tier]; }

  void setPitchesInTile(const int pitches_in_tile, bool tier)
  {
    pitches_in_tile_[tier] = pitches_in_tile;
  }
  int getPitchesInTile(bool tier) const { return pitches_in_tile_[tier]; }

  // same die area in both up and bottom

  int getXMin() const { return die_area_.xMin(); }
  int getYMin() const { return die_area_.yMin(); }

  void setXMin(int x) { die_area_.set_xlo(x); }
  void setYMin(int y) { die_area_.set_ylo(y); }

  int getXMax() const { return die_area_.xMax(); }
  int getYMax() const { return die_area_.yMax(); }

  int getNumLayers() const { return num_layers_; }

  const std::vector<int>& getSpacings() const { return spacings_; }
  const std::vector<int>& getMinWidths() const { return min_widths_; }

  void addSpacing(int value, int layer) { spacings_[layer] = value; }
  void addMinWidth(int value, int layer) { min_widths_[layer] = value; }

  const std::vector<int>& getHorizontalEdgesCapacities()
  {
    return horizontal_edges_capacities_;
  };
  const std::vector<int>& getVerticalEdgesCapacities()
  {
    return vertical_edges_capacities_;
  };

  void addHorizontalCapacity(int value, int layer)
  {
    horizontal_edges_capacities_[layer] = value;
  }
  void addVerticalCapacity(int value, int layer)
  {
    vertical_edges_capacities_[layer] = value;
  }

  void updateHorizontalEdgesCapacities(int layer, int reduction)
  {
    horizontal_edges_capacities_[layer] = reduction;
  };
  void updateVerticalEdgesCapacities(int layer, int reduction)
  {
    vertical_edges_capacities_[layer] = reduction;
  };

  odb::Point getMiddle();
  const odb::Rect& getGridArea() const;

  odb::Point getPositionOnGrid(const odb::Point& position, bool tier);

  void getBlockedTiles(const odb::Rect& obstruction,
                       odb::Rect& first_tile_bds,
                       odb::Rect& last_tile_bds,
                       odb::Point& first_tile,
                       odb::Point& last_tile);

  int computeTileReduce(const odb::Rect& obs,
                        const odb::Rect& tile,
                        int track_space,
                        bool first,
                        odb::dbTechLayerDir direction);

  interval<int>::type computeTileReduceInterval(const odb::Rect& obs,
                                                const odb::Rect& tile,
                                                int track_space,
                                                bool first,
                                                odb::dbTechLayerDir direction);

 private:
  odb::Rect die_area_;
  std::vector<int> tile_size_; // 2
  std::vector<int> x_grids_; // 2
  std::vector<int> y_grids_; // 2
  std::vector<bool> perfect_regular_x_; // 2
  std::vector<bool> perfect_regular_y_; // 2
  int num_layers_;
  std::vector<int> pitches_in_tile_{15, 15}; // 2
  std::vector<int> spacings_;
  std::vector<int> min_widths_;
  std::vector<int> horizontal_edges_capacities_;
  std::vector<int> vertical_edges_capacities_;
};

}  // namespace grt
