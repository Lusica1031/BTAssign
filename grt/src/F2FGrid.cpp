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

#include "Grid.h"

#include <complex>

namespace grt {

void F2FGrid::init(const odb::Rect& die_area,
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
                  const int num_layers)
{
  die_area_ = die_area;
  tile_size_.resize(2);
  x_grids_.resize(2);
  y_grids_.resize(2);
  perfect_regular_x_.resize(2);
  perfect_regular_y_.resize(2);

  tile_size_[0] = tile_size_bot;
  tile_size_[1] = tile_size_top;
  x_grids_[0] = x_grids_bot;
  x_grids_[1] = x_grids_top;
  y_grids_[0] = y_grids_bot;
  y_grids_[1] = y_grids_top;
  perfect_regular_x_[0] = perfect_regular_x_bot;
  perfect_regular_x_[1] = perfect_regular_x_top;
  perfect_regular_y_[0] = perfect_regular_y_bot;
  perfect_regular_y_[1] = perfect_regular_y_top;
  num_layers_ = num_layers;
  spacings_.resize(num_layers);
  min_widths_.resize(num_layers);
  horizontal_edges_capacities_.resize(num_layers);
  vertical_edges_capacities_.resize(num_layers);
}

void F2FGrid::clear()
{
  tile_size_.clear();
  x_grids_.clear();
  y_grids_.clear();
  perfect_regular_x_.clear();
  perfect_regular_y_.clear();
  spacings_.clear();
  min_widths_.clear();
  horizontal_edges_capacities_.clear();
  vertical_edges_capacities_.clear();
}

odb::Point F2FGrid::getPositionOnGrid(const odb::Point& position, bool tier)
{
  int x = position.x();
  int y = position.y();

  // Computing x and y center:
  int gcell_id_x = floor((float) ((x - die_area_.xMin()) / tile_size_[tier]));
  int gcell_id_y = floor((float) ((y - die_area_.yMin()) / tile_size_[tier]));

  if (gcell_id_x >= x_grids_[tier])
    gcell_id_x--;

  if (gcell_id_y >= y_grids_[tier])
    gcell_id_y--;

  int center_x
      = (gcell_id_x * tile_size_[tier]) + (tile_size_[tier] / 2) + die_area_.xMin();
  int center_y
      = (gcell_id_y * tile_size_[tier]) + (tile_size_[tier] / 2) + die_area_.yMin();

  return odb::Point(center_x, center_y);
}

void F2FGrid::getBlockedTiles(const odb::Rect& obstruction,
                           odb::Rect& first_tile_bds,
                           odb::Rect& last_tile_bds,
                           odb::Point& first_tile,
                           odb::Point& last_tile,
                           bool tier)
{
  odb::Point lower = obstruction.ll();  // lower bound of obstruction
  odb::Point upper = obstruction.ur();  // upper bound of obstruction

  lower
      = getPositionOnGrid(lower, tier);  // translate lower bound of obstruction to
                                   // the center of the tile where it is inside
  upper
      = getPositionOnGrid(upper, tier);  // translate upper bound of obstruction to
                                   // the center of the tile where it is inside

  // Get x and y indices of first blocked tile
  first_tile = {(lower.x() - (getTileSize(tier) / 2)) / getTileSize(tier),
                (lower.y() - (getTileSize(tier) / 2)) / getTileSize(tier)};

  // Get x and y indices of last blocked tile
  last_tile = {(upper.x() - (getTileSize(tier) / 2)) / getTileSize(tier),
               (upper.y() - (getTileSize(tier) / 2)) / getTileSize(tier)};

  odb::Point ll_first_tile = odb::Point(lower.x() - (getTileSize(tier) / 2),
                                        lower.y() - (getTileSize(tier) / 2));
  odb::Point ur_first_tile = odb::Point(lower.x() + (getTileSize(tier) / 2),
                                        lower.y() + (getTileSize(tier) / 2));

  odb::Point ll_last_tile = odb::Point(upper.x() - (getTileSize(tier) / 2),
                                       upper.y() - (getTileSize(tier) / 2));
  odb::Point ur_last_tile = odb::Point(upper.x() + (getTileSize(tier) / 2),
                                       upper.y() + (getTileSize(tier) / 2));

  if ((die_area_.xMax() - ur_last_tile.x()) / getTileSize(tier) < 1) {
    ur_last_tile.setX(die_area_.xMax());
  }
  if ((die_area_.yMax() - ur_last_tile.y()) / getTileSize(tier) < 1) {
    ur_last_tile.setY(die_area_.yMax());
  }

  first_tile_bds = odb::Rect(ll_first_tile, ur_first_tile);
  last_tile_bds = odb::Rect(ll_last_tile, ur_last_tile);
}

interval<int>::type F2FGrid::computeTileReduceInterval(
    const odb::Rect& obs,
    const odb::Rect& tile,
    int track_space,
    bool first,
    odb::dbTechLayerDir direction)
{
  int start_point, end_point;
  if (direction == odb::dbTechLayerDir::VERTICAL) {
    if (obs.xMin() >= tile.xMin() && obs.xMax() <= tile.xMax()) {
      start_point = obs.xMin();
      end_point = obs.xMax();
    } else if (first) {
      start_point = obs.xMin();
      end_point = tile.xMax();
    } else {
      start_point = tile.xMin();
      end_point = obs.xMax();
    }
  } else {
    if (obs.yMin() >= tile.yMin() && obs.yMax() <= tile.yMax()) {
      start_point = obs.yMin();
      end_point = obs.yMax();
    } else if (first) {
      start_point = obs.yMin();
      end_point = tile.yMax();
    } else {
      start_point = tile.yMin();
      end_point = obs.yMax();
    }
  }
  interval<int>::type reduce_interval(start_point, end_point);
  return reduce_interval;
}

int F2FGrid::computeTileReduce(const odb::Rect& obs,
                            const odb::Rect& tile,
                            int track_space,
                            bool first,
                            odb::dbTechLayerDir direction)
{
  int reduce = -1;
  if (direction == odb::dbTechLayerDir::VERTICAL) {
    if (obs.xMin() >= tile.xMin() && obs.xMax() <= tile.xMax()) {
      reduce = ceil(std::abs(obs.xMax() - obs.xMin()) / track_space);
    } else if (first) {
      reduce = ceil(std::abs(tile.xMax() - obs.xMin()) / track_space);
    } else {
      reduce = ceil(std::abs(obs.xMax() - tile.xMin()) / track_space);
    }
  } else {
    if (obs.yMin() >= tile.yMin() && obs.yMax() <= tile.yMax()) {
      reduce = ceil(std::abs(obs.yMax() - obs.yMin()) / track_space);
    } else if (first) {
      reduce = ceil(std::abs(tile.yMax() - obs.yMin()) / track_space);
    } else {
      reduce = ceil(std::abs(obs.yMax() - tile.yMin()) / track_space);
    }
  }

  return reduce;
}

odb::Point F2FGrid::getMiddle()
{
  return odb::Point((die_area_.xMin() + (die_area_.dx() / 2.0)),
                    (die_area_.yMin() + (die_area_.dy() / 2.0)));
}

const odb::Rect& F2FGrid::getGridArea() const
{
  return die_area_;
}

}  // namespace grt
