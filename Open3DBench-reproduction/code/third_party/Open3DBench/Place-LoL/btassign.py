#!/usr/bin/env python3
"""Run BTAssign-WL/BTAssign-R on an Open3DBench LoL placement."""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path

LOCAL_SITE = Path(__file__).resolve().parents[2] / "python"
if LOCAL_SITE.is_dir():
    sys.path.insert(0, str(LOCAL_SITE))

import numpy as np
from scipy.optimize import linear_sum_assignment


@dataclass
class Net:
    name: str
    top: tuple[tuple[float, float], ...]
    bottom: tuple[tuple[float, float], ...]

    @property
    def pins(self):
        return self.top + self.bottom

    @property
    def bbox(self):
        xs = [p[0] for p in self.pins]
        ys = [p[1] for p in self.pins]
        return min(xs), min(ys), max(xs), max(ys)


@dataclass
class Node:
    region: tuple[int, int, int, int]
    nets: list[str]
    children: list["Node"]


def parse_input(path: Path):
    lines = path.read_text().splitlines()
    libs, instances, nets = {}, {}, {}
    die = terminal_size = terminal_spacing = None
    i = 0
    while i < len(lines):
        words = lines[i].split()
        if not words:
            i += 1
            continue
        if words[0] == "LibCell":
            name, width, height, count = words[1], int(words[2]), int(words[3]), int(words[4])
            pins = {}
            for _ in range(count):
                i += 1
                pin = lines[i].split()
                pins[pin[1]] = (int(pin[2]), int(pin[3]))
            libs[name] = (width, height, pins)
        elif words[0] == "DieSize":
            die = tuple(map(int, words[1:5]))
        elif words[0] == "TerminalSize":
            terminal_size = (int(words[1]), int(words[2]))
        elif words[0] == "TerminalSpacing":
            terminal_spacing = int(words[1])
        elif words[0] == "NumInstances":
            for _ in range(int(words[1])):
                i += 1
                inst = lines[i].split()
                instances[inst[1]] = inst[2]
        elif words[0] == "NumNets":
            for _ in range(int(words[1])):
                i += 1
                while not lines[i].strip():
                    i += 1
                header = lines[i].split()
                net_name, count = header[1], int(header[2])
                pins = []
                for _ in range(count):
                    i += 1
                    pins.append(lines[i].split()[1])
                nets[net_name] = pins
        i += 1
    if None in (die, terminal_size, terminal_spacing):
        raise ValueError("Incomplete DieSize/TerminalSize/TerminalSpacing in input")
    return libs, instances, nets, die, terminal_size, terminal_spacing


def parse_placement(path: Path):
    placed, terminals = {}, {}
    die = None
    for raw in path.read_text().splitlines():
        words = raw.split()
        if not words:
            continue
        if words[0] == "TopDiePlacement":
            die = "top"
        elif words[0] == "BottomDiePlacement":
            die = "bottom"
        elif words[0] == "Inst":
            orient = words[4] if len(words) > 4 else "R0"
            placed[words[1]] = (die, int(words[2]), int(words[3]), orient)
        elif words[0] == "Terminal":
            terminals[words[1]] = (int(words[2]), int(words[3]))
    return placed, terminals


def orient_pin(x, y, width, height, px, py, orient):
    if orient in ("R0", "R360"):
        return x + px, y + py
    if orient == "R90":
        return x + height - py, y + px
    if orient == "R180":
        return x + width - px, y + height - py
    if orient == "R270":
        return x + py, y + width - px
    raise ValueError(f"Unsupported orientation: {orient}")


def build_nets(libs, instances, net_pins, placed, terminal_names):
    result = {}
    for net_name in terminal_names:
        top, bottom = [], []
        for ref in net_pins[net_name]:
            inst_name, pin_name = ref.split("/", 1)
            if inst_name not in placed:
                continue
            die, x, y, orient = placed[inst_name]
            width, height, pins = libs[instances[inst_name]]
            point = orient_pin(x, y, width, height, *pins[pin_name], orient)
            (top if die == "top" else bottom).append(point)
        if not top or not bottom:
            raise ValueError(f"{net_name} is not a placed two-die net")
        result[net_name] = Net(net_name, tuple(top), tuple(bottom))
    return result


def make_grid(die, terminal_size, spacing):
    xl, yl, xh, yh = die
    width, height = terminal_size
    pitch_x, pitch_y = width + spacing, height + spacing
    start_x = xl + spacing + width // 2
    start_y = yl + spacing + height // 2
    end_x = xh - spacing - (width + 1) // 2
    end_y = yh - spacing - (height + 1) // 2
    points = [
        (x, y)
        for y in range(start_y, end_y + 1, pitch_y)
        for x in range(start_x, end_x + 1, pitch_x)
    ]
    return points, pitch_x, pitch_y


def build_tree(region, names, nets, stop):
    xl, yl, xh, yh = region
    if xh - xl < stop or yh - yl < stop:
        return Node(region, list(names), [])
    xm, ym = (xl + xh) // 2, (yl + yh) // 2
    regions = [(xl, yl, xm, ym), (xm, yl, xh, ym), (xl, ym, xm, yh), (xm, ym, xh, yh)]
    groups = [[] for _ in regions]
    here = []
    for name in names:
        bx0, by0, bx1, by1 = nets[name].bbox
        child = next(
            (idx for idx, (cx0, cy0, cx1, cy1) in enumerate(regions)
             if cx0 <= bx0 and cy0 <= by0 and bx1 <= cx1 and by1 <= cy1),
            None,
        )
        (groups[child] if child is not None else here).append(name)
    children = [build_tree(r, group, nets, stop) for r, group in zip(regions, groups) if group]
    return Node(region, here, children)


def wl_cost(net: Net, terminal_points, terminal_size):
    terminals = [(point[0], point[1], None) for point in terminal_points]
    nodes = [(x, y, 0) for x, y in net.top] + [(x, y, 1) for x, y in net.bottom] + terminals
    count = len(nodes)
    used = [False] * count
    distance = [math.inf] * count
    distance[0] = 0.0
    total = 0.0
    for _ in range(count):
        node = min((i for i in range(count) if not used[i]), key=distance.__getitem__)
        if math.isinf(distance[node]):
            raise ValueError(f"Disconnected net {net.name}")
        used[node] = True
        total += distance[node]
        x0, y0, layer0 = nodes[node]
        for other, (x1, y1, layer1) in enumerate(nodes):
            if used[other] or (layer0 is not None and layer1 is not None and layer0 != layer1):
                continue
            distance[other] = min(distance[other], abs(x0 - x1) + abs(y0 - y1))
    return total


def prudy_value(point, nets, terminal_size):
    cx, cy = point
    total = 0.0
    for net in nets.values():
        xl, yl, xh, yh = net.bbox
        if xl <= cx <= xh and yl <= cy <= yh:
            total += 1.0 / max(xh - xl, terminal_size[0]) + 1.0 / max(yh - yl, terminal_size[1])
    return total


def even_sample(items, count):
    if count >= len(items):
        return list(items)
    step = len(items) / count
    return [items[min(int((i + 0.5) * step), len(items) - 1)] for i in range(count)]


def select_candidates(available, count):
    ordered = sorted(available, key=lambda p: (p[1], p[0]))
    count = min(count, len(ordered))
    return even_sample(ordered, count)


def assign_wl(root, nets, grid, terminal_size, budget, sample_factors):
    used = set()
    assignments = {name: [] for name in nets}

    def solve(node):
        pending = []
        for child in node.children:
            pending.extend(solve(child))
        pending.extend(node.nets)
        pending = list(dict.fromkeys(pending))
        if not pending:
            return []
        xl, yl, xh, yh = node.region
        for iteration in range(budget):
            current = [
                name
                for name in pending
                if len(assignments[name]) == iteration
                and (
                    iteration == 0
                    or 10 < len(nets[name].pins) + len(assignments[name]) < 100
                )
            ]
            if not current:
                continue
            available = [p for p in grid if p not in used and xl <= p[0] < xh and yl <= p[1] < yh]
            if len(available) < len(current):
                return [name for name in pending if not assignments[name]]
            factor = sample_factors[min(iteration, len(sample_factors) - 1)]
            candidates = select_candidates(
                available, max(len(current), len(current) * factor)
            )
            wl_costs = np.empty((len(current), len(candidates)), dtype=np.float64)
            for row, name in enumerate(current):
                wl = np.fromiter(
                    (
                        wl_cost(nets[name], assignments[name] + [point], terminal_size)
                        for point in candidates
                    ),
                    dtype=np.float64,
                )
                wl_costs[row] = wl
            rows, cols = linear_sum_assignment(wl_costs)
            for row, col in zip(rows, cols):
                point = candidates[col]
                if iteration > 0:
                    old_wl = wl_cost(nets[current[row]], assignments[current[row]], terminal_size)
                    new_wl = wl_cost(
                        nets[current[row]], assignments[current[row]] + [point], terminal_size
                    )
                    if (new_wl - old_wl) / old_wl > 0.5:
                        continue
                assignments[current[row]].append(point)
                used.add(point)
        return [name for name in pending if not assignments[name]]

    unresolved = solve(root)
    if unresolved:
        raise RuntimeError(f"Insufficient candidates for {len(unresolved)} nets")
    return assignments


def refine_route_aware(assignments, nets, grid, terminal_size, prudy):
    """Improve the WL solution without degrading either proxy metric."""
    refined = {name: list(points) for name, points in assignments.items()}
    used = {point for points in refined.values() for point in points}
    ordered_grid = sorted(grid, key=lambda point: (prudy[point], point[1], point[0]))

    for name in sorted(refined, key=lambda value: int(value[1:])):
        for index, old_point in enumerate(tuple(refined[name])):
            old_wl = wl_cost(nets[name], refined[name], terminal_size)
            used.remove(old_point)
            candidates = [
                point
                for point in ordered_grid
                if point not in used and prudy[point] <= prudy[old_point]
            ]
            wl = np.fromiter(
                (
                    wl_cost(
                        nets[name],
                        refined[name][:index] + [point] + refined[name][index + 1:],
                        terminal_size,
                    )
                    for point in candidates
                ),
                dtype=np.float64,
            )
            demand = np.fromiter((prudy[point] for point in candidates), dtype=np.float64)
            beta = wl.min() / demand.mean() if demand.mean() > 0 else 0.0
            cost = wl + beta * demand
            cost[wl > old_wl + 1e-6] = math.inf
            choice = int(cost.argmin())
            refined[name][index] = candidates[choice]
            used.add(candidates[choice])
    return refined


def legality(coords, die, terminal_size, spacing):
    pitch_x, pitch_y = terminal_size[0] + spacing, terminal_size[1] + spacing
    start_x = die[0] + spacing + terminal_size[0] // 2
    start_y = die[1] + spacing + terminal_size[1] // 2
    aligned = sum((x - start_x) % pitch_x == 0 and (y - start_y) % pitch_y == 0 for x, y in coords)
    violations = 0
    ordered = sorted(coords)
    for i, (x0, y0) in enumerate(ordered):
        for x1, y1 in ordered[i + 1:]:
            if x1 - x0 >= pitch_x:
                break
            if abs(y1 - y0) < pitch_y:
                violations += 1
    return {"grid_aligned": aligned, "unique": len(set(coords)), "spacing_violations": violations}


def metrics(assignments, nets, terminal_size, die, spacing):
    coords = [point for points in assignments.values() for point in points]
    wl = sum(wl_cost(nets[name], points, terminal_size) for name, points in assignments.items())
    demand = sum(prudy_value(point, nets, terminal_size) for point in coords)
    return {
        "mst_proxy_um": wl / 100.0,
        "prudy_sum": demand,
        "prudy_mean": demand / len(coords),
        "terminal_count": len(coords),
        **legality(coords, die, terminal_size, spacing),
    }


def write_placement(source: Path, target: Path, assignments):
    output = []
    for raw in source.read_text().splitlines():
        words = raw.split()
        if words and words[0] == "Terminal" and words[1] in assignments:
            x, y = assignments[words[1]][0]
            raw = f"Terminal {words[1]} {x} {y}"
        output.append(raw)
    target.write_text("\n".join(output) + "\n")


HBT_RE = re.compile(r"^\s*-\s+(HBT_(?:BOT|TOP)IN_(\d+))\s+(\S+)")
PLACED_RE = re.compile(r"^(\s*\+\s+PLACED\s+\(\s*)-?\d+\s+-?\d+(\s*\)\s+\S+\s*;.*)$")
COMPONENTS_RE = re.compile(r"^(\s*COMPONENTS\s+)(\d+)(\s*;.*)$")
HBT_CONNECTION_RE = re.compile(r"^(\s*)\(\s*(HBT_(?:BOT|TOP)IN_\d+)\s+(BOT|TOP)\s*\)(.*)$")


def write_assignments(target: Path, assignments):
    ordered = {
        name: assignments[name]
        for name in sorted(assignments, key=lambda value: int(value[1:]))
    }
    target.write_text(json.dumps(ordered, indent=2) + "\n")


def write_def(source: Path, target: Path, assignments, baseline):
    lines = source.read_text().splitlines()
    units = next(int(line.split()[3]) for line in lines if line.startswith("UNITS DISTANCE MICRONS"))
    scale = units // 100
    component_start = next(i for i, line in enumerate(lines) if line.startswith("COMPONENTS "))
    component_end = next(i for i, line in enumerate(lines) if line.strip() == "END COMPONENTS")
    source_component_count = sum(
        line.startswith("  - ") for line in lines[component_start + 1:component_end]
    )
    def_coords, component_info, active = {}, {}, None
    for line in lines:
        match = HBT_RE.match(line)
        if match:
            active = f"N{int(match.group(2)) + 1}"
            component_info[active] = (match.group(1), match.group(3))
            continue
        placed = PLACED_RE.match(line)
        if active and placed:
            nums = re.search(r"\(\s*(-?\d+)\s+(-?\d+)\s*\)", line)
            def_coords[active] = (int(nums.group(1)), int(nums.group(2)))
            active = None
    common = next(name for name in baseline if name in def_coords)
    offset = (def_coords[common][0] - baseline[common][0] * scale,
              def_coords[common][1] - baseline[common][1] * scale)
    duplicates = {
        component_info[name][0]: [
            f"{component_info[name][0]}_K{iteration}"
            for iteration in range(2, len(points) + 1)
        ]
        for name, points in assignments.items()
    }
    extra_count = sum(len(points) - 1 for points in assignments.values())
    output, active, in_components = [], None, False
    for line in lines:
        header = COMPONENTS_RE.match(line)
        if header:
            line = f"{header.group(1)}{source_component_count + extra_count}{header.group(3)}"
            in_components = True
            output.append(line)
            continue
        if in_components and line.strip() == "END COMPONENTS":
            for name in sorted(assignments, key=lambda value: int(value[1:])):
                component, macro = component_info[name]
                for iteration, (x, y) in enumerate(assignments[name][1:], start=2):
                    output.append(f"  - {component}_K{iteration} {macro}")
                    output.append(
                        f"    + PLACED ( {x * scale + offset[0]} {y * scale + offset[1]} ) N ;"
                    )
            output.append(line)
            in_components = False
            continue
        match = HBT_RE.match(line)
        if in_components and match:
            active = f"N{int(match.group(2)) + 1}"
            output.append(line)
            continue
        placed = PLACED_RE.match(line)
        if active in assignments and placed:
            x, y = assignments[active][0]
            line = f"{placed.group(1)}{x * scale + offset[0]} {y * scale + offset[1]}{placed.group(2)}"
            active = None
        output.append(line)
        connection = HBT_CONNECTION_RE.match(line)
        if connection:
            for duplicate in duplicates.get(connection.group(2), []):
                output.append(
                    f"{connection.group(1)}( {duplicate} {connection.group(3)} ){connection.group(4)}"
                )
    target.write_text("\n".join(output) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--placement", type=Path, required=True)
    parser.add_argument("--def", dest="def_file", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--partition-um", type=float, default=50.0)
    parser.add_argument("--sample-factor", type=int, default=8)
    parser.add_argument("--later-sample-factor", type=int, default=4)
    parser.add_argument("--budget", type=int, default=2)
    args = parser.parse_args()
    if args.budget < 1:
        parser.error("--budget must be at least 1")

    libs, instances, net_pins, die, terminal_size, spacing = parse_input(args.input)
    placed, baseline = parse_placement(args.placement)
    nets = build_nets(libs, instances, net_pins, placed, baseline)
    grid, pitch_x, pitch_y = make_grid(die, terminal_size, spacing)
    root = build_tree(die, sorted(nets, key=lambda n: int(n[1:])), nets, int(args.partition_um * 100))
    prudy = {point: prudy_value(point, nets, terminal_size) for point in grid}

    args.output_dir.mkdir(parents=True, exist_ok=True)
    report = {
        "input": str(args.input),
        "placement": str(args.placement),
        "three_d_nets": len(nets),
        "candidate_count": len(grid),
        "terminal_size_um": [v / 100 for v in terminal_size],
        "pitch_um": [pitch_x / 100, pitch_y / 100],
        "partition_um": args.partition_um,
        "budget": args.budget,
        "sample_factors": [args.sample_factor, args.later_sample_factor],
        "second_terminal_max_increment_ratio": 0.5,
        "baseline": metrics(
            {name: [point] for name, point in baseline.items()},
            nets,
            terminal_size,
            die,
            spacing,
        ),
    }
    wl_assignments = None
    for mode in ("wl", "r"):
        started = time.perf_counter()
        if mode == "wl":
            assignments = assign_wl(
                root,
                nets,
                grid,
                terminal_size,
                args.budget,
                [args.sample_factor, args.later_sample_factor],
            )
            wl_assignments = assignments
        else:
            assignments = refine_route_aware(
                wl_assignments, nets, grid, terminal_size, prudy
            )
        stem = args.placement.stem + f"_btassign_{mode}_k{args.budget}"
        if args.budget == 1:
            write_placement(args.placement, args.output_dir / f"{stem}.txt", assignments)
        write_assignments(args.output_dir / f"{stem}.json", assignments)
        write_def(args.def_file, args.output_dir / f"{stem}.def", assignments, baseline)
        report[mode] = metrics(assignments, nets, terminal_size, die, spacing)
        report[mode]["runtime_s"] = time.perf_counter() - started
        report[mode]["mst_improvement_pct"] = (
            (report["baseline"]["mst_proxy_um"] - report[mode]["mst_proxy_um"])
            / report["baseline"]["mst_proxy_um"] * 100
        )
        if mode == "r":
            report[mode]["strategy"] = "WL-seeded constrained WL + adaptive-beta P-RUDY"
    (args.output_dir / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
