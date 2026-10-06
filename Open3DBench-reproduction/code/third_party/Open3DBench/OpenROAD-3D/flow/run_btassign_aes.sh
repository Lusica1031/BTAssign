#!/usr/bin/env bash

set -euo pipefail

mode=${1:-}
stage=${2:-full}
if [[ "$mode" != "wl" && "$mode" != "r" ]] || [[ "$stage" != "grt" && "$stage" != "full" ]]; then
    echo "Usage: $0 <wl|r> [grt|full]" >&2
    exit 2
fi

if [[ "$stage" == "grt" ]]; then
    final_targets=do-5_1_grt
else
    final_targets="do-route do-finish"
fi

flow_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
place_lol_root=$(cd "$flow_dir/../../Place-LoL" && pwd)
input_def="$place_lol_root/results/btassign/aes_cadb1021_k2/aes_btassign_${mode}_k2.def"
design_config="designs/nangate45_3D/aes_cipher_top/config_btassign_${mode}_k2.mk"
openroad_exe=${OPENROAD_EXE:-$(command -v openroad || true)}

if [[ ! -f "$input_def" ]]; then
    echo "Missing K=2 DEF: $input_def" >&2
    exit 1
fi
if [[ -z "$openroad_exe" ]]; then
    echo "OpenROAD is not available. Run this script inside the Open3DBench evaluation container." >&2
    exit 1
fi

cd "$flow_dir"
make \
    DESIGN_CONFIG="$design_config" \
    LOLFLOW_FINAL_TARGETS="$final_targets" \
    OPENROAD_EXE="$openroad_exe" \
    do-lolflow
