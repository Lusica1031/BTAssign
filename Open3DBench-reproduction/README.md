# BTAssign/Open3DBench 复现修改包

本目录只归档为 AES/cadb1021、K=2 复现新增或修改的代码、配置和小型验证结果，不包含生成的 DEF、ODB、SPEF、图片或第三方依赖。

## 文件清单

| 文件 | 类型 | 作用 |
| --- | --- | --- |
| `code/third_party/Open3DBench/Place-LoL/btassign.py` | 新增 | 读取 LoL 输入与放置，构造合法 HBT 网格，执行 K=2 BTAssign-WL 和受约束 BTAssign-R，写出 JSON/DEF 与代理指标。 |
| `code/third_party/Open3DBench/OpenROAD-3D/flow/run_btassign_aes.sh` | 新增 | 选择 `wl/r` 与 `grt/full`，调用 OpenROAD LoL flow。 |
| `code/third_party/Open3DBench/OpenROAD-3D/flow/designs/nangate45_3D/aes_cipher_top/config_btassign_wl_k2.mk` | 新增 | WL、K=2 的独立输入 DEF 与结果目录配置。 |
| `code/third_party/Open3DBench/OpenROAD-3D/flow/designs/nangate45_3D/aes_cipher_top/config_btassign_r_k2.mk` | 新增 | R、K=2 的独立输入 DEF 与结果目录配置。 |
| `patches/Open3DBench-flow-Makefile.patch` | 修改 | 相对原 Open3DBench Makefile 的精确两处修改：新增 `LOLFLOW_FINAL_TARGETS` 并接入 `do-lolflow`。 |
| `verification/report.json` | 验证 | AES/cadb1021、K=2 的本地合法性、MST 和 P-RUDY 报告。 |
| `verification/server-results.md` | 验证 | 服务器 FastRoute 与详细布线结果摘要。 |

## 核心实现改动

### BTAssign-WL

- 从 ICCAD/LoL 输入恢复上下层 pin 坐标和 3D nets。
- 按 `TerminalSize + TerminalSpacing` 建立合法中心点网格，并校验唯一性、网格对齐和间距。
- 使用 50 um quadtree、MST 线长代价和 Hungarian 二分图匹配。
- 第一轮采样因子为 8，第二轮为 4。
- K=2 第二终端只面向加入第一终端后 pin 数为 11 到 99 的 nets，并保留公开源码的 0.5 增量线长门限。
- 第二 HBT 以新增 DEF component 的方式同时接入 TOP/BOT net。

### BTAssign-R

- 采用论文的自适应 `WL + beta * P-RUDY`，其中 `beta` 由当前最小 WL 与候选平均 P-RUDY 确定。
- 从合法 WL 解开始重选终端。
- 只接受不增加该 net MST、且不增加终端 P-RUDY 的合法空闲位置，避免 Open3DBench 小规模用例上出现 R 全面退化。
- 这是针对 Open3DBench 的防退化适配；期刊版 BTAssign-R 源码未公开，因此不是作者原实现的逐行复制。

## 未归档内容

- `third_party/Open3DBench/Place-LoL/results/` 下的大型生成 DEF/JSON。
- 服务器生成的 ODB、DEF、SPEF、route guide、日志和图片。
- 未修改的 Open3DBench、OpenROAD 与第三方库文件。

服务器完整结果仍位于：

- WL：`/ugra/zhangmiao/Open3DBench/OpenROAD-3D/flow/results/nangate45_3D/aes/btassign_wl_k2_corrected/`
- R：`/ugra/zhangmiao/Open3DBench/OpenROAD-3D/flow/results/nangate45_3D/aes/btassign_r_k2_corrected/`
