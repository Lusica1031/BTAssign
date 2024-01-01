# BTAssign

This is a repo to implement F2F bonding terminal assignment algorithms into [FastRoute](https://github.com/The-OpenROAD-Project/OpenROAD/tree/master/src/grt)

#### Additional Package Requirement
ortools: [Install](https://developers.google.com/optimization/install/cpp)  (Recommend: binary installation)

#### Folder Structure
grt - BTAssign, which is integrated into FastRoute.
parser - make the netlist generated by Cadence Genus can be used in the contest binaries.

#### BTAssign
The key implementation can refer to grt/src/quadtree.cpp.

#### Usage
Replace the [FastRoute](https://github.com/The-OpenROAD-Project/OpenROAD/tree/master/src/grt) in OpenROAD, and directly use the OpenROAD framework.
