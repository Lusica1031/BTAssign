################################################################################
## Authors: Siting Liu
##          (Advisor: Bei Yu)
################################################################################

add_library(F2FRoute
  src/F2FRoute.cpp
  src/RSMT.cpp
  src/RipUp.cpp
  src/maze.cpp
  src/maze3D.cpp
  src/route.cpp
  src/utility.cpp
)

target_include_directories(F2FRoute
  PUBLIC
    include
  PRIVATE
    # Ugly but necessary to get GRoute.h and avoid a circular dependency
    # with FastRoute.  Remove once this data moves to OpenDB
    ../../include
  )

target_link_libraries(F2FRoute
  PRIVATE
    utl
    stt
    odb
    Boost::boost
)
