add_executable(main
        utils.hpp
        graph.h
        graph.cpp
        lower_bounds.h
        lower_bounds.cpp
        reductions.h
        reductions.cpp
        algo.h
        algo.cpp
        main.cpp
)
set_target_properties(main PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
