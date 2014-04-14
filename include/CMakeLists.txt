include_directories (.)

set (base_src
    defaults.hpp
    tools.hpp
    num_io.hpp
    grid_base.hpp 
    matsubara_grid.hpp 
    real_grid.hpp 
    kmesh.hpp 
    container.hpp container.hxx
    #GridObject.hpp GridObject.hxx
    )

add_library(libgftools OBJECT ${base_src})

install(DIRECTORY . DESTINATION include/gftools
          FILES_MATCHING PATTERN "*.hpp" PATTERN "*.hxx"
        )
