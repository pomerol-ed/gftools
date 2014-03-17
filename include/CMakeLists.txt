include_directories (.)

set (base_src
    defaults.hpp
    tools.hpp
    grid_base.hpp 
    matsubara_grid.hpp 
    RealGrid.hpp 
    KMesh.hpp 
    Container.hpp Container.hxx
    GridObject.hpp GridObject.hxx
    )

add_library(libgftools OBJECT ${base_src})

install(DIRECTORY . DESTINATION include/gftools
          FILES_MATCHING PATTERN "*.hpp" PATTERN "*.hxx"
        )
