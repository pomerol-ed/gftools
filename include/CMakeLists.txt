include_directories (.)

set (base_src
    Defaults.hpp
    Tools.hpp
    Grid.hpp 
    MatsubaraGrid.hpp 
    RealGrid.hpp 
    KMesh.hpp 
    Container.hpp Container.hxx
    GridObject.hpp GridObject.hxx
    )

add_library(libgftools OBJECT ${base_src})

install(DIRECTORY . DESTINATION include/gftools
          FILES_MATCHING PATTERN "*.hpp *.hxx")
