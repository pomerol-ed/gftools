include_directories (.)

set (base_src
    Common.h
    Grid.h Grid.hpp 
    Container.h Container.hpp
    GridObject.h GridObject.hpp
    )

add_library(gftoolslib OBJECT ${base_src})

install(DIRECTORY . DESTINATION include/fk
          FILES_MATCHING PATTERN "*.h *.hpp")
