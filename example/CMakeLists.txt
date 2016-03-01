include_directories (. ..)

find_package (Boost COMPONENTS program_options) 

# Here all the examples are set. The source is file is assumed to be ${example}.cpp
set (examples 
example01
df_fk_2d
)

foreach (example ${examples})
    set(example_src ${example}.cpp)
    add_executable(${example} ${example_src} $<TARGET_OBJECTS:libgftools>)
    target_link_libraries(${example} ${Boost_LIBRARIES})
endforeach(example)
