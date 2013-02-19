include_directories (. ../include)

# Here all the examples are set. The source is file is assumed to be ${example}.cpp
set (examples 
example01
)

foreach (example ${examples})
    set(example_src ${example}.cpp)
    add_executable(${example} ${example_src} $<TARGET_OBJECTS:libgftools>)
endforeach(example)
