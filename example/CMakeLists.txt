include_directories (. ..)

# Here all the examples are set. The source is file is assumed to be ${example}.cpp
set (examples 
example01
#df1
#example02
#example03
)

foreach (example ${examples})
    set(example_src ${example}.cpp)
    add_executable(${example} ${example_src} $<TARGET_OBJECTS:libgftools>)
endforeach(example)
