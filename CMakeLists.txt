cmake_minimum_required(VERSION 3.8)
project(all2nr)

add_subdirectory(./mylib)
add_executable(main main.cpp)
add_executable(maincsrel maincsrel.cpp)
add_executable(augment augment.cpp)
add_executable(augmentcsrel augmentcsrel.cpp)

target_link_libraries(main reliabilitylib)
target_link_libraries(maincsrel reliabilitylib)
target_link_libraries(augment reliabilitylib)
target_link_libraries(augmentcsrel reliabilitylib)