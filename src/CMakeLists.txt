add_library(RCAC RCAC.cpp)

add_executable(RCAC_binary main.cpp)
add_executable(Test Test_RCAC_functions.cpp)
#add_executable(maintest.cpp)

target_link_libraries(RCAC_binary RCAC)
target_link_libraries(Test RCAC)