ADD_LIBRARY(easylife easylife.c)

#message(STATUS "Current source dir is ${CMAKE_CURRENT_SOURCE_DIR}")

target_include_directories(easylife PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})