add_library(easylife easylife/easylife.c)

#message(STATUS "Current source dir is ${CMAKE_CURRENT_SOURCE_DIR}")

target_include_directories(easylife
	PUBLIC
		$<BUILD_INTERFACE:${easylife_SOURCE_DIR}>
)