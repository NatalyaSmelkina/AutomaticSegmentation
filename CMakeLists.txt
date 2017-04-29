cmake_minimum_required(VERSION 3.2)

project(AutomaticSegmentation)

# Find ITK
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

# Find AGTK
find_package(AGTK REQUIRED)
include(${AGTK_USE_FILE})

#Find openCL
find_package (OpenCL REQUIRED)

message("ITK_DIR ${ITK_DIR}")
message("AGTK_INCLUDE_DIRS ${AGTK_INCLUDE_DIRS}")

add_executable(extract-large-airways extract-large-airways.cxx)
target_link_libraries(extract-large-airways ${AGTK_LIBRARIES} ${ITK_LIBRARIES} ${OPENCL_LIBRARIES})

include_directories(${CMAKE_SOURCE_DIR}/Utils)
include_directories(${CMAKE_SOURCE_DIR}/Core)

add_subdirectory(${CMAKE_SOURCE_DIR}/Utils)
add_subdirectory(${CMAKE_SOURCE_DIR}/Core)