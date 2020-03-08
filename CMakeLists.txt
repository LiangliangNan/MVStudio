cmake_minimum_required(VERSION 3.1)

get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})


# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

################################################################################

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

### Configuration
set(MVSTUDIO_ROOT         ${CMAKE_CURRENT_LIST_DIR})
set(MVSTUDIO_BUILD_DIR    ${CMAKE_CURRENT_BINARY_DIR})
set(MVSTUDIO_THIRD_PARTY  ${MVSTUDIO_ROOT}/3rd_party)
set(MVSTUDIO_INCLUDE_DIR  ${MVSTUDIO_ROOT})

ADD_SUBDIRECTORY("3rd_party/clapack")
ADD_SUBDIRECTORY("3rd_party/cminpack")
ADD_SUBDIRECTORY("3rd_party/cmvs-pmvs")
ADD_SUBDIRECTORY("3rd_party/graclus")
ADD_SUBDIRECTORY("3rd_party/jpeg")
ADD_SUBDIRECTORY("3rd_party/nlopt")
ADD_SUBDIRECTORY("3rd_party/pba")
ADD_SUBDIRECTORY("3rd_party/png")
ADD_SUBDIRECTORY("3rd_party/glog")
ADD_SUBDIRECTORY("3rd_party/sba")
ADD_SUBDIRECTORY("3rd_party/tinycthread")
ADD_SUBDIRECTORY("3rd_party/zlib")
ADD_SUBDIRECTORY("3rd_party/glew")
ADD_SUBDIRECTORY("3rd_party/rply")
ADD_SUBDIRECTORY("3rd_party/SiftGPU")
ADD_SUBDIRECTORY("algo")
ADD_SUBDIRECTORY("math")
ADD_SUBDIRECTORY("image")
ADD_SUBDIRECTORY("mvglib")
ADD_SUBDIRECTORY("sfm")
ADD_SUBDIRECTORY("MVStudio")

ADD_SUBDIRECTORY("easy3d")
ADD_SUBDIRECTORY("Tutorial_303_Viewer_Qt")

ADD_SUBDIRECTORY("test")


## Resources
# Copy resources dirs into our shadow build directory
file( COPY ${CMAKE_CURRENT_SOURCE_DIR}/resources DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} )
