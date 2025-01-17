#
# This file was genared by ./proj2cmake and will be overwritten on it's next run!
# Please put all configurations in the cmake_conf/*.cmake files.
#

INCLUDE("3rd_png.cmake")

ADD_LIBRARY(3rd_png STATIC
            ${3rd_image_png_SRC})

target_include_directories(3rd_png PRIVATE ${MVSTUDIO_THIRD_PARTY}/zlib ${CMAKE_BINARY_DIR}/3rd_party/zlib)
set_target_properties(3rd_png PROPERTIES FOLDER "3rd_party")