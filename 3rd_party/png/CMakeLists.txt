#
# This file was genared by ./proj2cmake and will be overwritten on it's next run!
# Please put all configurations in the cmake_conf/*.cmake files.
#

cmake_minimum_required(VERSION 2.8)


INCLUDE("3rd_png.cmake")

ADD_LIBRARY(png STATIC
            ${3rd_image_png_SRC})

TARGET_LINK_LIBRARIES(3rd_image_png
            ${3rd_image_png_DEPS}
            ${3rd_image_png_ADDITIONAL_DEPS}
            ${SOLUTION_STATIC_LIB_DEPS}
            ${SOLUTION_GENERAL_DEPS})

