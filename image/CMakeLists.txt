#
# This file was genared by ./proj2cmake and will be overwritten on it's next run!
# Please put all configurations in the cmake_conf/*.cmake files.
#

cmake_minimum_required(VERSION 2.8)


INCLUDE("image.cmake")

ADD_LIBRARY(image SHARED
            ${image_SRC})

TARGET_LINK_LIBRARIES(image
            ${image_DEPS}
            ${image_ADDITIONAL_DEPS}
            ${SOLUTION_SHARED_LIB_DEPS}
            ${SOLUTION_GENERAL_DEPS})

set(PNG_INCLUDE_DIR ${MVSTUDIO_THIRD_PARTY}/png)
set(JPEG_INCLUDE_DIR ${MVSTUDIO_THIRD_PARTY}/jpeg ${MVSTUDIO_BUILD_DIR}/3rd_party/jpeg)
target_include_directories(image PRIVATE ${PNG_INCLUDE_DIR} ${JPEG_INCLUDE_DIR})