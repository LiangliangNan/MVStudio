#
# This file was genared by ./proj2cmake and will be overwritten on it's next run!
# Please put all configurations in the cmake_conf/*.cmake files.
#

cmake_minimum_required(VERSION 2.8)


INCLUDE("test.cmake")

ADD_EXECUTABLE(test
            ${test_SRC})

TARGET_LINK_LIBRARIES(test
            ${test_DEPS}
            ${test_ADDITIONAL_DEPS}
            ${SOLUTION_APP_DEPS}
            ${SOLUTION_GENERAL_DEPS})

