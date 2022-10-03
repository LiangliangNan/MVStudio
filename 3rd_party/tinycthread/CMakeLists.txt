#
# This file was genared by ./proj2cmake and will be overwritten on it's next run!
# Please put all configurations in the cmake_conf/*.cmake files.
#

cmake_minimum_required(VERSION 2.10)


INCLUDE("tinycthread.cmake")


ADD_LIBRARY(tinycthread STATIC
            ${3rd_tinycthread_SRC})

TARGET_LINK_LIBRARIES(tinycthread
            ${3rd_tinycthread_DEPS}
            ${3rd_tinycthread_ADDITIONAL_DEPS}
            ${SOLUTION_STATIC_LIB_DEPS}
            ${SOLUTION_GENERAL_DEPS})

set_target_properties(tinycthread PROPERTIES FOLDER "3rd_party")