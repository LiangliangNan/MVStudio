#
# This file was genared by ./proj2cmake and will be overwritten on it's next run!
# Please put all configurations in the cmake_conf/*.cmake files.
#

INCLUDE("sba-1.5.cmake")

ADD_LIBRARY(3rd_sba STATIC
            ${3rd_sba_SRC})

TARGET_LINK_LIBRARIES(3rd_sba
            ${3rd_sba_DEPS}
            ${3rd_sba_ADDITIONAL_DEPS}
            ${SOLUTION_STATIC_LIB_DEPS}
            ${SOLUTION_GENERAL_DEPS})

set_target_properties(3rd_sba PROPERTIES FOLDER "3rd_party")