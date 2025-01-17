#
# This file was genared by ./proj2cmake and will be overwritten on it's next run!
# Please put all configurations in the cmake_conf/*.cmake files.
#

INCLUDE("tinycthread.cmake")


ADD_LIBRARY(3rd_tinycthread STATIC
            ${3rd_tinycthread_SRC})

set_target_properties(3rd_tinycthread PROPERTIES FOLDER "3rd_party")