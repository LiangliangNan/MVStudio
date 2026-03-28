#
# This file was genared by ./proj2cmake and will be overwritten on it's next run!
# Please put all configurations in the cmake_conf/*.cmake files.
#

include("sba-1.5.cmake")
add_library(3rd_sba STATIC ${3rd_sba_SRC})
target_link_libraries(3rd_sba 3rd_blas 3rd_f2c 3rd_clapack)
set_target_properties(3rd_sba PROPERTIES FOLDER "3rd_party")