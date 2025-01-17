# *************************************************************************
#      Copyright (C) 2015 Liangliang Nan <liangliang.nan@gmail.com>
#      https://3d.bk.tudelft.nl/liangliang/
#
#      This file is part of Easy3D. If it is useful in your research/work,
#      I would be grateful if you show your appreciation by citing it:
#      ------------------------------------------------------------------
#           Liangliang Nan.
#           Easy3D: a lightweight, easy-to-use, and efficient C++ library
#           for processing and rendering 3D data.
#           Journal of Open Source Software, 6(64), 3255, 2021.
#      ------------------------------------------------------------------
#
#      Easy3D is free software; you can redistribute it and/or modify
#      it under the terms of the GNU General Public License Version 3
#      as published by the Free Software Foundation.
#
#      Easy3D is distributed in the hope that it will be useful,
#      but WITHOUT ANY WARRANTY; without even the implied warranty of
#      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#      GNU General Public License for more details.
#
#      You should have received a copy of the GNU General Public License
#      along with this program. If not, see <http://www.gnu.org/licenses/>.
# *************************************************************************


# ------------------------------------------------------------------------------
# This file sets up Qt5 for CMake. When Qt5 was setup successfully, 'Qt5_FOUND'
# will be set. If Qt5 is not found, it will stop the configuration and show an
# error message.
# If Qt5 is found, it will set QtLibs to the corresponding Qt libraries, e.g.,
# Qt5Core, Qt5Gui, Qt5Widgets, Qt5OpenGL, Qt5Xml, etc.
#
# To use Qt5 in your program, you only need to include this file and specifying
# Qt libraries to link against, e.g.,
#       ------------------------------------------------------------------------
#           project(${PROJECT_NAME})
#           include( ../cmake/UseQt5.cmake )
#           add_executable(${PROJECT_NAME}, main.cpp)
#           target_link_libraries(${PROJECT_NAME} ${QtLibs})
#       ------------------------------------------------------------------------
# NOTE: 'UseQt5.cmake' must be included after you define your project but before
#       'add_executable()' or 'add_library()'.
#
#   The recommended way to specify libraries and headers with CMake is to use the
#   target_link_libraries command. This command automatically adds appropriate
#   include directories, compile definitions, the position-independent-code flag,
#   and links to the qtmain.lib library on Windows.
# ------------------------------------------------------------------------------

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
# Instruct CMake to run uic automatically when needed.
set(CMAKE_AUTOUIC ON)
# Instruct CMake to run rcc automatically when needed.
set(CMAKE_AUTORCC ON)

# This will find the Qt files.

find_package(Qt5 COMPONENTS Core Widgets OpenGL Xml QUIET)
if (Qt5_FOUND)
    message(STATUS "Found Qt5 version: ${Qt5Core_VERSION}")
    message(STATUS "Qt5 directory: ${Qt5_DIR}")
    set(QtLibs Qt5::Core Qt5::Widgets Qt5::OpenGL Qt5::Xml)
else()
    message(FATAL_ERROR "${PROJECT_NAME} requires Qt but Qt was not found. You can set 'Qt5_DIR' to the "
            "directory containing 'Qt5Config.cmake' or 'qt5-config.cmake' in 'CMakeCache.txt'. "
            "Optionally, you can set the Qt5 root directory 'QT5_ROOT_PATH' to the directory "
            "containing the 'bin' folder.")
endif()
