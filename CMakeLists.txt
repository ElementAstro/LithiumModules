cmake_minimum_required(VERSION 3.20)

project(lithium_modules C CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/")
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake_modules/")

include(CheckIncludeFileCXX)

# add_subdirectory(lithium_image)

# add_subdirectory(lithium_search)

add_subdirectory(hydrogen_core)

add_subdirectory(hydrogen_client)

if(NOT WIN32)
    add_subdirectory(hydrogen_loader)
endif()