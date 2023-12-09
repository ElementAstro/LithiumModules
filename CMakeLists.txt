cmake_minimum_required(VERSION 3.20)

project(lithium_modules C CXX)

add_subdirectory(lithium_image)

add_subdirectory(lithium_search)

add_subdirectory(hydrogen_core)

if(NOT WIN32)
    add_subdirectory(hydrogen_loader)
endif()