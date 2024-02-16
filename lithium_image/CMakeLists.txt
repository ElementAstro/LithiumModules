cmake_minimum_required(VERSION 3.20)
project(lithium_image C CXX)

include(CheckIncludeFileCXX)

option(ENABLE_DEBUG "Enable Debug Mode" OFF)
if(ENABLE_DEBUG)
    set(ENABLE_DEBUG_FLAG "1")
endif()

option(ENABLE_OPENCV_FLAG "Enable OpenCV" OFF)
if (ENABLE_OPENCV_FLAG)
    find_package(OpenCV REQUIRED)
    message("-- Found OpenCV ${OpenCV_VERSION}: ${OpenCV_LIBRARIES}")
else()
    message("-- OpenCV Not Found. Will not build opencv image module.")
endif()


option(ENABLE_CIMG_FLAG "Enable CImg" On)
if(ENABLE_CIMG_FLAG)
    set(ENABLE_CIMG "1")
endif()

option(ENABLE_STB_FLAG "Enable STB" OFF)
if(ENABLE_STB_FLAG)
    set(ENABLE_STB "1")
endif()

configure_file(config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

CHECK_INCLUDE_FILE_CXX(format HAS_STD_FORMAT)

if(NOT HAS_STD_FORMAT)
    find_package(fmt REQUIRED)
endif()

list(APPEND ${PROJECT_NAME}_LIBS
    atom-component
    loguru
    ${CFITSIO_LIBRARIES}
)

if(ENABLE_OPENCV_FLAG)
    list(APPEND ${PROJECT_NAME}_LIBS
        ${OpenCV_LIBRARIES}
    )
endif()
# Sources
list(APPEND ${PROJECT_NAME}_SOURCES
    main.cpp
)

if(ENABLE_OPENCV_FLAG)
    list(APPEND ${PROJECT_NAME}_SOURCES
        opencv/image.cpp
    )
endif()

if(ENABLE_CIMG_FLAG)
    list(APPEND ${PROJECT_NAME}_SOURCES
        cimg/image.cpp
        cimg/draw.cpp
        cimg/gaussian.cpp
    )
endif()

if(ENABLE_OPENCV_FLAG)
    list(APPEND ${PROJECT_NAME}_HEADERS
        opencv/image.hpp
    )
endif()

if(ENABLE_CIMG_FLAG)
    list(APPEND ${PROJECT_NAME}_HEADERS
        cimg/image.hpp
    )
endif()

# Private Headers
list(APPEND ${PROJECT_NAME}_PRIVATE_HEADERS

)

# Build Object Library
add_library(${PROJECT_NAME}_OBJECT OBJECT)
set_property(TARGET ${PROJECT_NAME}_OBJECT PROPERTY POSITION_INDEPENDENT_CODE 1)

target_link_libraries(${PROJECT_NAME}_OBJECT loguru atomstatic)

target_sources(${PROJECT_NAME}_OBJECT
    PUBLIC
    ${${PROJECT_NAME}_HEADERS}
    PRIVATE
    ${${PROJECT_NAME}_SOURCES}
    ${${PROJECT_NAME}_PRIVATE_HEADERS}
)

target_link_libraries(${PROJECT_NAME}_OBJECT ${${PROJECT_NAME}_LIBS})

add_library(${PROJECT_NAME}static STATIC)

target_link_libraries(${PROJECT_NAME}static ${PROJECT_NAME}_OBJECT ${${PROJECT_NAME}_LIBS})
target_link_libraries(${PROJECT_NAME}static ${CMAKE_THREAD_LIBS_INIT})
target_include_directories(${PROJECT_NAME}static PUBLIC .)

set_target_properties(${PROJECT_NAME}static PROPERTIES
    VERSION ${CMAKE_HYDROGEN_VERSION_STRING}
    SOVERSION ${HYDROGEN_SOVERSION}
    OUTPUT_NAME ${PROJECT_NAME} # this same name like shared library - backwards compatibility
)

install(TARGETS ${PROJECT_NAME}static
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

add_library(${PROJECT_NAME}shared SHARED)

target_link_libraries(${PROJECT_NAME}shared ${PROJECT_NAME}_OBJECT ${${PROJECT_NAME}_LIBS})
target_link_libraries(${PROJECT_NAME}shared ${CMAKE_THREAD_LIBS_INIT})
target_include_directories(${PROJECT_NAME}shared PUBLIC .)
set_target_properties(${PROJECT_NAME}shared PROPERTIES
    VERSION ${CMAKE_HYDROGEN_VERSION_STRING}
    SOVERSION ${HYDROGEN_SOVERSION}
    OUTPUT_NAME ${PROJECT_NAME} # this same name like shared library - backwards compatibility
)
install(TARGETS ${PROJECT_NAME}shared
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)