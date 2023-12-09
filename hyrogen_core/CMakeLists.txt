cmake_minimum_required(VERSION 3.13)
project(hydrogencore C CXX)

find_package(Threads REQUIRED)

add_subdirectory(base)

include_directories(.)
include_directories(base)
include_directories(property)
include_directories(utils)

# Sources
list(APPEND ${PROJECT_NAME}_SOURCES
    property/hydrogenproperties.cpp
    property/hydrogenproperty.cpp
    property/hydrogenpropertybasic.cpp
    property/hydrogenpropertyview.cpp
    property/hydrogenpropertytext.cpp
    property/hydrogenpropertynumber.cpp
    property/hydrogenpropertyswitch.cpp
    property/hydrogenpropertylight.cpp
    property/hydrogenpropertyblob.cpp

    # connection/ttybase.cpp
)

# Headers
list(APPEND ${PROJECT_NAME}_HEADERS

    property/hydrogenproperties.h
    property/hydrogenproperty.h
    property/hydrogenpropertybasic.h
    property/hydrogenpropertyview.h
    property/hydrogenpropertytext.h
    property/hydrogenpropertynumber.h
    property/hydrogenpropertyswitch.h
    property/hydrogenpropertylight.h
    property/hydrogenpropertyblob.h

    property/hydrogenwidgetview.h
    property/hydrogenwidgettraits.h
)

# Private Headers
list(APPEND ${PROJECT_NAME}_PRIVATE_HEADERS
    property/hydrogenproperty_p.h
    property/hydrogenproperties_p.h
    property/hydrogenpropertyblob_p.h
    property/hydrogenpropertyswitch_p.h
    property/hydrogenpropertylight_p.h
    property/hydrogenpropertytext_p.h
    property/hydrogenpropertynumber_p.h
    property/hydrogenpropertybasic_p.h
    # TODO
)

# Build Object Library
add_library(${PROJECT_NAME}_OBJECT OBJECT)
set_property(TARGET ${PROJECT_NAME}_OBJECT PROPERTY POSITION_INDEPENDENT_CODE 1)

target_compile_definitions(${PROJECT_NAME}_OBJECT PRIVATE "-DHAVE_LIBNOVA")

if(WIN32)
target_link_libraries(${PROJECT_NAME}_OBJECT setupapi wsock32 ws2_32 shlwapi iphlpapi)
endif()

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