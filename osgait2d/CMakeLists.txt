project(DynamicWalker)

cmake_minimum_required(VERSION 2.6)

LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
FIND_PACKAGE(OpenSim REQUIRED)

SET(CMAKE_CONFIGURATION_TYPES "RelWithDebInfo;Release"
      CACHE STRING "Semicolon separated list of supported configuration types.
      Only supports Debug, Release, MinSizeRel, and RelWithDebInfo; anything
      else will be ignored." FORCE )

SET(CMAKE_CXX_FLAGS "-std=c++11")
SET(KIT DynamicWalker)
SET(UKIT DYNAMICWALKER)

set(DYNAMICWALKER_BUILD_DOC FALSE CACHE BOOL
    "Create API documentation (requires Doxygen).")

# Doxygen API documentation
if(DYNAMICWALKER_BUILD_DOC)
    find_package(Doxygen REQUIRED)
    if(DOXYGEN_FOUND)
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile
            ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile @ONLY)
        add_custom_target(doc
            ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
            COMMENT "Generating API documentation with Doxygen" VERBATIM)
        install(DIRECTORY ${PROJECT_BINARY_DIR}/doc/html/ DESTINATION doc)
    endif(DOXYGEN_FOUND)
endif(DYNAMICWALKER_BUILD_DOC)

INCLUDE_DIRECTORIES(${OPENSIMSIMBODY_INCLUDE_DIRS})

MESSAGE("${OPENSIMBODY_LIBRARIES}")
ADD_SUBDIRECTORY(DynamicWalker)
ADD_SUBDIRECTORY(tests)
