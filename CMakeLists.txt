# Created by the script cgal_create_cmake_script (and then adapted manually).
# This is the CMake script for compiling a CGAL application.
# cmake ../ -DCMAKE_BUILD_TYPE=Debug

cmake_minimum_required(VERSION 3.1...3.20)
project(ifc2indoorgml)

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()
if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

## To add expensive tests
# add_definitions("-DCGAL_CHECK_EXPENSIVE")
# add_definitions("-Wall -Wextra")

## For profilling with gprof
# add_definitions("-pg")
# SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")

## To add an is_valid test after each operation (only in debug mode)
# add_definitions("-DCGAL_CMAP_TEST_VALID_REMOVALS")
# add_definitions("-DCGAL_CMAP_TEST_VALID_CONTRACTIONS")
# add_definitions("-DCGAL_CMAP_TEST_VALID_INSERTIONS")

#add_definitions("-std=c++11")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Option allowing to profile each operation of the demo (cout times on stdout).
add_definitions(-DCGAL_PROFILE_LCC_DEMO)


##################
find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt5)

find_package(Qt5 QUIET COMPONENTS Script OpenGL Svg)


#############################################################################
### IFC++, Carve (comes with IFC++) and OpenSceneGraph (OSG)
############################################################################
set(CMAKE_FIND_DEBUG_MODE TRUE)

#find_path(IFCPP_INCLUDE_DIRS)
set(IFCPP_SOURCE_DIR $ENV{IFCPP_DIR})

include_directories(
    ${IFCPP_SOURCE_DIR}/IfcPlusPlus/src
    ${IFCPP_SOURCE_DIR}/external/
    ${IFCPP_SOURCE_DIR}/external/Carve
    ${IFCPP_SOURCE_DIR}/external/Carve/src
    ${IFCPP_SOURCE_DIR}/external/Carve/src/lib
    ${IFCPP_SOURCE_DIR}/external/Carve/src/common
    ${IFCPP_SOURCE_DIR}/external/Carve/src/include
    $ENV{OSG_INCLUDE_DIR}
    $ENV{OSG_INCLUDE_REL_DIR}
    $ENV{OSG_INCLUDE_DEBUG_DIR}
)

#link_directories(
#    $ENV{IFCPP_LIB_DEBUG_DIR}
#    $ENV{IFCPP_LIB_DIR}
#    $ENV{IFCPP_CARVE_LIB_DIR}
#    $ENV{IFCPP_CARVE_LIB_DEBUG_DIR}
#)

IF(CMAKE_BUILD_TYPE EQUAL "Release")
    message("RELEASE MODE")
    find_library(IFCPP_LIBRARIES IfcPlusPlus $ENV{IFCPP_LIB_DIR})
    find_library(Carve_LIBRARIES carve $ENV{IFCPP_CARVE_LIB_DIR})
#    find_library(OSG_LIBRARIES osg osgUtil osgText $ENV{OSG_LIB_DIR})
ELSE()
    message("DEBUG MODE " $ENV{IFCPP_LIB_DEBUG_DIR})
    find_library(IFCPP_LIBRARIES IfcPlusPlusd $ENV{IFCPP_LIB_DEBUG_DIR})
    find_library(Carve_LIBRARIES carved $ENV{IFCPP_CARVE_LIB_DEBUG_DIR})
#    find_library(OSG_LIBRARIES osgd osgUtild osgTextd $ENV{OSG_LIB_DEBUG_DIR})
ENDIF()

find_package(OpenSceneGraph COMPONENTS osg osgUtil osgText)

set(CMAKE_FIND_DEBUG_MODE FALSE)

set(IFCPP_FOUND TRUE)

if(NOT (IFCPP_LIBRARIES AND Carve_LIBRARIES) )
    message("IFC++ is missing! IFC files will not be handled.")
    set(IFCPP_FOUND FALSE)
elseif(NOT OPENSCENEGRAPH_FOUND )
    message("OpenSceneGraph is missing! IFC files will not be handled.")
    set(IFCPP_FOUND FALSE)
else()
    add_definitions(-DIFCPP_ON)
    message("IFC++ and OpenSceneGraph sucessfully found!")
endif()


if(NOT
   (CGAL_Qt5_FOUND
    AND Qt5_FOUND
    AND (IFCPP_LIBRARIES AND Carve_LIBRARIES)))

  message(STATUS "NOTICE: This code requires CGAL, "
                 "Qt5, and IFC++ to compile and run.")

else()

  add_definitions(-DCGAL_USE_BASIC_VIEWER -DQT_NO_KEYWORDS)

  # ui file, created wih Qt Designer
  qt5_wrap_ui(uis
      MainWindow.ui

      # IFC++ components
      modules/ifc2lcc/UI/IFC_loader.ui
      modules/ifc2lcc/UI/load_spaces.ui

      # IndoorGML
      modules/lcc2indoorgml/UI/addthematiclayer.ui
  )

  # qrc files (resources files, that contain icons, at least)
  qt5_add_resources(CGAL_Qt5_RESOURCE_FILES ./ifc2indoorgml.qrc)


#  file(GLOB MODULES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/modules/"
#                             "${CMAKE_CURRENT_SOURCE_DIR}/modules/*")
#  foreach(SUB_DIR ${MODULES})
#    IF(IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/modules/${SUB_DIR}")
#       add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/modules/${SUB_DIR}")
#     endif()
#  endforeach()


  add_executable(ifc2indoorgml
    typedefs.h
    LCC_SpecialOps.h
    main.cpp
    Viewer.cpp
    MainWindow.cpp
    LCC_SpecialOps.cpp

    # IFC++ components
    modules/ifc2lcc/IFC2LCC.h
    modules/ifc2lcc/IFC_Importer.h
    modules/ifc2lcc/IFC_loader_ui.h
    modules/ifc2lcc/IFC2LCC.cpp
    modules/ifc2lcc/IFC_Importer.cpp
    modules/ifc2lcc/IFC_loader_ui.cpp

    # IndoorGML component
    modules/lcc2indoorgml/IndoorGML.h
    modules/lcc2indoorgml/IndoorGML_reader.h
    modules/lcc2indoorgml/IndoorGML_writer.h
    modules/lcc2indoorgml/LCC2IndoorGML_Ops.h
    modules/lcc2indoorgml/LCC2IndoorGML.h
    modules/lcc2indoorgml/LCC2IndoorGML.cpp
    modules/lcc2indoorgml/UI/addthematiclayer.h
    modules/lcc2indoorgml/UI/addthematiclayer.cpp
    # rapidXml for IndoorGML I/O
    modules/lcc2indoorgml/rapidxml/rapidxml.hpp
    modules/lcc2indoorgml/rapidxml/rapidxml_print.hpp
    modules/lcc2indoorgml/rapidxml/rapidxml_utils.hpp
    modules/lcc2indoorgml/rapidxml/rapidxml_iterators.hpp

    ${uis}
    ${CGAL_Qt5_RESOURCE_FILES}
    ${CGAL_Qt5_MOC_FILES}
    )

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS ifc2indoorgml)

  target_include_directories(ifc2indoorgml PUBLIC
      ${OPENSCENEGRAPH_INCLUDE_DIRS}
  )

  target_link_libraries(ifc2indoorgml PUBLIC
    CGAL::CGAL CGAL::CGAL_Qt5 Qt5::Gui Qt5::OpenGL
    # IFC++ components
    ${OPENSCENEGRAPH_LIBRARIES} ${IFCPP_LIBRARIES} ${Carve_LIBRARIES}
  )


#  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
#  cgal_add_compilation_test(ifc2indoorgml)
endif()
