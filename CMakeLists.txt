cmake_minimum_required(VERSION 3.14)

project(RHITA_geom LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(QT NAMES Qt6 Qt5 REQUIRED COMPONENTS Core)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Core)

add_compile_definitions(QT_NO_KEYWORDS)

add_executable(rhita
  main.cpp
  geometry/ObjectPool.h
  geometry/mesh.cpp
  geometry/mesh.h
  geometry/node.cpp
  geometry/node.h
  geometry/element.cpp
  geometry/element.h
  geometry/czinsertiontool.cpp
  geometry/czinsertiontool.h
  geometry/cohesivezone.cpp
  geometry/cohesivezone.h

  geometry2d/node2d.h
  geometry2d/node2d.cpp
  geometry2d/mesh2d.cpp
  geometry2d/mesh2d.h
  geometry2d/element2d.cpp
  geometry2d/element2d.h
  geometry2d/czinsertiontool2d.cpp
  geometry2d/czinsertiontool2d.h
  geometry2d/cohesivezone2d.cpp
  geometry2d/cohesivezone2d.h

  generator.h
  generator.cpp
  generatoralt1.h
  generatoralt1.cpp
  generator3d.h
  generator3d.cpp
)

target_link_libraries(rhita Qt${QT_VERSION_MAJOR}::Core)
target_include_directories(rhita PRIVATE geometry geometry2d)

# Gmsh
target_include_directories(rhita PRIVATE "C:\\libraries\\gmsh-4.9.5-source\\api")
target_link_directories(rhita PRIVATE "C:\\libraries\\gmsh-build\\Release")
target_link_libraries(rhita gmsh)

# Eigen
find_package(Eigen3 REQUIRED NO_MODULE HINTS "C:\\libraries\\eigen-build")
target_link_libraries(rhita Eigen3::Eigen)

# spdlog
if (WIN32)
    target_link_libraries(rhita spdlog)
    target_include_directories(rhita PRIVATE "C:\\libraries\\spdlog\\include")
    target_link_directories(rhita PRIVATE "C:\\libraries\\spdlog-build\\Release")
else()
    target_link_libraries(rhita spdlog::spdlog)
endif()
