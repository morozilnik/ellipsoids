cmake_minimum_required(VERSION 3.0)
project (Diploma)

SET (PETUNIN_ELLIPSE_SRC
	helpers.hpp
	helpers.cpp
	Petunin.hpp
	Petunin.cpp
	HarPelet.hpp
	HarPelet.cpp
	Malandain.hpp
	Malandain.cpp)

add_library("PetuninEllipses"
	${PETUNIN_ELLIPSE_SRC})

SET(EIGEN_INCLUDE_DIR_HINTS "$ENV{CVNAR}/cvnar-external-libs/win32-x86/include")
LIST(APPEND CMAKE_MODULE_PATH "$ENV{CVNAR}/cvnar-external-libs/win32-x86/CMake")
find_package(Eigen REQUIRED)
include_directories(${EIGEN_INCLUDE_DIR})

add_executable("TestEllipses" test.cpp)
target_link_libraries("TestEllipses" "PetuninEllipses")