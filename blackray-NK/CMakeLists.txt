cmake_minimum_required(VERSION 3.16)

project(blackray_deformed_kerr VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(BLACKRAY_NATIVE "Enable -march=native" ON)

find_package(OpenMP REQUIRED)

set(BLACKRAY_CORE_SOURCES
	src/metric.cpp
	src/diffeqs.cpp
	src/isco.cpp
	src/raytracer.cpp
	src/convolver.cpp
	src/blackray_c_api.cpp
)

add_library(blackray_core SHARED ${BLACKRAY_CORE_SOURCES})
target_include_directories(blackray_core
	PUBLIC
		${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_compile_options(blackray_core PRIVATE -O3)
target_link_libraries(blackray_core PUBLIC OpenMP::OpenMP_CXX)

if(BLACKRAY_NATIVE)
	target_compile_options(blackray_core PRIVATE -march=native)
endif()

add_executable(blackray_cli src/blackray_cli.cpp)
target_compile_options(blackray_cli PRIVATE -O3)
target_link_libraries(blackray_cli PRIVATE blackray_core)

if(BLACKRAY_NATIVE)
	target_compile_options(blackray_cli PRIVATE -march=native)
endif()

set_target_properties(blackray_core PROPERTIES
	OUTPUT_NAME blackray_core
	POSITION_INDEPENDENT_CODE ON
)

install(TARGETS blackray_core blackray_cli
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
)
