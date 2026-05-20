#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.30)
project(oPs LANGUAGES CXX)

#----------------------------------------------------------------------------
# ROOT dependency
find_package(ROOT REQUIRED COMPONENTS Minuit2 MathCore Spectrum)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# Locate sources and header# (CONFIGURE_DEPENDS makes CMake re-run if new files are added)s for this project
file(GLOB sources CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/src/*.h)

#----------------------------------------------------------------------------
# Add the executable
#
add_executable(oPs ${PROJECT_SOURCE_DIR}/src/main.cc ${sources} ${headers} )

#----------------------------------------------------------------------------
# Include directories
target_include_directories(oPs PRIVATE
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/external/cxxopts/include
    ${ROOT_INCLUDE_DIRS}
)

#----------------------------------------------------------------------------
# Link libraries
target_link_libraries(oPs ${ROOT_LIBRARIES})


