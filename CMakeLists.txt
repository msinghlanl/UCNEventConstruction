cmake_minimum_required(VERSION 3.0)
project(Coincidence)

# Include the header files directory
include_directories(include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Add all your source files here
set(SOURCES 
    src/CoincidenceMultiplePmt.C
)

set(USE_ROOT TRUE)
set(CMAKE_CXX_STANDARD 17)

if (USE_ROOT)
    find_package(ROOT REQUIRED COMPONENTS Spectrum) # TODO add required components
    add_definitions(-DPLOT)

    include(${ROOT_USE_FILE})
endif ()

# Create the executable from the source files
add_executable(UCNConstruction ${SOURCES})
target_link_libraries(UCNConstruction ${ROOT_LIBRARIES})
target_link_libraries(UCNConstruction PRIVATE ${CMAKE_DL_LIBS})


# target_link_libraries(MyExecutable Boost::boost)