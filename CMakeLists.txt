cmake_minimum_required(VERSION 3.10)
project(HTSP)

set(CMAKE_CXX_STANDARD 14)

# boost
# ADD: Add "program_options" to COMPONENTS (critical fix for parameter parsing)
find_package(Boost 1.67.0 COMPONENTS system filesystem program_options REQUIRED)
if (Boost_FOUND)
    message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARY_DIRS}")
    message(STATUS "Boost_VERSION: ${Boost_VERSION}")
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
endif ()

add_definitions(-w) # omit warning information

add_executable(DPH2H main.cpp BasicFun.cpp PH2H.cpp)
add_executable(DH2H mainH2H.cpp BasicFunH2H.cpp H2H.cpp)

if (Boost_FOUND)
    if(UNIX AND NOT APPLE)
        # ADD: Link boost_program_options for Linux
        target_link_libraries(DPH2H pthread boost_system boost_thread boost_program_options)# For linux
        target_link_libraries(DH2H pthread boost_system boost_thread boost_program_options)# For linux
    endif()
    if(APPLE)
        # ADD: Link boost_program_options-mt for macOS (matches boost_thread-mt naming)
        target_link_libraries(DPH2H boost_system boost_thread-mt boost_program_options-mt)# For mac
        target_link_libraries(DH2H boost_system boost_thread-mt boost_program_options-mt)# For mac
    endif()
endif ()