cmake_minimum_required(VERSION 2.8.4)
set(CMAKE_VERBOSE_MAKEFILE on)

set(NFSIM_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/../nfsim/include")

set(INCLUDE_DIRS
    ${NFSIM_INCLUDE_DIR}/
    ${NFSIM_INCLUDE_DIR}/nauty24 
    ${NFSIM_INCLUDE_DIR}/NFutil 
    ${NFSIM_INCLUDE_DIR}/NFutil/MTrand 
    ${NFSIM_INCLUDE_DIR}/NFtest/transcription 
    ${NFSIM_INCLUDE_DIR}/NFtest/tlbr 
    ${NFSIM_INCLUDE_DIR}/NFtest/simple_system 
    ${NFSIM_INCLUDE_DIR}/NFtest/agentcell/cell 
    ${NFSIM_INCLUDE_DIR}/NFtest/agentcell 
    ${NFSIM_INCLUDE_DIR}/NFscheduler 
    ${NFSIM_INCLUDE_DIR}/NFreactions/transformations 
    ${NFSIM_INCLUDE_DIR}/NFreactions/reactions 
    ${NFSIM_INCLUDE_DIR}/NFreactions/reactantLists 
    ${NFSIM_INCLUDE_DIR}/NFreactions/mappings 
    ${NFSIM_INCLUDE_DIR}/NFreactions 
    ${NFSIM_INCLUDE_DIR}/NFoutput 
    ${NFSIM_INCLUDE_DIR}/NFinput 
    ${NFSIM_INCLUDE_DIR}/NFinput/TinyXML 
    ${NFSIM_INCLUDE_DIR}/NFinput/json 
    ${NFSIM_INCLUDE_DIR}/NFfunction/muParser 
    ${NFSIM_INCLUDE_DIR}/NFfunction 
    ${NFSIM_INCLUDE_DIR}/NFcore 
    ${NFSIM_INCLUDE_DIR}/NFcore/reactionSelector 
    ${NFSIM_INCLUDE_DIR}/NFcore/moleculeLists 
)


project(nfsim_c)

add_definitions("-fPIC")

set(SRC_FILES
    src/nfsim_c.cpp
    src/nfsim_c.h
    src/nfsim_c_structs.h
    src/nfsim_c_structs.cpp
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${INCLUDE_DIRS})

add_definitions("-Wno-deprecated-declarations")
add_definitions(-std=c++11)


# -- nfsimCInterface_static --

add_library(${PROJECT_NAME}_static STATIC ${SRC_FILES})


# -- nfsimCInterface (dynamic) --

add_library(${PROJECT_NAME} SHARED ${SRC_FILES})
target_link_libraries(${PROJECT_NAME} NFsim_static_fpic)
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE 1)
