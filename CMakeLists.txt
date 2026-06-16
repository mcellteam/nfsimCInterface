cmake_minimum_required(VERSION 3.25...4.3)

project(nfsim_c CXX)

set(SRC_FILES
    src/nfsim_c.cpp
    src/nfsim_c.h
    src/nfsim_c_structs.h
    src/nfsim_c_structs.cpp
)

# static (no link deps) and shared (links the PIC nfsim) C interface variants
add_library(nfsim_c_static STATIC ${SRC_FILES})
add_library(nfsim_c        SHARED ${SRC_FILES})
add_library(NFsimC::NFsimC ALIAS nfsim_c_static)
target_link_libraries(nfsim_c PUBLIC NFsim::NFsim_pic)

foreach(_t nfsim_c_static nfsim_c)
  # PUBLIC: consumers include "nfsim_c.h" from this src directory.
  target_include_directories(${_t} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
  # PRIVATE: nfsim's own headers, pulled straight from the NFsim target so the
  # long include list is not duplicated here.
  target_include_directories(${_t} PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}
      $<TARGET_PROPERTY:NFsim::NFsim,INTERFACE_INCLUDE_DIRECTORIES>)
  target_compile_features(${_t} PUBLIC cxx_std_11)
  target_compile_options(${_t} PRIVATE
      $<$<NOT:$<PLATFORM_ID:Windows>>:-Wno-return-type;-Wno-sign-compare;-Wno-deprecated-declarations>)
  set_target_properties(${_t} PROPERTIES POSITION_INDEPENDENT_CODE ON)
endforeach()
