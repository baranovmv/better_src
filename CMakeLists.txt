cmake_minimum_required(VERSION 3.27)
project(src)

set(CMAKE_CXX_STANDARD 20)

function(checkout_submodule DIRECTORY)
  if(NOT EXISTS ${DIRECTORY}/.git)
    if(NOT GIT_EXECUTABLE)
      find_package(Git REQUIRED)
    endif()
    message(STATUS "Checking out ${DIRECTORY}...")
    execute_process(
      COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive .
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      RESULT_VARIABLE GIT_EXIT_CODE
      )
    if(NOT GIT_EXIT_CODE EQUAL "0")
      message(FATAL_ERROR "Failed to check out ${DIRECTORY}")
    endif()
  endif()
endfunction(checkout_submodule)
checkout_submodule(
  ${CMAKE_CURRENT_SOURCE_DIR}/external/ctypesgen
)
set(CTYPESGEN ${CMAKE_CURRENT_SOURCE_DIR}/external/ctypesgen)


set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_INSTALL_RPATH "$ORIGIN")
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

add_definitions(
  -D_POSIX_C_SOURCE=200809L
)

if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  add_compile_options(
    # enable
    -Wcast-qual
    -Wctor-dtor-privacy
    -Wformat-security
    -Wformat=2
    -Wlogical-op
    -Woverlength-strings
    -Wpointer-arith
    -Wstrict-null-sentinel
    -Wuninitialized
    # disable
    -Wno-psabi
    -Wno-stringop-overflow
    -Wno-system-headers
    -Wno-unused-function
    -Wno-unused-parameter
  )
endif()

add_compile_options(
  # enable assertions (even in release build
  -pthread
  -Werror
  #    -march=tigerlake
  #-march=core-avx2#tigerlake
  #    -march=x86-64
#  -march=armv6
#  --with-float=hard
#  --with-fpu=vfp
  -funwind-tables
  -fno-omit-frame-pointer
)
add_link_options(
  -pthread
)

if(${CMAKE_BUILD_TYPE} STREQUAL Release)
  add_compile_options(
    -UNDEBUG
    -O3
    -ggdb
  )
elseif (${CMAKE_BUILD_TYPE} STREQUAL Debug)
  add_compile_options(
    -DNDEBUG
    -O0
    -ggdb
  )
  add_compile_options(
    -fsanitize=address
    -fsanitize=undefined
  )
  add_link_options(
    -fsanitize=address
    -fsanitize=undefined
    -static-libasan
  )
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
add_library(src SHARED converter.h src.cpp)
get_target_property(src_DIR src LIBRARY_OUTPUT_DIRECTORY)

set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
include_directories(${INCLUDE_DIR})

find_package(Python3 REQUIRED)
add_custom_command(
  COMMENT "Generating python binding"
  DEPENDS ${INCLUDE_DIR}/src.h ${src}
  OUTPUT ${PROJECT_SOURCE_DIR}/tests/src.py
  COMMAND ${Python3_EXECUTABLE} ${CTYPESGEN}/run.py -o ${PROJECT_SOURCE_DIR}/tests/src.py
        -I ${INCLUDE_DIR} -L ${src_DIR} -l $<TARGET_FILE:src> -D __GNUC__ ${INCLUDE_DIR}/src.h
)
add_custom_target(tool_bin ALL
  DEPENDS ${PROJECT_SOURCE_DIR}/tests/src.py
)

add_executable(bench_src tests/benchmark.cpp)
target_link_libraries(bench_src src)

add_executable(bench_speex tests/benchmark_speex.cpp)
target_link_libraries(bench_speex speexdsp)
