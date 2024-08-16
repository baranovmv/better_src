cmake_minimum_required(VERSION 3.27)
project(src)

set(CMAKE_CXX_STANDARD 20)

add_library(src SHARED converter.h
        src.cpp)
