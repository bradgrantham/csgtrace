cmake_minimum_required(VERSION 3.5)


## Project

project(csgtrace)

## Project targets

add_executable(csgtrace main.cpp)
set_property(TARGET csgtrace PROPERTY CXX_STANDARD 17)

