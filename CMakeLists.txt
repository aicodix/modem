cmake_minimum_required(VERSION 2.8.12)

project(modem LANGUAGES CXX)

include_directories(../dsp)
include_directories(../code)

add_executable(decode decode.cc)
add_executable(encode encode.cc)
