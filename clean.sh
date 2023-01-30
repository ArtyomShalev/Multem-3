#!/bin/bash
find . -name "CMakeFiles" -exec rm -rf "{}" \;
find . -name "CMakeCache.txt" -exec rm -rf "{}" \;
find . -name "Makefile" -exec rm -r "{}" \;
find . -name "cmake_install.cmake" -exec rm -r "{}" \;
find . -name "CTestTestfile.cmake" -exec rm -r "{}" \;

