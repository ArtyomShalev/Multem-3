# faddeeva_fortran

This project provides the [MIT Faddeeva package](http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package) 
interface to Fortran programming language. 

## Description

Efficient and correct calculation of the special functions is relevant in every field of science. 
S.G. Johnson MIT implementation of various error functions is the best
implementation according to [REVIEW OF CPU AND GPU FADDEEVA IMPLEMENTATIONS](https://inspirehep.net/literature/1470416)

## How to use the project

We recommend to use Linux operating system to run Faddeeva_fortran.
This project is not independent. It provides the Fortran module with
various error functions of arbitrary complex arguments. 

### Prerequisites
* cmake 3.16.3 or higher
* GNU Fortran/gcc (tested for Ubuntu 20.04.1) or other fortran/c++ compilers

### Installation and linking to your project
See [example project](https://github.com/ArtyomShalev/fortran_faddeeva_usage_example.git)
or follow the instructions below
<pre><code>
git clone https://github.com/ArtyomShalev/faddeeva_fortran.git
</pre></code>
or as git submodule
<pre><code>
cd >your_main_project<
git submodule add https://github.com/ArtyomShalev/faddeeva_fortran.git 
</pre></code>
Add lines to your main CMakeLists.txt file
<pre><code>
add_subdirectory(faddeeva_fortran)
target_link_libraries(>your targets< faddeeva_fortran::faddeeva_fortran) 
</pre></code>

### Local installation
<pre><code>
git clone https://github.com/ArtyomShalev/faddeeva_fortran.git
cd faddeeva_fortran
cmake .
make
</pre></code>

## Documentation
To generate documentation in html and latex formats 
doxygen should be installed and one needs:
<pre><code>
cd docs
doxygen config
</pre></code>
The generated HTML documentation can be viewed by pointing a HTML browser to the index.html file in the
html directory. The generated LATEX documentation must first be compiled by a LATEX compiler.
## Tests

All functions from Faddeeva package are covered with tests. 
Test complex arguments were borrowed from original source Faddeeva.cc
file except Inf and NaN according they are not supported in Fortran 90. Reference values were computed by WolframAlpha.

### Launching the tests
<pre><code>
cmake .
make
ctest
</pre></code>

## About

This project was done with the assistance of [Konstantin 
Ladutenko](https://github.com/kostyfisik). 