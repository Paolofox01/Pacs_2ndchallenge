### Second challenge of the APSC course 

## Content

`algebra.hpp` -> definition of a class that stores sparse matrixes in dynamic or compressed form, either with Row-major or Column-major ordering

`Matrix_Impl.hpp` -> implementaation of the members of algebra.hpp

`main.cpp` -> reads a matrix in .mtx format and multiplies it with a vector of the right dimension fillend with ones using both the dynamic and compressed form, printing the time needed for each computation

`matrix.json` -> contains the name of the .mtx file to be read

`Doxifile` -> contains the instructions needed to create the documentation using Doxygen

`Makefile` -> contains the intructions needed to compile the code

## How to compile and run the code

modify the PACS_ROOT path in the Makefile and then run 

````
make
````

to run the code using the matrix written in the json file just run

````
./main
````

if you want you can also use another .mtx file without modyfing the .json by running

````
./main filename.mtx
````

## How to create the documentiation with Doxygen

run 

````
make doc
````

then open the index.html file in /doc/html to view the documentation 


