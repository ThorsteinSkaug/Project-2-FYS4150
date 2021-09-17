all: Compile Link Execute

Compile:
	g++ -c -std=c++11 project2_3.cpp -larmadillo

Link:
	g++ -std=c++11 project2_3.o -o project2_3.exe  -larmadillo 

Execute:
	./project2_3.exe
