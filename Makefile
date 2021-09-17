all: Compile Link Execute

Compile:
	g++ -c -std=c++11 project2_4.cpp -larmadillo

Link:
	g++ -std=c++11 project2_4.o -o project2_4.exe  -larmadillo

Execute:
	./project2_4.exe
