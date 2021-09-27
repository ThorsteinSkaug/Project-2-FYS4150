3: Compile3 Link3 Execute3
4: Compile4 Link4 Execute4

Compile3:
	g++ -c -std=c++11 project2_problem3.cpp -larmadillo
Compile4:
	g++ -c -std=c++11 project2_problem4567.cpp -larmadillo

Link3:
	g++ -std=c++11 project2_problem3.o -o project2_problem3.exe  -larmadillo
Link4:
	g++ -std=c++11 project2_problem4567.o -o project2_problem4567.exe  -larmadillo

Execute3:
	./project2_problem3.exe
Execute4:
	./project2_problem4567.exe
