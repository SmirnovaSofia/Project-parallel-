main: main.o
	g++ -o main main.cpp

main.o: main.cpp
	g++ -c main.cpp -o main.o 

float: main

double: main
	g++ -DAHAHA main.cpp
	mv a.out main

clean:
	rm *.o main

.PHONY: clean