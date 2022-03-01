
all:
	g++ -O3 -std=c++17  -g -march=native md.cpp -o md

clean:
	rm -f md
