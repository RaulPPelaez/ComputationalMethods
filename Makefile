
all:
	g++ -O3 -std=c++17  -g -march=native md.cpp -o md -ltbb

clean:
	rm -f md
