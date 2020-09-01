rm bin/a.out
g++ --std=c++11 -O2 *.cpp -o bin/a.out -I/usr/local/include -L/usr/local/lib -lgsl -lgslcblas -lm
cd bin
./a.out>../std.txt&
