untitled.ppm: run
	feh $@

run: raytracer
	./$<

%: %.cpp
	c++ -o $@ -std=c++11 -O3 -Wall $<

%_gcc: %.cpp
	gcc -o $@ -x c -std=c99 -lm -O3 -Wall $<

format:
	clang-format -i *.cpp
