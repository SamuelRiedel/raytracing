untitled.ppm: run
	feh $@

run: raytracer
	./$<

%: %.cpp
	c++ -o $@ -O3 -Wall $<

format:
	clang-format -i *.cpp
