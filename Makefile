# local compile
dsim: dsim.cpp mvec.h
	g++ dsim.cpp -o dsim -lSDL -lSDL_ttf -O3 -I/usr/include/SDL/ -lSDL_gfx -std=c++17 -lgd -fconcepts

# compile for the web with Emscripten
dsim.html: dsim.cpp mvec.h
	em++ -s TOTAL_MEMORY=65536000 dsim.cpp -o dsim.html -lSDL -O3 -I/usr/include/SDL/ -std=c++1z

