CXX    = g++
CFLAGS = -Wall -std=c++0x -g

all : libpngwriter testOctree testMatrix MainYee Main

libpngwriter    :
	cd pngwriter; make all

testOctree  :
	g++ $@.cpp -o $@.exe $(CFLAGS)

testMatrix  :
	g++ $@.cpp -o $@.exe $(CFLAGS)

MainYee : libpngwriter
	mkdir -p output
	mpic++ $@.cpp -o $@.exe $(CFLAGS) -DNO_FREETYPE -I ./pngwriter/ -L ./pngwriter/ -lz -lpngwriter -lpng

Main    :
	mpic++ $@.cpp -o $@.exe $(CFLAGS)

clean   :
	cd pngwriter; make clean
	rm -f testOctree.exe
	rm -f testOctree
	rm -f testMatrix.exe
	rm -f testMatrix
	rm -f MainYee.exe
	rm -f MainYee
	rm -f Main.exe
	rm -f Main
	rm -f output/*
	rm -f Ex.webm
	rm -f Ey.webm
	rm -f Ez.webm
	rm -f Hx.webm
	rm -f Hy.webm
	rm -f Hz.webm

animation:
	ffmpeg -i output/Ex_%05d.png Ex.webm
	ffmpeg -i output/Ey_%05d.png Ey.webm
	ffmpeg -i output/Ez_%05d.png Ez.webm
	ffmpeg -i output/Hx_%05d.png Hx.webm
	ffmpeg -i output/Hy_%05d.png Hy.webm
	ffmpeg -i output/Hz_%05d.png Hz.webm