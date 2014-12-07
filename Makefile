CXX    = g++
CFLAGS = -Wall -std=c++0x -g

all : libpngwriter testOctree testMatrix MainYee Main

libpngwriter    :
	cd pngwriter; make all

testOctree  : testOctree.cpp
	g++ $@.cpp -o $@.exe $(CFLAGS)

testMatrix  : testMatrix.cpp
	g++ $@.cpp -o $@.exe $(CFLAGS)

MainYee : libpngwriter MainYee.cpp
	mkdir -p output
	mpic++ $@.cpp -o $@.exe $(CFLAGS) -DNO_FREETYPE -I ./pngwriter/ -L ./pngwriter/ -lz -lpngwriter -lpng

Main    : Main.cpp
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
	ffmpeg -i output/Ex_%05d.png -c:v libvpx -crf 4 Ex.webm
	ffmpeg -i output/Ey_%05d.png -c:v libvpx -crf 4 Ey.webm
	ffmpeg -i output/Ez_%05d.png -c:v libvpx -crf 4 Ez.webm
	ffmpeg -i output/Hx_%05d.png -c:v libvpx -crf 4 Hx.webm
	ffmpeg -i output/Hy_%05d.png -c:v libvpx -crf 4 Hy.webm
	ffmpeg -i output/Hz_%05d.png -c:v libvpx -crf 4 Hz.webm