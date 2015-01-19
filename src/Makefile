CXX     = g++
INCLUDE = -I . -I ./libs
CFLAGS  = $(INCLUDE) -Wall -std=c++0x -g -O0

all : libpngwriter testOctree testMatrix MainYee Main

libpngwriter:
	cd pngwriter; make all

TVector : libs/math/TVector.tpp libs/math/TVector.h

TBaseMatrix : libs/math/TBaseMatrix.h

Octree : libs/octree/OctreeNode.h libs/octree/OctreeNode.tpp libs/octree/OctreeToSvg.h libs/octree/OctreeToSvg.tpp

testOctree  : examples/testOctree.cpp TVector Octree
	rm -f $@.exe
	g++ $< -o $@.exe $(CFLAGS)

testOctree2 : examples/testOctree2.cpp TVector Octree
	rm -f $@.exe
	g++ $< -o $@.exe $(CFLAGS)

testMatrix  : examples/testMatrix.cpp TBaseMatrix TVector
	rm -f $@.exe
	g++ $< -o $@.exe $(CFLAGS)

MainYee : yeesolver/MainYee.cpp yeesolver/YeeCell.h libpngwriter
	mkdir -p output
	rm -f $@.exe
	mpic++ $< -o $@.exe $(CFLAGS) -DNO_FREETYPE -I ./pngwriter -L ./pngwriter -lz -lpngwriter -lpng

MainYeeOctree : yeesolver/MainYeeOctree.cpp yeesolver/YeeCell.h Octree libpngwriter
	mkdir -p output
	rm -f $@.exe
	mpic++ $< -o $@.exe $(CFLAGS) -DNO_FREETYPE -I ./pngwriter -L ./pngwriter -lz -lpngwriter -lpng

testSimBox : examples/testSimBox/Main.cpp examples/testSimBox/Communicator.h examples/testSimBox/Cell.h TVector TBaseMatrix
	rm -f $@.exe
	mpic++ $< -o $@.exe $(CFLAGS)

clean:
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