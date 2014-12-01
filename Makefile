CXX    = g++
CFLAGS = -Wall -std=c++0x

all : libpngwriter testOctree testMatrix MainYee Main

libpngwriter    :
	cd pngwriter; make all

testOctree  :
	g++ $@.cpp -o $@.exe $(CFLAGS)

testMatrix  :
	g++ $@.cpp -o $@.exe $(CFLAGS)

MainYee : libpngwriter
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
