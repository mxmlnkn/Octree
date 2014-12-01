########################## PNGwriter #########################################
#
#   Website: Main:             http://pngwriter.sourceforge.net/
#            GitHub.com:       https://github.com/pngwriter/pngwriter
#            Sourceforge.net:  http://sourceforge.net/projects/pngwriter/
#
#
#    Author:                    Paul Blackburn https://github.com/individual61
#                               Axel Huebl https://github.com/ax3l
#
#    Email:                     individual61@users.sourceforge.net
#
#    Version:                   0.5.5 (June 2014)
#
#    Description:               Library that allows plotting a 48 bit
#                               PNG image pixel by pixel, which can
#                               then be opened with a graphics program.
#
#    License:                   GNU General Public License
#                               (C) 2002-2014 Paul Blackburn
#                               (C) 2013-2014 Axel Huebl
#
##############################################################################

######################## IMPORTANT ############################################
# This makefile is meant to help beginning programmers learn about simple 
# Makefiles. The compilation uses the object file, and not the library,
# because I believe it to be a more 'hands-on' approach.
###############################################################################

CXX      = g++
CXXFLAGS = -O3 -Wall -DNO_FREETYPE
INC      = -I./
LIBS     = -L./ -lz -lpngwriter -lpng

OBJECTS  = pngwriter.o
EXAMPLES = pngtest lyapunov diamond blackwhite

all: libpngwriter.a $(EXAMPLES)

libpngwriter.a: $(OBJECTS)
    # pack .a archive (library containing object file(s))
	ar rv $@ $^
    # create and store index into that archive
	ranlib $@
    # $@ and $^ are left and right side of :

pngwriter.o: pngwriter.cc pngwriter.h
	$(CXX) $(CXXFLAGS) $(INC) -g -c -o pngwriter.o pngwriter.cc

pngtest: libpngwriter.a pngtest.cc
	$(CXX) $(CXXFLAGS) $(INC) pngtest.cc -o pngtest $(LIBS) 

lyapunov: libpngwriter.a lyapunov.cc
	$(CXX) $(CXXFLAGS) $(INC) lyapunov.cc -o lyapunov $(LIBS) 

diamond: libpngwriter.a diamond.cc
	$(CXX) $(CXXFLAGS) $(INC) diamond.cc -o diamond $(LIBS) 

blackwhite: libpngwriter.a blackwhite.cc
	$(CXX) $(CXXFLAGS) $(INC) blackwhite.cc -o blackwhite $(LIBS) 

clean	:    
	rm -f $(EXAMPLES) *~ arcoiris.png copiaburro.png one.png two.png
	rm -f out.png triangles.png
	rm -f $(OBJECTS) libpngwriter.a