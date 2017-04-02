# CSWITCHES is a list of all switches passed to the C compiler.  I strongly
#   recommend using the best level of optimization.  I also strongly
#   recommend timing each level of optimization to see which is the
#   best.  For instance, when I had a DEC Alpha using DEC's optimizing
#   compiler, the -O2 switch generated a notably faster version of Triangle
#   than the -O3 switch.  Go figure.
#
# By default, Triangle and Show Me use double precision floating point
#   numbers.  If you prefer single precision, use the -DSINGLE switch.
#   Double precision uses more memory, but improves the resolution of
#   the meshes you can generate with Triangle.  It also reduces the
#   likelihood of a floating exception due to overflow.  Also, it is
#   much faster than single precision on many architectures.  I recommend
#   double precision unless you want to generate a mesh for which you do
#   not have enough memory to use double precision.
#
# If yours is not a Unix system, use the -DNO_TIMER switch to eliminate the
#   Unix-specific timer code.  Also, don't try to compile Show Me; it only
#   works with X Windows.
#
# To get the exact arithmetic to work right on an Intel processor, use the
#   -DCPU86 switch with Microsoft C, or the -DLINUX switch with gcc running
#   on Linux.  The floating-point arithmetic might not be robust otherwise.
#   Please see http://www.cs.cmu.edu/~quake/robust.pc.html for details.
#
# If you are modifying Triangle, I recommend using the -DSELF_CHECK switch
#   while you are debugging.  Defining the SELF_CHECK symbol causes
#   Triangle to include self-checking code.  Triangle will execute more
#   slowly, however, so be sure to remove this switch before compiling a
#   production version.
#
# If the size of the Triangle binary is important to you, you may wish to
#   generate a reduced version of Triangle.  The -DREDUCED switch gets rid
#   of all features that are primarily of research interest.  Specifically,
#   defining the REDUCED symbol eliminates the -i, -F, -s, and -C switches.
#   The -DCDT_ONLY switch gets rid of all meshing algorithms above and beyond
#   constrained Delaunay triangulation.  Specifically, defining the CDT_ONLY
#   symbol eliminates the -r, -q, -a, -u, -D, -S, and -s switches.  The
#   REDUCED and CDT_ONLY symbols may be particularly attractive when Triangle
#   is called by another program that does not need all of Triangle's
#   features; in this case, these switches should appear as part of
#   "TRILIBDEFS" below.
#
# On some systems, you may need to include -I/usr/local/include and/or
#   -L/usr/local/lib in the compiler options to ensure that the X include
#   files and libraries that Show Me needs are found.  If you get errors
#   like "Can't find include file X11/Xlib.h", you need the former switch.
#   Try compiling without them first; add them if that fails.
#
# An example CSWITCHES line is:
#
#   CSWITCHES = -O -DNO_TIMER -DLINUX -I/usr/X11R6/include -L/usr/X11R6/lib

# TRILIBDEFS is a list of definitions used to compile an object code version
#   of Triangle (triangle.o) to be called by another program.  The file
#   "triangle.h" contains detailed information on how to call triangle.o.
#
# The -DTRILIBRARY should always be used when compiling Triangle into an
#   object file.
#
# An example TRILIBDEFS line is:
#
#   TRILIBDEFS = -DTRILIBRARY -DREDUCED -DCDT_ONLY
#=========================================================

SRC = $(PWD)/src/
BIN = $(PWD)/bin/
DEBUG = $(PWD)/debug/
INC = $(PWD)/include/
LIB = $(PWD)/lib/

CC = gcc

CSWITCHES = -DLINUX -DSELF_CHECK -I/usr/X11R6/include -L/usr/X11R6/lib
TRILIBDEFS = -DTRILIBRARY

DBG = -ggdb
OPT = -O -DNDEBUG

RM = /bin/rm -f



# +-----------------------------+
# | Release configuration |
# +-----------------------------+
shared: $(LIB)triangle.o $(LIB)allocate.o $(LIB)tesselate.o $(SRC)tesselate.c
				@echo "Compiling   $@"
				@$(CC) -shared $(OPT) $(CSWITCHES) -I$(INC) -o \
					$(LIB)libtesselate.so \
					$(LIB)tesselate.o \
					$(LIB)allocate.o \
					$(LIB)triangle.o -lm

trilib: $(BIN)triangle.o $(BIN)tesselate

$(BIN)tesselate: $(SRC)tesselate.c $(BIN)allocate.o $(BIN)triangle.o
				@echo "Compiling   $@"
				@$(CC) $(OPT) $(CSWITCHES) -I$(INC) -o \
					$(BIN)tesselate \
					$(SRC)tesselate.c \
					$(BIN)allocate.o \
					$(BIN)triangle.o -lm

$(LIB)tesselate.o: $(SRC)tesselate.c $(LIB)allocate.o $(LIB)triangle.o
				@echo "Compiling   $@"
				@$(CC) -fpic $(OPT) $(CSWITCHES) -I$(INC) -c -o \
					$(LIB)tesselate.o \
					$(SRC)tesselate.c


$(BIN)triangle.o: $(SRC)triangle.c $(INC)triangle.h
				@echo "Compiling   $@"
				@$(CC) $(OPT) $(CSWITCHES) $(TRILIBDEFS) -I$(INC) -c -o \
					$(BIN)triangle.o \
					$(SRC)triangle.c

$(LIB)triangle.o: $(SRC)triangle.c $(INC)triangle.h
				@echo "Compiling   $@"
				@$(CC) -fpic $(OPT) $(CSWITCHES) $(TRILIBDEFS) -I$(INC) -c -o \
					$(LIB)triangle.o \
					$(SRC)triangle.c

$(BIN)allocate.o: $(SRC)allocate.c $(INC)allocate.h
				@echo "Compiling   $@"
				@$(CC) $(OPT) $(CSWITCHES) $(TRILIBDEFS) -I$(INC) -c -o \
					$(BIN)allocate.o \
					$(SRC)allocate.c

$(LIB)allocate.o: $(SRC)allocate.c $(INC)allocate.h
				@echo "Compiling   $@"
				@$(CC) -fpic $(OPT) $(CSWITCHES) $(TRILIBDEFS) -I$(INC) -c -o \
					$(LIB)allocate.o \
					$(SRC)allocate.c

# +---------------------------+
# | Debug configuration |
# +---------------------------+
debug: $(DEBUG)triangle.o $(DEBUG)tesselate

$(DEBUG)tesselate: $(SRC)tesselate.c $(DEBUG)triangle.o
				@echo "Compiling debug in configuration  $@"
				@$(CC) $(DBG) $(CSWITCHES) -I$(INC) -o \
					$(DEBUG)tesselate \
					$(SRC)tesselate.c \
					$(DEBUG)triangle.o -lm

$(DEBUG)triangle.o: $(SRC)triangle.c $(INC)triangle.h
				@echo "Compiling debug in configuration  $@"
				@$(CC) $(DBG) $(CSWITCHES) $(TRILIBDEFS) -I$(INC) -c -o \
					$(DEBUG)triangle.o \
					$(SRC)triangle.c

clean:
	@echo "Clearing directories..."
	@$(RM) $(LIB)*.* $(BIN)*.* $(DEBUG)*.*
