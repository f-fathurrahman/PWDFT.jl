MAKE = make
AR = ar

CC = gcc
CC_INC = 
CC_OPTS = -Wall -O3 -fPIC -cpp

LIBS = -lm

C_SRC = \
arithmetic.c  determination.c  mathfunc.c    primitive.c         spacegroup.c    spin.c \
cell.c        hall_symbol.c    niggli.c      refinement.c        spg_database.c  symmetry.c \
debug.c       kgrid.c          overlap.c     sitesym_database.c  spglib.c        test.c \
delaunay.c    kpoint.c         pointgroup.c  site_symmetry.c     spglib_f.c      tetrahedron_method.c

OBJ = $(C_SRC:.c=.o) 

%.o : %.c
	$(CC) $(CC_INC) $(CC_OPTS) -c -o $(*F).o $<

# Targets
lib: $(OBJ)
	ar rcs libspglib.so *.o

clean:
	rm -rf *.o *.mod libspglib.so *.x


