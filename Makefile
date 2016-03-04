# define variables
VPATH   = ./
HDRDIR  = /home/li12242/Documents/workspace/MIDG/include
LIBDIR	= /home/li12242/Documents/workspace/MIDG/lib
MPIDIR = /home/li12242/Public/Tom/mpich2-gnu-3.2/

# adjust this for your system

# set options for this machine
# specify which compilers to use for c, fortran and linking
CC	= /home/li12242/Public/Tom/mpich2-gnu-3.2/bin/mpicc
LD	= /home/li12242/Public/Tom/mpich2-gnu-3.2/bin/mpicxx

# compiler flags to be used (set to compile with debugging on)
CFLAGS = -Dp_N=4 -DNDG2d -I/opt/local/include -I/usr/include/malloc -I$(HDRDIR) -O3

# link flags to be used 
LDFLAGS	= -I$(HDRDIR) -L. -L$(LIBDIR) -O3 -I$(MPIDIR)/include -L$(MPIDIR)/lib

# libraries to be linked in
LIBS	=  -lparmetis -lmetis  -lm 

# types of files we are going to construct rules for
.SUFFIXES: .c .f .cu

# rule for .c files
.c.o:
	$(CC) $(CFLAGS) -o $*.o -c $*.c

# list of objects to be compiled
OBJS    = \
	src/Mesh2d.o\
	src/Utils.o\
	src/LoadBalance2d.o\
	src/FacePair2d.o\
	src/ParallelPairs.o\
	src/BuildMaps2d.o\
	src/StartUp2d.o\
	src/MaxwellsRun2d.o\
	src/MaxwellsMPI2d.o\
	src/MaxwellsDriver2d.o\
	src/MaxwellsRHS2d.o\
	src/InitCPU2d.o

all:$(OBJS)
	$(LD)  $(LDFLAGS) -o MaxwellsCPU2d $(OBJS) $(LIBS)


# what to do if user types "make clean"
clean :
	rm -r $(OBJS)
	rm MaxwellsCPU2d