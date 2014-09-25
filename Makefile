MPI_CC=mpicxx
MPI_FLAGS=-O3

MPI_LD_LIB=
MPI_INCL=
MPI_LIBS=

CC=g++
CC_FLAGS=-O3

SOURCE=lbm_r2.cxx
TARGET=pois3D


$(TARGET):  $(SOURCE) lbm_lib.o vtk_lib.o
	$(MPI_CC) -o $(TARGET) $(SOURCE) $(MPI_FLAGS) lbm_lib.o vtk_lib.o

lbm_lib.o:  lbm_lib.cxx
	$(CC) -c lbm_lib.cxx $(CC_FLAGS) 

vtk_lib.o:  vtk_lib.cxx
	$(CC) -c vtk_lib.cxx

clean:
	rm *.o $(TARGET)