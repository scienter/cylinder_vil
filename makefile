EXEC = show
#CC = mpicc
CC = /opt/openmpi/bin/mpicc

OBJS = main.o parameterSetting.o findparam.o boundary.o saveFile.o fieldSolve.o loadLaser.o loadPlasma.o pml.o utility.o fieldShareX.o saveFieldHDF.o saveDensityHDF.o saveParticleHDF.o saveDumpHDF.o solveCharge.o clean.o materials.o ionization.o interpolation.o particlePush.o rearrangeParticles.o removeEdge.o movingDomain.o particleShareX.o filter.o updateCurrent.o restoreDumpHDF.o redist.o 

CFLAGS= -I/opt/openmpi/include -I/usr/include/gsl

#CFLAGS= -I/home/scienter/gsl/include
#LDFLAGS= -L/home/scienter/gsl/lib
#INCDIR = /home/scienter/gsl/include /home/scienter/gsl/lib
#INCL = constants.h laser.h mesh.h particle.h plasma.h
#LIBS = -lm -lhdf5 -lgsl -lgslcblas
LIBS  = -Wl,-rpath,/opt/intel/lib/intel64 -Wl,-rpath,/opt/openmpi/lib -L/opt/hdf5/lib -Wl,-rpath,/opt/zlib/lib -lpython2.7 -lgslcblas -lgsl -lmpi -lhdf5 -lpthread -lm -lz


$(EXEC):$(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) $(LIBS) -o $(EXEC)
$(OBJS):$(INCL)
clean:
	@rm -f *.o *~ $(EXEC)
