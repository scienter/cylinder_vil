EXEC = show
CC = h5pcc
#CC = mpicc
OBJS = main.o parameterSetting.o findparam.o boundary.o saveFile.o fieldSolve.o loadLaser.o loadPlasma.o pml.o utility.o fieldShareX.o saveFieldHDF.o saveDensityHDF.o saveParticleHDF.o saveDumpHDF.o solveCharge.o clean.o materials.o ionization.o interpolation.o particlePush.o rearrangeParticles.o removeEdge.o movingDomain.o particleShareX.o filter.o updateCurrent.o restoreDumpHDF.o redist.o calConservation.o 


#-------- for Beam server ----------# 
#CFLAGS = -I/opt/gsl/2.6/include -I/opt/hdf5/1.8.20/include
#LDFLAGS = -L/opt/gsl/2.6/lib -L/opt/hdf5/1.8.20/lib

#-------- for PAL ----------#
CFLAGS = -I/home/scienter/gsl2.6/include -I/usr/lib/x86_64-linux-gnu/hdf5/openmpi/include
LDFLAGS = -L/home/scienter/gsl2.6/lib -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi

#-------- for home computer ----------#
#CFLAGS= -I/opt/gsl/2.6/include
#LDFLAGS= -L/opt/gsl/2.6/lib





INCL = constants.h laser.h mesh.h particle.h plasma.h
LIBS = -lm -lhdf5 -lgsl -lgslcblas




$(EXEC):$(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) $(LIBS) -o $(EXEC)
$(OBJS):$(INCL)
clean:
	@rm -f *.o *~ $(EXEC)
