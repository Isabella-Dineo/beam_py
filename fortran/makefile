######################################################################
#
#
#########################################################################
# macro definitions
FFLAGS = -fno-backslash 
# UBUNTU VERSION 
FORT = gfortran
PGPLOT = -L/usr/lib -lpgplot -L/usr/lib/x86_64-linux-gnu -lpng -lz
LDFLAGS=-L/usr/lib/x86_64-linux-gnu -lX11
CFITSIO=-L/usr/lib/x86_64-linux-gnu -lcfitsio
# MAC VERSION 
#FORT = gfortran
#PGPLOT = -L/opt/local/lib -lpgplot -L/sw/lib -lpng -lz
#LDFLAGS=-L/usr/X11R6/lib -lX11 -L/sw/lib -laquaterm -Wl,-framework -Wl,Foundation
#CFITSIO=-L/opt/local/lib -lcfitsio
#
OBJ = beam_sub.o plotbeam.o plotlos.o math.o mapphi.o plotprof.o patchbeam.o \
	ran1.o plotpa.o getparams.o getlos.o addnoise.o deriv.o classify.o \
	profparams.o rvm.o scatter.o beamfits.o

#
# target definition list
#
all     : beam genpop trueStokes gendunc

beam:	beam.o $(OBJ)
	@echo 'Linking...'
	@$(FORT) $(FFLAGS) beam.o $(OBJ) -o beam $(CFITSIO) $(LDFLAGS) $(PGPLOT)
	@echo 'Done. Created: beam'
beam2fits:	beam2fits.o beamfits.o
	@echo 'Linking...'
	@$(FORT) $(FFLAGS) beam2fits.o beamfits.o -o beam2fits $(CFITSIO) $(LDFLAGS) 
	@echo 'Done. Created: beamfits'
genpop: genpop.o dohist.o $(OBJ)
	@$(FORT) $(FFLAGS) genpop.o dohist.o $(OBJ) -o genpop $(LDFLAGS)  $(PGPLOT)
gendunc: gendunc.o dohist.o $(OBJ)
	@$(FORT) $(FFLAGS) gendunc.o dohist.o $(OBJ) -o gendunc $(LDFLAGS)  $(PGPLOT)
trueStokes: trueStokes.o $(OBJ)
	@$(FORT) $(FFLAGS) trueStokes.o $(OBJ) -o trSt $(LDFLAGS)  $(PGPLOT)

beamfits.o:		beamfits.f
	$(FORT)  $(FFLAGS) -c beamfits.f
beam2fits.o:		beam2fits.f
	$(FORT)  $(FFLAGS) -c beam2fits.f
beam.o:		beam.f
	$(FORT)  $(FFLAGS) -c beam.f
beam_sub.o:		beam_sub.f
	$(FORT)  $(FFLAGS) -c beam_sub.f
genpop.o:		genpop.f
	$(FORT)  $(FFLAGS) -c genpop.f
gendunc.o:		gendunc.f
	$(FORT)  $(FFLAGS) -c gendunc.f
trueStokes.o:		trueStokes.f
	$(FORT)  $(FFLAGS) -c trueStokes.f
plotbeam.o:	plotbeam.f
	$(FORT)  $(FFLAGS) -c plotbeam.f
plotlos.o:	plotlos.f
	$(FORT)  $(FFLAGS) -c plotlos.f
math.o:	math.f
	$(FORT)  $(FFLAGS) -c math.f
mapphi.o:	mapphi.f
	$(FORT)  $(FFLAGS) -c mapphi.f
plotprof.o:	plotprof.f
	$(FORT)  $(FFLAGS) -c plotprof.f
patchbeam.o:	patchbeam.f
	$(FORT)  $(FFLAGS) -c patchbeam.f
ran1.o:	ran1.f
	$(FORT)  $(FFLAGS) -c ran1.f
plotpa.o:	plotpa.f
	$(FORT)  $(FFLAGS) -c plotpa.f
getparams.o:	getparams.f
	$(FORT)  $(FFLAGS) -c getparams.f
getlos.o:	getlos.f
	$(FORT)  $(FFLAGS) -c getlos.f
addnoise.o:	addnoise.f
	$(FORT)  $(FFLAGS) -c addnoise.f
scatter.o:	scatter.f
	$(FORT)  $(FFLAGS) -c scatter.f
deriv.o:	deriv.f
	$(FORT)  $(FFLAGS) -c deriv.f
classify.o:	classify.f
	$(FORT)  $(FFLAGS) -c classify.f
dohist.o:	dohist.f
	$(FORT)  $(FFLAGS) -c dohist.f
profparams.o:	profparams.f
	$(FORT)  $(FFLAGS) -c profparams.f
rvm.o:  rvm.f
	$(FORT)  $(FFLAGS) -c rvm.f
clean:	
	rm *.o *~
