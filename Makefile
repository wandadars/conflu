### Makefile for CONFLU ########################################################

### Macros #####################################################################
SHELL= /bin/sh

BINDIR= $(HOME)/bin
DISTPREFIX= conflu
DISTNAME= $(subst *,$(shell date '+%m%d%y%H%M'),$(DISTPREFIX)_*.tar)
EXEC= conflu
LDFLAGS= 
LIBS= 
RIFFRAFF= *.bak *.bck *.ckp
RM= rm -f
SYSTEM= $(shell uname -s)

### Platform-specific stuff ####################################################
MACHINE= $(shell uname -n)


### Streamline Numerics Machine 0 ###
ifeq ($(MACHINE),SN-HPC-0)
  FC= mpif90
  FFLAGS= -g 
  #FFLAGS= -O 
  INCLUDES= 
  LDFLAGS+= 
  LIBS+= 
  MODEXT= mod  
endif 

### gator machine at uf ###
ifeq ($(MACHINE),gator1.ufhpc)
  FC= mpif90
  FFLAGS= -g -qsuffix=f=f90
  #FFLAGS= -O 
  INCLUDES= 
  LDFLAGS+= 
  LIBS+= 
  MODEXT= mod  
endif 

### gator machine at uf ###
ifeq ($(MACHINE),gator2.ufhpc)
  FC= mpif90
  FFLAGS= -g -qsuffix=f=f90
  #FFLAGS= -O 
  INCLUDES= 
  LDFLAGS+= 
  LIBS+= 
  MODEXT= mod  
endif 

### submit machine at uf ###
ifeq ($(MACHINE),submit2.ufhpc)
  FC= mpif90
  FFLAGS= -g -qsuffix=f=f90
  #FFLAGS= -O 
  INCLUDES= 
  LDFLAGS+= 
  LIBS+= 
  MODEXT= mod  
endif 

ifeq ($(MACHINE),submit.ufhpc)
  FC= mpif90
  FFLAGS= -g -qsuffix=f=f90
  #FFLAGS= -O 
  INCLUDES= 
  LDFLAGS+= 
  LIBS+= 
  MODEXT= mod  
endif 

### My mac workstation, intel compiler ###
ifeq ($(MACHINE),crocco-mae-ufl-edu.local)
  FC= ifort
  FFLAGS= -g n free
  #FFLAGS= -g -e n -f free -YEXT_NAMES=LCS -YEXT_SFX=_
  #FFLAGS= -O -f free
  INCLUDES= 
  LDFLAGS+= 
  LIBS+= 
  MODEXT= mod  
endif 

### My workstation, Absoft compiler ###
ifeq ($(MACHINE),titov.csar.uiuc.edu)
  FC= f90 
  FFLAGS= -g -e n -f free
  #FFLAGS= -g -e n -f free -YEXT_NAMES=LCS -YEXT_SFX=_
  #FFLAGS= -O -f free
  INCLUDES= 
  LDFLAGS+= 
  LIBS+= 
  MODEXT= mod  
endif 

### My laptop, Absoft compiler ###
ifeq ($(MACHINE),popovich.csar.uiuc.edu)
  FC= f90 
  FFLAGS= -g -en -f free
  #FFLAGS= -g -e n -f free -YEXT_NAMES=LCS -YEXT_SFX=_
  #FFLAGS= -O -f free
  INCLUDES= 
  LDFLAGS+= 
  LIBS+= 
  MODEXT= mod  
endif 

### Turing at UIUC CSE, Absoft compiler ###
ifeq ($(MACHINE),turing.cs.uiuc.edu)
  FC= f90
  #FFLAGS= -g -e n -f free
  #FFLAGS= -g -e n -f free -YEXT_NAMES=LCS -YEXT_SFX=_
  FFLAGS= -O -f free 
  INCLUDES= 
  LDFLAGS+= 
  LIBS+= 
  MODEXT= mod
endif 

### Posic at UIUC NCSA ###
ifeq ($(MACHINE),ntsc1169.ncsa.uiuc.edu)
  FC= pgf90 
  FFLAGS= -g -Mbounds -Mchkptr -Mfreeform -Mstandard
  INCLUDES= 
  LDFLAGS+= 
  LIBS+= 
endif 

### New turing ###
ifeq ($(findstring turing,$(MACHINE)),turing)
  FC= xlf90 
  FFLAGS= -g -qsuffix=f=f90
  #FFLAGS= -O 
  INCLUDES= 
  LDFLAGS+= 
  LIBS+= 
  MODEXT= mod  
endif 

### Pattern rules ##############################################################
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

### Source and object section ##################################################
SRC= modGlobals.f90 modError.f90 modGrid.f90 modHashTable.f90 modSortSearch.f90\
     CENTAUR2Generic.f90 conflu.f90 COBALT2CENTAUR.f90\
     distortVertices.f90 generic2CENTAUR.f90 ijk2l.f90 PLOT3D2CENTAUR.f90\
     PLOT3D2Generic.f90 readGridCarpentier.f90 readGridCENTAUR.f90\
     readGridCOBALT.f90 readGridGAMBIT2d.f90 readGridLeyland.f90\
     readGridPLOT3D.f90 readGridPLOT3DBinary.f90 readGridPlourde.f90\
     readGridSCREAM.f90 readGridTRIANGLE.f90 writeGridCENTAUR.f90\
     writeGridCENTAURASCII.f90 writeGridPLOT3DBinary.f90
OBJ= $(SRC:.f90=.o)

### Target section #############################################################
$(EXEC): $(OBJ)
	$(FC) $(OBJ) $(LDFLAGS) -o $@ $(LIBS)

.PHONY: clean install
clean:
	$(RM) $(OBJ) *.$(MODEXT) $(EXEC) $(RIFFRAFF) 

clear:
	$(RM) $(OBJ) *.$(MODEXT) $(RIFFRAFF)         
	
install: $(EXEC)
	cp $(EXEC) $(BINDIR)/$(EXEC)

dist:	
	$(MAKE) clean
	tar -cvf $(DISTNAME) * 
	gzip $(DISTNAME)

### Dependencies section #######################################################

