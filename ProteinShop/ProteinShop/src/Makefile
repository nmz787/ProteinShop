######################################################################
# makefile boilerplate:
######################################################################

# What operating system type are we building on?
# (Currently supported options are Linux and IRIX64 and Darwin)
SYSTEM = $(shell uname)

# What directory are include files in?
INCDIR = include

# What directory are library files in?
ifeq ($(SYSTEM),IRIX64)
  # Use n32 library files with CC:
  LIBDIR = lib32
else
  # Use standard library directory for other OSes:
  LIBDIR = lib
endif

# What directory are executables/build tools in?
BINDIR = bin

# Where is the Fltk installation located?
# The compiler expects include files under $(FLTK_BASE)/$(INCDIR),
# libraries under $(FLTK_BASE)/$(LIBDIR), and any needed build tools
# under $(FLTK_BASE)/$(BINDIR). If you have a different setup,
# you need to hand-edit the compiler flags below.
ifeq ($(SYSTEM),Darwin)
FLTK_BASE = /usr/local
else
FLTK_BASE = /usr
endif

# Which Fltk libraries to link against?
# On some installations, Fltk OpenGL functions are in a separate
# library, libfltk_gl. If the linker complains about undefined
# symbols, try adding -lfltk_gl to the line below before -lfltk.
FLTK_LIBS = -L$(FLTK_BASE)/$(LIBDIR) -lfltk_gl -lfltk

# What include directories do we need?
CIFLAGS = -I.

# What additional compiler flags do we want?
#CFLAGS = -g0 -O3 -DNDEBUG
#CFLAGS = -g -O0  -Wall
CFLAGS = -g

######################################################################
# System-specific settings:
######################################################################
# Add -DFLTK_NEW -DFLTK_DOUBLE at the end of $CSYSFLAGS 
# if your FLTK version is greater than 1.0.11

ifeq ($(SYSTEM),Darwin)
  #
  # For compilation under Mac OSX:
  #
  
  # Compiler version and flags
  CCOMP = g++
  CSYSFLAGS = -D__APPLE__ -D__LITTLE_ENDIAN=4321 -D__BIG_ENDIAN=1234 -DFLTK_NEW -DFLTK_DOUBLE
  CMAKEEXEFLAGS = 

  # More library paths:
  GL_BASE = /usr
  X11_BASE = /usr/X11R6
  
  # Libraries to link against:
  MAC_LIBS =  -framework AGL -framework OpenGL -framework Carbon 
  GL_LIBS = $(MAC_LIBS) -L$(GL_BASE)/$(LIBDIR) -lGL -L$(X11_BASE)/$(LIBDIR) -lX11
  # Possible steps after linking...
  POSTBUILD	= /Developer/Tools/Rez -t APPL -o
endif

ifeq ($(SYSTEM),Linux)
  #
  # For compilation under Linux:
  #
  
  # Compiler version and flags
  CCOMP = g++
  CSYSFLAGS = -D__LINUX__ -D__LITTLE_ENDIAN=1234 -D__BIG_ENDIAN=4321 -DFLTK_NEW -DFLTK_DOUBLE -DGL_GLEXT_PROTOTYPES
  CMAKEEXEFLAGS = -rdynamic
  
  # Protect against g++ version 4.1.0 and up:
  ifeq ($(shell expr `g++ -dumpversion` ">=" "4.1.0"),1)
    CSYSFLAGS += -ffriend-injection
  endif
  
  # More library paths:
  GL_BASE = /usr
  X11_BASE = /usr/X11R6
  
  # Libraries to link against:
  GL_LIBS = -L$(GL_BASE)/$(LIBDIR) -lGL
endif

ifeq ($(SYSTEM),IRIX64)
  #
  # For compilation under IRIX64:
  #
  
  # Compiler version and flags
  CCOMP = CC
  CSYSFLAGS = -LANG:std -woff 1355 -no_auto_include -no_prelink -ptused
  CSYSFLAGS += -D__SGI_IRIX__ -D__LITTLE_ENDIAN=4321 -D__BIG_ENDIAN=1234 -DNONSTANDARD_TEMPLATES
  CMAKEEXEFLAGS = 
  
  # More library paths:
  GL_BASE = /usr
  X11_BASE = /usr
  
  # Libraries to link against:
  GL_LIBS = -L$(GL_BASE)/$(LIBDIR) -lGL -L$(X11_BASE)/$(LIBDIR) -lX11
endif

INSTALL = mv ProteinShop ../bin; cd ../bin 

######################################################################
# Default C++ compiler command:
######################################################################

CCFLAGS = $(CSYSFLAGS) $(CIFLAGS) $(CFLAGS) -I$(FLTK_BASE)/$(INCDIR) 

%.o: %.cpp
	$(CCOMP) -c -o $@ $(CCFLAGS) $<

######################################################################
# Dependencies and targets:
######################################################################

ALL = ProteinShop
all: $(ALL)

clean:
	rm -f *.o *.so Makefile.bak core.* ; true
	rm -rf ii_files ; true
	rm -f Misc/*.o ; true
	rm -rf Misc/ii_files ; true
	rm -f Math/*.o ; true
	rm -rf Math/ii_files ; true
	rm -f Geometry/*.o ; true
	rm -rf Geometry/ii_files ; true
	rm -f $(ALL) ; true
	rm -f ProteinGui.h ProteinGui.cpp ; true


######################################################################
# Miscellaneous library classes:
######################################################################

FILEIO_H = FileIO.h

ENDIANNESS_H = Endianness.h

PRIORITYHEAP_H = PriorityHeap.h

POOLALLOCATOR_H = PoolAllocator.h

HASHTABLE_H = $(POOLALLOCATOR_H) HashTable.h

SERIALPORT_H = SerialPort.h
SerialPort.o: $(SERIALPORT_H) SerialPort.cpp

TCPSOCKET2_H = TCPSocket2.h
TCPSocket2.o: $(TCPSOCKET2_H) TCPSocket2.cpp

TIMER_H = Timer.h
Timer.o: $(TIMER_H) Timer.cpp

CONFIGURATIONFILE_H = ConfigurationFile.h
ConfigurationFile.o: $(CONFIGURATIONFILE_H) ConfigurationFile.cpp

DENSEMATRIX_H = DenseMatrix.h
DenseMatrix.o: $(DENSEMATRIX_H) DenseMatrix.cpp

UTILITY_H = Misc/Utility.h

MISC_O = SerialPort.o \
         TCPSocket2.o \
         Timer.o \
         ConfigurationFile.o \
         DenseMatrix.o

######################################################################
# Mathematics library classes:
######################################################################

MMATH_H = Math/Math.h

MCONSTANTS_H = Math/Constants.h
Math/Constants.o: $(MCONSTANTS_H) Math/Constants.cpp

MRANDOM_H = Math/Random.h
Math/Random.o: $(MRANDOM_H) Math/Random.cpp

MATH_O = Math/Constants.o \
         Math/Random.o \
         Math/Interval.o

######################################################################
# Geometry library classes:
######################################################################

GCOMPONENTARRAY_H = $(MMATH_H) Geometry/ComponentArray.h
Geometry/ComponentArray.o: $(GCOMPONENTARRAY_H) Geometry/ComponentArray.cpp

GVECTOR_H = $(MMATH_H) $(GCOMPONENTARRAY_H) Geometry/Vector.h
Geometry/Vector.o: $(GVECTOR_H) Geometry/Vector.cpp

GPOINT_H = $(MMATH_H) $(GCOMPONENTARRAY_H) $(GVECTOR_H) Geometry/Point.h
Geometry/Point.o: $(GPOINT_H) Geometry/Point.cpp

GAFFINECOMBINER_H = $(GCOMPONENTARRAY_H) $(GPOINT_H) Geometry/AffineCombiner.h

GHVECTOR_H = $(MMATH_H) $(GCOMPONENTARRAY_H) $(GVECTOR_H) $(GPOINT_H) Geometry/HVector.h
Geometry/HVector.o: $(GHVECTOR_H) Geometry/HVector.cpp

GMATRIXHELPERFUNCTIONS_H = $(UTILITY_H) Geometry/MatrixHelperFunctions.h

GMATRIX_H = $(GCOMPONENTARRAY_H) $(GVECTOR_H) $(GPOINT_H) $(GHVECTOR_H) Geometry/Matrix.h
Geometry/Matrix.o: $(MMATH_H) $(GMATRIX_H) $(GMATRIXHELPERFUNCTIONS_H) Geometry/Matrix.cpp

GROTATION_H = $(MMATH_H) $(GVECTOR_H) $(GPOINT_H) $(GHVECTOR_H) Geometry/Rotation.h
Geometry/Rotation.o: $(MCONSTANTS_H) $(GMATRIX_H) $(GROTATION_H) Geometry/Rotation.cpp

GTRANSLATIONTRANSFORMATION_H = $(GVECTOR_H) $(GPOINT_H) $(GHVECTOR_H) Geometry/TranslationTransformation.h
Geometry/TranslationTransformation.o: $(GTRANSLATIONTRANSFORMATION_H) Geometry/TranslationTransformation.cpp

GROTATIONTRANSFORMATION_H = $(GVECTOR_H) $(GPOINT_H) $(GHVECTOR_H) $(GROTATION_H) Geometry/RotationTransformation.h
Geometry/RotationTransformation.o: $(GROTATIONTRANSFORMATION_H) Geometry/RotationTransformation.cpp

GORTHONORMALTRANSFORMATION_H = $(GVECTOR_H) $(GPOINT_H) $(GHVECTOR_H) $(GROTATION_H) $(GMATRIX_H) Geometry/OrthonormalTransformation.h
Geometry/OrthonormalTransformation.o: $(GTRANSLATIONTRANSFORMATION_H) $(GROTATIONTRANSFORMATION_H) $(GORTHONORMALTRANSFORMATION_H) Geometry/OrthonormalTransformation.cpp

GUNIFORMSCALINGTRANSFORMATION_H = $(GVECTOR_H) $(GPOINT_H) $(GHVECTOR_H) Geometry/UniformScalingTransformation.h
Geometry/UniformScalingTransformation.o: $(GUNIFORMSCALINGTRANSFORMATION_H) Geometry/UniformScalingTransformation.cpp

GORTHOGONALTRANSFORMATION_H = $(GVECTOR_H) $(GPOINT_H) $(GHVECTOR_H) $(GROTATION_H) Geometry/OrthogonalTransformation.h
Geometry/OrthogonalTransformation.o: $(GTRANSLATIONTRANSFORMATION_H) $(GROTATIONTRANSFORMATION_H) $(GORTHONORMALTRANSFORMATION_H) $(GUNIFORMSCALINGTRANSFORMATION_H) $(GORTHOGONALTRANSFORMATION_H) Geometry/OrthogonalTransformation.cpp

GSCALINGTRANSFORMATION_H = $(GCOMPONENTARRAY_H) $(GVECTOR_H) $(GPOINT_H) $(GHVECTOR_H) Geometry/ScalingTransformation.h
Geometry/ScalingTransformation.o: $(GUNIFORMSCALINGTRANSFORMATION_H) $(GSCALINGTRANSFORMATION_H) Geometry/ScalingTransformation.cpp

GAFFINETRANSFORMATION_H = $(GVECTOR_H) $(GPOINT_H) $(GHVECTOR_H) $(GMATRIX_H) $(GROTATION_H) Geometry/AffineTransformation.h
Geometry/AffineTransformation.o: $(GTRANSLATIONTRANSFORMATION_H) $(GROTATIONTRANSFORMATION_H) $(GORTHONORMALTRANSFORMATION_H) $(GUNIFORMSCALINGTRANSFORMATION_H) $(GORTHOGONALTRANSFORMATION_H) $(SCALINGTRANSFORMATION_H) $(GAFFINETRANSFORMATION_H) Geometry/AffineTransformation.cpp

GPROJECTIVETRANSFORMATION_H = $(GVECTOR_H) $(GPOINT_H) $(GHVECTOR_H) $(GMATRIX_H) $(GROTATION_H) Geometry/ProjectiveTransformation.h
Geometry/ProjectiveTransformation.o: $(GTRANSLATIONTRANSFORMATION_H) $(GROTATIONTRANSFORMATION_H) $(GORTHONORMALTRANSFORMATION_H) $(GUNIFORMSCALINGTRANSFORMATION_H) $(GORTHOGONALTRANSFORMATION_H) $(SCALINGTRANSFORMATION_H) $(GAFFINETRANSFORMATION_H) $(GPROJECTIVETRANSFORMATION_H) Geometry/ProjectiveTransformation.cpp

GRAY_H = $(GVECTOR_H) $(GPOINT_H) Geometry/Ray.h

GHITRESULT_H = $(MCONSTANTS_H) Geometry/HitResult.h

GSOLIDHITRESULT_H = $(GHITRESULT_H) Geometry/SolidHitResult.h

GPLANE_H = $(MMATH_H) $(GCOMPONENTARRAY_H) $(GVECTOR_H) $(GPOINT_H) $(GMATRIX_H) $(GRAY_H) $(GSOLIDHITRESULT_H) Geometry/Plane.h

GSPHERE_H = $(MMATH_H) $(GVECTOR_H) $(GPOINT_H) $(GORTHONORMALTRANSFORMATION_H) $(GORTHOGONALTRANSFORMATION_H) $(GRAY_H) $(GSOLIDHITRESULT_H) Geometry/Sphere.h

GCYLINDER_H = $(MMATH_H) $(MCONSTANTS_H) $(GVECTOR_H) $(GPOINT_H) $(GORTHONORMALTRANSFORMATION_H) $(GORTHOGONALTRANSFORMATION_H) $(GRAY_H) $(GSOLIDHITRESULT_H) Geometry/Cylinder.h

GCONE_H = $(MMATH_H) $(MCONSTANTS_H) $(GVECTOR_H) $(GPOINT_H) $(GORTHONORMALTRANSFORMATION_H) $(GORTHOGONALTRANSFORMATION_H) $(GRAY_H) $(GSOLIDHITRESULT_H) Geometry/Cone.h

GBOX_H = $(GCOMPONENTARRAY_H) $(GPOINT_H) $(GVECTOR_H) $(GRAY_H) $(GSOLIDHITRESULT_H) Geometry/Box.h
Geometry/Box.o: $(MCONSTANTS_H) $(GBOX_H) Geometry/Box.cpp

GPOLYGON_H = $(GPOINT_H) $(GVECTOR_H) $(GPLANE_H) Geometry/Polygon.h
Geometry/Polygon.o: $(UTILITY_H) $(GCOMPONENTARRAY_H) $(GAFFINECOMBINER_H) $(GPOLYGON_H) Geometry/Polygon.cpp

GSPLINECURVE_H = $(GVECTOR_H) $(GPOINT_H) Geometry/SplineCurve.h
Geometry/SplineCurve.o: $(GSPLINECURVE_H) Geometry/SplineCurve.cpp

GSPLINEPATCH_H = $(GCOMPONENTARRAY_H) $(GVECTOR_H) $(GPOINT_H) Geometry/SplinePatch.h
Geometry/SplinePatch.o: $(GSPLINEPATCH_H) Geometry/SplinePatch.cpp

GRANDOM_H = Geometry/Random.h
Geometry/Random.o: $(MMATH_H) $(MRANDOM_H) $(GPOINT_H) $(GVECTOR_H) $(GBOX_H) $(GRANDOM_H) Geometry/Random.cpp

GEOMETRY_O = Geometry/ComponentArray.o \
             Geometry/Vector.o \
             Geometry/Point.o \
             Geometry/HVector.o \
             Geometry/Matrix.o \
             Geometry/Rotation.o \
             Geometry/TranslationTransformation.o \
             Geometry/RotationTransformation.o \
             Geometry/OrthonormalTransformation.o \
             Geometry/UniformScalingTransformation.o \
             Geometry/OrthogonalTransformation.o \
             Geometry/ScalingTransformation.o \
             Geometry/AffineTransformation.o \
             Geometry/ProjectiveTransformation.o \
             Geometry/Box.o \
	     Geometry/Polygon.o \
             Geometry/SplineCurve.o \
             Geometry/SplinePatch.o \
             Geometry/Random.o

######################################################################
# OpenGL support library classes:
######################################################################

GLTEMPLATES_H = GLTemplates.h

GLTEXTURES_H = GLTextures.h

GLTYPES_H = $(GLTEMPLATES_H) GLTypes.h

GLVERTEX_H = $(GLTYPES_H) GLVertex.h

GLCOLOROPS_H = $(GLTYPES_H) GLColorOps.h

GLCOLORMAP_H = $(GLTYPES_H) GLColorMap.h
GLColorMap.o: $(FILEIO_H) $(ENDIANNESS_H) $(GLCOLORMAP_H) GLColorMap.cpp

GLMATERIAL_H = $(GLTYPES_H) GLMaterial.h
GLMaterial.o: $(GLMATERIAL_H) GLMaterial.cpp

GLTEXT_H =  GLText.h
GLText.o: CCFLAGS += -I$(FLTK_BASE)/$(INCDIR) 
GLText.o: $(GLTEXT_H) GLText.cpp

GLLIGHTSOURCE_H = $(GLTYPES_H) GLLightSource.h
GLLightSource.o: $(GLLIGHTSOURCE_H) GLLightSource.cpp

GLCONTEXTDATA_H = $(HASHTABLE_H) GLContextData.h

GLLINEILLUMINATOR_H = $(GLCONTEXTDATA_H) $(GLTYPES_H) $(GLMATERIAL_H) GLLineIlluminator.h
GLLineIlluminator.o: $(MMATH_H) $(GCOMPONENTARRAY_H) $(GMATRIX_H) $(GLCOLOROPS_H) $(GLLINEILLUMINATOR_H) GLLineIlluminator.cpp

GLGEOMETRY_H = $(GLTYPES_H) $(GCOMPONENTARRAY_H) $(GVECTOR_H) $(GPOINT_H) $(GHVECTOR_H) GLGeometry.h

GLTRANSFORMATIONS_H = $(GROTATION_H) $(GORTHONORMALTRANSFORMATION_H) $(GORTHOGONALTRANSFORMATION_H) $(GAFFINETRANSFORMATION_H) $(GPROJECTIVETRANSFORMATION_H) GLTransformations.h $(GLGEOMETRY_H) GLTransformations.cpp
GLTransformations.o: $(GLGEOMETRY_H) $(GLTRANSFORMATIONS_H) GLTransformations.cpp

GLSUPPORT_O = GLColorMap.o \
              GLMaterial.o \
              GLLightSource.o \
              GLText.o \
              GLLineIlluminator.o \
              GLTransformations.o

######################################################################
# Inverse kinematics classes:
######################################################################

IK_H = $(GVECTOR_H) $(GPOINT_H) $(GORTHONORMALTRANSFORMATION_H) $(GAFFINETRANSFORMATION_H) $(DENSEMATRIX_H) IK.h
IK_O = $(MMATH_H) $(MCONSTANTS_H) IK.o

######################################################################
# Protein visualization classes:
######################################################################

DRAGBOX_H = $(GVECTOR_H) $(GPOINT_H) $(GORTHONORMALTRANSFORMATION_H) $(GPROJECTIVETRANSFORMATION_H) DragBox.h
DragBox.o: $(MMATH_H) $(MCONSTANTS_H) $(GLTYPES_H) $(GLGEOMETRY_H) $(DRAGBOX_H) DragBox.cpp

MDGEOMETRY_H = $(GVECTOR_H) $(GPOINT_H) $(GRAY_H) MDGeometry.h

DISTANCERANGE_H = $(MDGEOMETRY_H) DistanceRange.h

ATOM_H = $(GSPHERE_H) $(MDGEOMETRY_H) Atom.h
Atom.o: $(ATOM_H) Atom.cpp

PROTEIN_H = $(MCONSTANTS_H) $(DRAGBOX_H) $(DISTANCERANGE_H) $(MDGEOMETRY_H) $(ATOM_H) Protein.h
Protein.o: $(MMATH_H) $(MCONSTANTS_H) $(GAFFINECOMBINER_H) $(GMATRIX_H) $(GAFFINETRANSFORMATION_H) $(GCYLINDER_H) $(PROTEIN_H) Protein.cpp

PROTEINFACTORY_H = $(PROTEIN_H) ProteinFactory.h
ProteinFactory.o: $(PROTEIN_H) $(PROTEINFACTROY_H) ProteinFactory.cpp

PARSEPDBFILE_H = ParsePdbFile.h
ParsePdbFile.o: $(PROTEIN_H) $(PARSEPDBFILE_H) ParsePdbFile.cpp

PROTEINRENDERER_H = $(CONFIGURATIONFILE_H) $(GBOX_H) $(GSPLINECURVE_H) $(GLVERTEX_H) $(GLCOLORMAP_H) $(GLMATERIAL_H) $(GLCONTEXTDATA_H) $(GLLINEILLUMINATOR_H) $(GLTEXT_H) $(PROTEIN_H) ProteinRenderer.h
GLSPHERE_CPP = $(MMATH_H) $(GLTEMPLATES_H) $(GLTYPES_H) GLSphere.cpp
ProteinRenderer.o: $(MCONSTANTS_H) $(GLGEOMETRY_H) $(PROTEINRENDERER_H) ProteinGui.h $(GLSPHERE_CPP) ProteinRenderer.cpp

CREATEPROTEIN_H = CreateProtein.h
ReadStandards.o: $(PROTEIN_H) $(CREATEPROTEIN_H) ReadStandards.cpp
ReadPredictionFile.o: $(PROTEIN_H) $(CREATEPROTEIN_H) ReadPredictionFile.cpp
SetDihedrals.o: $(PROTEIN_H) $(CREATEPROTEIN_H) SetDihedrals.cpp
BuildBeta.o: CCFLAGS += -I$(FLTK_BASE)/$(INCDIR)
BuildBeta.o: $(PROTEIN_H) BuildBeta.cpp
CreateProtein.o: $(PROTEIN_H) $(PARSEPDBFILE_H) $(CREATEPROTEIN_H) CreateProtein.cpp

PROTEIN_O = DragBox.o \
            Atom.o \
            Protein.o \
            ProteinFactory.o \
            ParsePdbFile.o \
            ProteinRenderer.o \
            ReadStandards.o \
            ReadPredictionFile.o \
            SetDihedrals.o \
            BuildBeta.o \
            CreateProtein.o

######################################################################
# Energy computation and visualization:
######################################################################

ENERGYAPI_H = EnergyAPI.h
EnergyAPI.o: $(PROTEIN_H) $(PARSEPDBFILE_H) $(ENERGYAPI_H) EnergyAPI.cpp

VOLUMERENDERER_H = $(GCOMPONENTARRAY_H) $(GHVECTOR_H) $(GMATRIX) $(GPOINT_H) $(GPROJECTIVETRANSFORMATION) $(GVECTOR_H) $(GLCONTEXTDATA_H) VolumeRenderer.h
VolumeRenderer.o: $(GHVECTOR_H) $(GPROJECTIVETRANSFORMATION_H) $(GLTEXTURES_H) $(PRIORITYHEAP_H) $(VOLUMERENDERER_H) VolumeRenderer.cpp

COLORFUNCTION_H = $(GLTYPES_H) ColorFunction.h
ColorFunction.o:  $(COLORFUNCTION_H) ColorFunction.cpp

ENERGYRENDERER_H = $(GLOBALS_H) $(COLORFUNCTION_H) EnergyRenderer.h
EnergyRenderer.o: $(ENERGYRENDERER_H) EnergyRenderer.cpp

ATOMCLASSIFIER_H = $(PROTEIN_H) AtomClassifier.h
AtomClassifier.o: $(ATOMCLASSIFIER_H) AtomClassifier.cpp

PROTEINRECORD_H = Buffer.h
ProteinRecord.o: $(PROTEINRECORD_H) ProteinRecord.cpp

ENERGYCOMPUTATION_O = EnergyAPI.o       \
                      EnergyRenderer.o  \
                      ColorFunction.o   \
                      AtomClassifier.o  \
                      VolumeRenderer.o  \
                      ProteinRecord.o

######################################################################
# Stride:
######################################################################

STRIDE_H = Stride.h
Stride.o: $(STRIDE_H) Stride.cpp

STRIDE_O = Stride.o

######################################################################
# Client/server communication for computational steering:
######################################################################

CLIENTSERVERPIPE_H = $(TCPSOCKET2_H) $(PROTEIN_H) ClientServerPipe.h

PROTEINCLIENT_H = ProteinClient.h
ProteinClient.o: $(CLIENTSERVERPIPE_H) $(PROTEINCLIENT_H) ProteinClient.cpp
	$(CCOMP) -c -o ProteinClient.o $(CCFLAGS) -DCLIENT_SIDE ProteinClient.cpp

CLIENTSERVER_O = ProteinClient.o

######################################################################
# Driver programs:
######################################################################

UNDOBUFFER_H = $(GORTHONORMALTRANSFORMATION_H) $(PROTEIN_H) UndoBuffer.h
UndoBuffer.o: $(GAFFINETRANSFORMATION_H) $(UNDOBUFFER_H) UndoBuffer.cpp

GLOBALS_H = $(CONFIGURATIONFILE_H) $(GBOX_H) $(GLTYPES_H) $(GLTEXT_H) $(DRAGBOX_H) $(PROTEIN_H) $(PROTEINCLIENT_H) $(PROTEINMANIPULATOR_H) $(PROTEINRENDERER_H) $(ENERGYAPI_H) $(UNDOBUFFER_H) Globals.h
Globals.o: CCFLAGS += -I$(FLTK_BASE)/$(INCDIR)
Globals.o: $(MMATH_H) $(GAFFINETRANSFORMATION_H) $(GLTRANSFORMATIONS_H) $(GLCONTEXTDATA_H) $(GLOBALS_H) Globals.cpp

SPACEBALL_H = $(SERIALPORT_H) SpaceBall.h
SpaceBall.o: $(SPACEBALL_H) SpaceBall.cpp

MYSPACEBALL_H = $(CONFIGURATIONFILE_H) $(SPACEBALL_H) MySpaceBall.h
MySpaceBall.o: $(MMATH_H) $(MYSPACEBALL_H) MySpaceBall.cpp

MYFLGLWINDOW_H = $(GVECTOR_H) $(GPOINT_H) $(GRAY_H) $(GORTHONORMALTRANSFORMATION_H) $(GPROJECTIVETRANSFORMATION_H) $(GLTYPES_H) $(GLLIGHTSOURCE_H) $(GLCONTEXTDATA_H) MyFlGlWindow.h
#MyFlGlWindow.o: CCFLAGS += -I$(FLTK_BASE)/$(INCDIR)
MyFlGlWindow.o: $(MMATH_H) $(GLTEMPLATES_H) $(GLGEOMETRY_H) $(GLTRANSFORMATIONS_H) $(MYSPACEBALL_H) $(MYFLGLWINDOW_H) $(GLOBALS_H) MyFlGlWindow.cpp

MYROLLER_H = MyRoller.h
MyRoller.o: CCFLAGS += -I$(FLTK_BASE)/$(INCDIR)
MyRoller.o: $(MYROLLER_H) $(GLOBALS_H) MyRoller.cpp

SEQUENCEWIDGET_H = SequenceWidget.h
SequenceWidget.o: CCFLAGS += -I$(FLTK_BASE)/$(INCDIR)
SequenceWidget.o: $(PROTEIN_H) $(PROTEINRENDERER_H) $(SEQUENCEWIDGET_H) SequenceWidget.cpp

RAMACHANDRANPLOT_H = $(GLTYPES_H) $(MDGEOMETRY_H) RamachandranPlot.h
RamachandranPlot.o: CCFLAGS += -I$(FLTK_BASE)/$(INCDIR)
RamachandranPlot.o: $(MCONSTANTS_H) $(GLTEMPLATES_H) $(RAMACHANDRANPLOT_H) $(GLOBALS_H) RamachandranPlot.cpp

ProteinGui.h ProteinGui.cpp: ProteinGui.fl
	$(FLTK_BASE)/$(BINDIR)/fluid -c -h ProteinGui.h -o ProteinGui.cpp ProteinGui.fl

ProteinGui.o: CCFLAGS += -I$(FLTK_BASE)/$(INCDIR)
ProteinGui.o: $(MYFLGLWINDOW_H) $(MYROLLER_H) $(SEQUENCEWIDGET_H) ProteinGui.h ProteinGui.cpp

ProteinFltk.o: CCFLAGS += -I$(FLTK_BASE)/$(INCDIR)
ProteinFltk.o: $(MMATH_H) $(MCONSTANTS_H) $(PARSEPDBFILE_H) $(CREATEPROTEIN_H) ProteinGui.h $(GLOBALS_H) ProteinFltk.cpp

DRIVER_O = UndoBuffer.o \
           ProteinInteractor.o \
           Globals.o \
           SpaceBall.o \
           MySpaceBall.o \
           MyFlGlWindow.o \
           MyRoller.o \
           SequenceWidget.o \
           RamachandranPlot.o \
           ProteinGui.o \
           ProteinFltk.o

######################################################################
# Fortran language components:
######################################################################


FC = g77
FFLAGS = -w
FORTRAN_O = fitsq.o


######################################################################
# ProteinShop binary:
######################################################################

PROTEINSHOP_O = $(MISC_O)               \
                $(MATH_O)               \
                $(GEOMETRY_O)           \
                $(GLSUPPORT_O)          \
                $(IK_O)                 \
                $(PROTEIN_O)            \
                $(ENERGYCOMPUTATION_O)  \
                $(CLIENTSERVER_O)       \
                $(STRIDE_O)             \
                $(DRIVER_O)             \
                $(FORTRAN_O)

PROTEINSHOP_DEPS = $(PROTEINSHOP_O)

PROTEINSHOP_LIBS = $(FLTK_LIBS) \
                   $(GL_LIBS)   \
                   -lpthread    \
                   -lm -ldl

# For SGI's CC:
ifeq ($(CCOMP),CC)
  PROTEINSHOP_LIBS += -lCio
endif

ifeq ($(SYSTEM),Darwin)
ProteinShop: $(PROTEINSHOP_DEPS)
	$(CCOMP) $(CMAKEEXEFLAGS) -o ProteinShop $(PROTEINSHOP_O) $(PROTEINSHOP_LIBS)
	$(POSTBUILD) $@ $(FLTK_BASE)/$(INCDIR)/FL/mac.r
	$(INSTALL)
else
ProteinShop: $(PROTEINSHOP_DEPS)
	$(CCOMP) $(CMAKEEXEFLAGS) -o ProteinShop $(PROTEINSHOP_O) $(PROTEINSHOP_LIBS)
	$(INSTALL)
endif


# DO NOT DELETE
