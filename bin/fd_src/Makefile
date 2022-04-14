CXX  = g++
CXXFLAGS = -O3
CXXFLAGS += -w
#LDFLAGS  = -static
# flags for gprof, it is recommended by 
# https://users.cs.duke.edu/~ola/courses/programming/gprof.html to
# leave out the -g option and include the -O (or -O2 or -O3) option
# to avoid biased timing information.
# To profile the program, first run the full program to generate gmon.out
# and then uses gprof to get statistics
#CXXFLAGS += -pg
#CXXFLAGS += -ggdb
#CXXFLAGS += -Wall

OBJ  = CommonPara.o ParseInput.o InputData.o BasicFunc.o GeometryCalc.o Getopt.o Fragments.o SegCluster.o ProgramFunctions.o ParsePDB.o PrintFunc.o ParseGeometryFiles.o ParseFoldDesignMovementFiles.o FoldDesignMovement.o ParseFoldDesignEnergyFiles.o FoldDesignEnergyFunction.o 
PROG = FoldDesign

#### compile main programs ####

all: ${PROG} 

FoldDesign: ${OBJ} main.cpp ProgramFunctions.h Getopt.h 
	${CXX} ${OBJ} main.cpp ${CXXFLAGS} -o $@ ${LDFLAGS}

#### compile object files ####

CommonPara.o: CommonPara.cpp CommonPara.h
	${CXX} -c CommonPara.cpp ${CXXFLAGS} -o $@

BasicFunc.o: BasicFunc.cpp CommonPara.h BasicFunc.h
	${CXX} -c BasicFunc.cpp ${CXXFLAGS} -o $@


ParseInput.o: ParseInput.cpp CommonPara.h ParseInput.h
	${CXX} -c ParseInput.cpp ${CXXFLAGS} -o $@


ParsePDB.o: ParsePDB.cpp CommonPara.h ParsePDB.h BasicFunc.h
	${CXX} -c ParsePDB.cpp ${CXXFLAGS} -o $@


Getopt.o: Getopt.cpp Getopt.h
	${CXX} -c Getopt.cpp ${CXXFLAGS} -o $@

Fragments.o: Fragments.cpp CommonPara.h InputData.h ParsePDB.h BasicFunc.h Fragments.h SegCluster.h
	${CXX} -c Fragments.cpp ${CXXFLAGS} -o $@


InputData.o: InputData.cpp CommonPara.h ParseInput.h InputData.h
	${CXX} -c InputData.cpp ${CXXFLAGS} -o $@


ParseFoldDesignEnergyFiles.o: ParseFoldDesignEnergyFiles.cpp CommonPara.h BasicFunc.h GeometryCalc.h InputData.h ParseFoldDesignEnergyFiles.h
	${CXX} -c ParseFoldDesignEnergyFiles.cpp ${CXXFLAGS} -o $@


ProgramFunctions.o: ProgramFunctions.cpp CommonPara.h ParseInput.h ParsePDB.h Fragments.h PrintFunc.h ProgramFunctions.h InputData.h ParseGeometryFiles.h ParseFoldDesignEnergyFiles.h FoldDesignEnergyFunction.h 
	${CXX} -c ProgramFunctions.cpp ${CXXFLAGS} -o $@


SegCluster.o: SegCluster.cpp CommonPara.h BasicFunc.h ParsePDB.h SegCluster.h
	${CXX} -c SegCluster.cpp ${CXXFLAGS} -o $@


ParseFoldDesignMovementFiles.o: ParseFoldDesignMovementFiles.cpp CommonPara.h GeometryCalc.h FoldDesignEnergyFunction.h BasicFunc.h ParsePDB.h InputData.h ParseFoldDesignMovementFiles.h PrintFunc.h
	${CXX} -c ParseFoldDesignMovementFiles.cpp ${CXXFLAGS} -o $@


FoldDesignMovement.o: FoldDesignMovement.cpp CommonPara.h PrintFunc.h FoldDesignMovement.h InputData.h ParsePDB.h BasicFunc.h ParseFoldDesignMovementFiles.h FoldDesignEnergyFunction.h GeometryCalc.h
	${CXX} -c FoldDesignMovement.cpp ${CXXFLAGS} -o $@


FoldDesignEnergyFunction.o: FoldDesignEnergyFunction.cpp CommonPara.h BasicFunc.h FoldDesignEnergyFunction.h ParseFoldDesignEnergyFiles.h ParseInput.h GeometryCalc.h
	${CXX} -c FoldDesignEnergyFunction.cpp ${CXXFLAGS} -o $@ 


ParseGeometryFiles.o: ParseGeometryFiles.cpp CommonPara.h GeometryCalc.h ParseGeometryFiles.h
	${CXX} -c ParseGeometryFiles.cpp ${CXXFLAGS} -o $@

GeometryCalc.o: GeometryCalc.cpp CommonPara.h BasicFunc.h GeometryCalc.h ParseGeometryFiles.h
	${CXX} -c GeometryCalc.cpp ${CXXFLAGS} -o $@


#### auxiliary operations ####

install: ${PROG}
	cp ${PROG} ../

clean:
	rm *.o ${PROG}
