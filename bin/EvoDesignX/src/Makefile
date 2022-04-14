CXX  = g++
CXXFLAGS = -O3
CXXFLAGS += -w
LDFLAGS  = -static
LDFLAGS += -lgsl
OBJ  = Atom.o AtomParamsSet.o Chain.o DEE.o DesignSite.o EEF1Ligand.o EnergyFunction.o EnergyMatrix.o EnergyOptimization.o ErrorTracker.o EvoAminoName.o EvoGetPhiPsi.o EvoGetSeq2SA.o EvoGetSeq2SS.o Evolution.o EvoSeqAlign.o EvoUtility.o FlexibleBackbone.o GeometryCalc.o Getopt.o kabsch.o ProgramFunction.o Residue.o ResidueTopology.o RotamerBuilder.o Rotamer.o RotamerOptimizer.o Sequence.o SmallMol.o Structure.o Utility.o WeightOpt.o 

PROG = designx5

all: ${PROG}

designx5: ${OBJ} Main.cpp ProgramFunction.h Getopt.h RotamerBuilder.h Structure.h EnergyMatrix.h Evolution.h DEE.h EnergyOptimization.h WeightOpt.h EEF1Ligand.h EvoAminoName.h FlexibleBackbone.h
	${CXX} ${OBJ} Main.cpp ${CXXFLAGS} -o $@ ${LDFLAGS}


Atom.o: Atom.cpp Atom.h ErrorTracker.h Utility.h GeometryCalc.h
	${CXX} -c Atom.cpp ${CXXFLAGS} -o $@

AtomParamsSet.o: AtomParamsSet.cpp AtomParamsSet.h Atom.h
	${CXX} -c AtomParamsSet.cpp ${CXXFLAGS} -o $@

Chain.o: Chain.cpp Chain.h Residue.h
	${CXX} -c Chain.cpp ${CXXFLAGS} -o $@

DEE.o: DEE.cpp DEE.h SmallMol.h EnergyMatrix.h
	${CXX} -c DEE.cpp ${CXXFLAGS} -o $@

DesignSite.o: DesignSite.cpp DesignSite.h Rotamer.h
	${CXX} -c DesignSite.cpp ${CXXFLAGS} -o $@

EEF1Ligand.o: EEF1Ligand.cpp EEF1Ligand.h Atom.h
	${CXX} -c EEF1Ligand.cpp ${CXXFLAGS} -o $@

EnergyFunction.o: EnergyFunction.cpp EnergyFunction.h Residue.h Rotamer.h
	${CXX} -c EnergyFunction.cpp ${CXXFLAGS} -o $@

EnergyMatrix.o: EnergyMatrix.cpp EnergyMatrix.h Utility.h Structure.h ErrorTracker.h EnergyFunction.h
	${CXX} -c EnergyMatrix.cpp ${CXXFLAGS} -o $@

EnergyOptimization.o: EnergyOptimization.cpp EnergyOptimization.h Sequence.h Evolution.h FlexibleBackbone.h DEE.h EnergyMatrix.h SmallMol.h Sequence.h 
	${CXX} -c EnergyOptimization.cpp ${CXXFLAGS} -o $@

ErrorTracker.o: ErrorTracker.cpp ErrorTracker.h 
	${CXX} -c ErrorTracker.cpp ${CXXFLAGS} -o $@

EvoAminoName.o: EvoAminoName.cpp EvoAminoName.h
	${CXX} -c EvoAminoName.cpp ${CXXFLAGS} -o $@




EvoGetPhiPsi.o: EvoGetPhiPsi.cpp EvoGetPhiPsi.h EvoUtility.h
	${CXX} -c EvoGetPhiPsi.cpp ${CXXFLAGS} -o $@

EvoGetSeq2SA.o: EvoGetSeq2SA.cpp EvoGetSeq2SA.h EvoUtility.h
	${CXX} -c EvoGetSeq2SA.cpp ${CXXFLAGS} -o $@

EvoGetSeq2SS.o: EvoGetSeq2SS.cpp EvoGetSeq2SS.h Utility.h
	${CXX} -c EvoGetSeq2SS.cpp ${CXXFLAGS} -o $@

Evolution.o: Evolution.cpp Evolution.h EvoSeqAlign.h EvoGetSeq2SA.h EvoGetSeq2SS.h EvoGetPhiPsi.h EvoUtility.h ErrorTracker.h EvoAminoName.h
	${CXX} -c Evolution.cpp ${CXXFLAGS} -o $@

EvoSeqAlign.o: EvoSeqAlign.cpp EvoSeqAlign.h EvoAminoName.h EvoUtility.h
	${CXX} -c EvoSeqAlign.cpp ${CXXFLAGS} -o $@

EvoUtility.o: EvoUtility.cpp EvoUtility.h 
	${CXX} -c EvoUtility.cpp ${CXXFLAGS} -o $@

FlexibleBackbone.o: FlexibleBackbone.cpp FlexibleBackbone.h RotamerBuilder.h RotamerOptimizer.h Residue.h Utility.h kabsch.h EnergyFunction.h Structure.h
	${CXX} -c FlexibleBackbone.cpp ${CXXFLAGS} -o $@

GeometryCalc.o: GeometryCalc.cpp GeometryCalc.h ErrorTracker.h 
	${CXX} -c GeometryCalc.cpp ${CXXFLAGS} -o $@

Getopt.o: Getopt.cpp Getopt.h
	${CXX} -c Getopt.cpp ${CXXFLAGS} -o $@

kabsch.o: kabsch.cpp kabsch.h 
	${CXX} -c kabsch.cpp ${CXXFLAGS} -o $@ ${LDFLAGS}

ProgramFunction.o: ProgramFunction.cpp ProgramFunction.h RotamerBuilder.h RotamerOptimizer.h Evolution.h EnergyFunction.h Structure.h
	${CXX} -c ProgramFunction.cpp ${CXXFLAGS} -o $@

Residue.o: Residue.cpp Residue.h AtomParamsSet.h ResidueTopology.h
	${CXX} -c Residue.cpp ${CXXFLAGS} -o $@

ResidueTopology.o: ResidueTopology.cpp ResidueTopology.h Utility.h Atom.h
	${CXX} -c ResidueTopology.cpp ${CXXFLAGS} -o $@

RotamerBuilder.o: RotamerBuilder.cpp RotamerBuilder.h EnergyFunction.h Structure.h
	${CXX} -c RotamerBuilder.cpp ${CXXFLAGS} -o $@

Rotamer.o: Rotamer.cpp Rotamer.h Atom.h Residue.h
	${CXX} -c Rotamer.cpp ${CXXFLAGS} -o $@

RotamerOptimizer.o: RotamerOptimizer.cpp RotamerOptimizer.h EnergyFunction.h Structure.h
	${CXX} -c RotamerOptimizer.cpp ${CXXFLAGS} -o $@

Sequence.o: Sequence.cpp Sequence.h Structure.h EnergyMatrix.h EnergyFunction.h 
	${CXX} -c Sequence.cpp ${CXXFLAGS} -o $@

SmallMol.o: SmallMol.cpp SmallMol.h DesignSite.h Rotamer.h
	${CXX} -c SmallMol.cpp ${CXXFLAGS} -o $@

Structure.o: Structure.cpp Structure.h EnergyFunction.h Chain.h Rotamer.h DesignSite.h SmallMol.h
	${CXX} -c Structure.cpp ${CXXFLAGS} -o $@

Utility.o: Utility.cpp Utility.h ErrorTracker.h
	${CXX} -c Utility.cpp ${CXXFLAGS} -o $@

WeightOpt.o: WeightOpt.cpp WeightOpt.h ErrorTracker.h Utility.h EnergyFunction.h 

install: ${PROG}
	cp ${PROG} ../

clean:
	rm *.o ${PROG}




