/*******************************************************************************************************************************
Copyright (c) 2020 Xiaoqiang Huang (tommyhuangthu@foxmail.com, xiaoqiah@umich.edu)

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation 
files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, 
modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the 
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE 
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR 
IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
********************************************************************************************************************************/

#ifndef STRUCTURE_H
#define STRUCTURE_H

#include "Chain.h"
#include "Rotamer.h"
#include "DesignSite.h"
#include "SmallMol.h"


typedef struct _Structure{
  char name[MAX_LENGTH_STRUCTURE_NAME+1];
  Chain* chains;
  DesignSite *designSites;
  int chainNum;
  int designSiteCount;
} Structure;

typedef struct _StructureFlex{
  ChainFlex *chainsFlex;
  int chainNum;
} StructureFlex;

int StructureCreate(Structure* pThis);
int StructureFlexCreate(StructureFlex* pThis);
int StructureDestroy(Structure* pThis);
int StructureFlexDestroy(StructureFlex* pThis);
char* StructureGetName(Structure* pThis);
int StructureSetName(Structure* pThis, char* newName);
int StructureGetChainCount(Structure* pThis);
int StructureFlexGetChainCount(StructureFlex* pThis);
Chain* StructureGetChain(Structure* pThis, int index);
ChainFlex* StructureFlexGetChain(StructureFlex* pThis, int index);
Chain* StructureFindChainByName(Structure* pThis, char* chainName);
int StructureFindChainIndex(Structure* pThis, char* chainName, int* index);
int StructureFlexFindChainIndex(StructureFlex* pThis, char* chainName, int* index);
int StructureFindSmallMol(Structure* pThis, Residue** ppSmallMol);
int StructureAddChain(Structure* pThis, Chain* newChain);
int StructureDeleteChain(Structure* pThis, char* chainName);
int StructureShowInPDBFormat(Structure* pThis, BOOL showHydrogen, FILE* pFile);

//deal with design sites
int StructureGetDesignSiteCount(Structure* pThis);
DesignSite* StructureGetDesignSite(Structure* pThis, int index);
DesignSite* StructureFindDesignSite(Structure* pThis, int chainIndex, int resiIndex);
int StructureShowDesignSites(Structure* pThis, FILE* pFile);
int ProteinSiteAddDesignSite(Structure* pThis, int chainIndex, int resiIndex);
int StructureRemoveAllDesignSites(Structure* pThis);
int ProteinSiteRemoveDesignSite(Structure* pThis, int chainIndex, int resiIndex);

//other functions
int ChainComputeResiduePosition(Structure *pStructure, int chainIndex);
int StructureComputeResiduePosition(Structure *pStructure);
int StructureGetAminoAcidComposition(Structure* pStructure, int *aas);

//checker and debuggers
int StructureShowAtomParameter(Structure* pStructure);
int StructureShowBondInformation(Structure* pStructure);

int StructureConfig(Structure *pStructure, char* pdbFile, AtomParamsSet* pAtomParams, ResiTopoSet* pResiTopos);
int StructureFlexConfig(Structure *pStructure,StructureFlex *pStructureFlex);
int StructureConfigLigand(Structure* pStructure, char* mol2file, AtomParamsSet* pAtomParams, ResiTopoSet* pTopos);

DesignSite * StructureFindDesignSiteByChainName(Structure *pStructure, char *chainName, int posInChain);
int DesignSiteIndexFind(Structure *pStructure, char *chainName, int posInSeg);
int DesignSiteIndexGet(Structure *pStructure,int chainIndex, int resiIndex);
int StructureCalcProteinResidueSidechainTorsion(Structure* pThis, ResiTopoSet* pResiTopos);

int StructureCalcPhiPsi(Structure* pStructure);
int StructureFlexCalcAngles(Structure *pStructure,StructureFlex *pStructureFlex);
int StructureFlexCalcOmega(Structure *pStructure,StructureFlex *pStructureFlex);
int StructureFlexCalcBondLengths(Structure *pStructure,StructureFlex *pStructureFlex);
int StructureCopy(Structure* pThis, Structure* pOther);

#endif // STRUCTURE_H
