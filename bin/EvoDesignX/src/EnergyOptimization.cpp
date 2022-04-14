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

#pragma warning(disable:4244)
#pragma warning(disable:4305)
#include "EnergyOptimization.h"
#include "Sequence.h"
#include "Evolution.h"
#include "FlexibleBackbone.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

extern BOOL FLAG_EVOLUTION;
extern BOOL FLAG_PHYSICS;
extern BOOL FLAG_MONOMER;
extern BOOL FLAG_PPI;
extern BOOL FLAG_PROT_LIG;
extern BOOL FLAG_ENZYME;
extern BOOL FLAG_MOL2;
extern BOOL FLAG_EVOPHIPSI;

extern char BEST_SEQ[MAX_LENGTH_ONE_LINE_IN_FILE+1];
extern char BEST_STRUCT[MAX_LENGTH_ONE_LINE_IN_FILE+1];
extern char BEST_MOL2[MAX_LENGTH_FILE_NAME+1];

extern int TOT_SEQ_LEN;
extern int NTRAJ;

extern char PDBID[MAX_LENGTH_FILE_NAME+1];
extern char DES_CHAINS[10];
extern double WGT_PROFILE;

int SitePairConsDeploy(CataConsSitePairArray *pSitePairArray, Structure *pStructure)
{
  int i;
  for(i=0; i<CataConsSitePairArrayGetCount(pSitePairArray); i++){
    CataConsSitePair *pSitePair   = CataConsSitePairArrayGet(pSitePairArray, i);
    DesignSite *pFirstDesignSite  = StructureFindDesignSiteByChainName(pStructure, pSitePair->firstSiteChainName, pSitePair->firstSitePosInChain);
    DesignSite *pSecondDesignSite = StructureFindDesignSiteByChainName(pStructure, pSitePair->secondSiteChainName, pSitePair->secondSitePosInChain);
    CataConsSitePairDeploy(pSitePair, DesignSiteGetRotamers(pFirstDesignSite), DesignSiteGetRotamers(pSecondDesignSite));
  }
  return Success;
}

int EnergyMatrixUpdateForCataConsNew(EnergyMatrix *pMatrix, RotamerList* pList, CataConsSitePair *pSitePair, Structure *pStructure)
{
  DesignSite *pFirstDesignSite  = StructureFindDesignSiteByChainName(pStructure, pSitePair->firstSiteChainName, pSitePair->firstSitePosInChain);
  DesignSite *pSecondDesignSite = StructureFindDesignSiteByChainName(pStructure, pSitePair->secondSiteChainName, pSitePair->secondSitePosInChain);
  int firstDesignSiteIndex      = DesignSiteIndexFind(pStructure, pSitePair->firstSiteChainName, pSitePair->firstSitePosInChain);
  int secondDesignSiteIndex     = DesignSiteIndexFind(pStructure, pSitePair->secondSiteChainName, pSitePair->secondSitePosInChain);
  EnergyMatrixBlock* pBlockII = EnergyMatrixGetBlock(pMatrix, firstDesignSiteIndex, firstDesignSiteIndex);
  EnergyMatrixBlock* pBlockKK = EnergyMatrixGetBlock(pMatrix, secondDesignSiteIndex, secondDesignSiteIndex);
  for(int j=0; j<pBlockII->RotamerCountSiteI; j++){
    int trueIndexIJ;
    Rotamer *pRotamerJ;
    RotamerOriginalIndexGet(pList, firstDesignSiteIndex, j, &trueIndexIJ);
    pRotamerJ = RotamerSetGet(DesignSiteGetRotamers(pFirstDesignSite), trueIndexIJ);
    for(int s=0; s<pBlockKK->RotamerCountSiteK; s++){
      int trueIndexKS;
      Rotamer *pRotamerS;
      RotamerOriginalIndexGet(pList, secondDesignSiteIndex, s, &trueIndexKS);
      pRotamerS = RotamerSetGet(DesignSiteGetRotamers(pSecondDesignSite), trueIndexKS);
      if(CataConsSitePairCheck(pSitePair, pRotamerJ, pRotamerS) == FALSE){
        if(firstDesignSiteIndex < secondDesignSiteIndex){
          *EnergyMatrixGet(pMatrix, firstDesignSiteIndex, secondDesignSiteIndex, j, s) += 1e4;
        }
        else{
          *EnergyMatrixGet(pMatrix, secondDesignSiteIndex, firstDesignSiteIndex, s, j) += 1e4;
        }
      }
    }
  }
  return Success;
}

int EnergyMatrixUpdateForCataConsArrayNew(EnergyMatrix *pMatrix, RotamerList* pList, CataConsSitePairArray *pSitePairArray, Structure *pStructure)
{
  int i;
  for(i=0; i<CataConsSitePairArrayGetCount(pSitePairArray); i++){
    CataConsSitePair *pSitePair = CataConsSitePairArrayGet(pSitePairArray, i);
    EnergyMatrixUpdateForCataConsNew(pMatrix, pList, pSitePair, pStructure);
  }
  return Success;
}




int DesignSiteGetRotamerTypeAndCount(RotamerList* pList, Structure* pStructure, StringArray** ppRotamerType, IntArray** ppRotamerCount){
  int totalRemainRotamerCount=0;
  int designSiteCount=StructureGetDesignSiteCount(pStructure);
  StringArray* pStringArray = *ppRotamerType;
  IntArray*	 pIntArray	  = *ppRotamerCount;

  printf("show rotamer distribution at each design position\n");
  for(int i=0; i<designSiteCount; i++){
    int remainRotamerCount=0;
    DesignSite* pDesignSite = StructureGetDesignSite(pStructure, i);
    RotamerSet* pRotamerSet = DesignSiteGetRotamers(pDesignSite);
    IntArrayCreate(&pIntArray[i], 0);
    StringArrayCreate(&pStringArray[i]);
    for(int j=0; j<pList->rotamerCount[i]; j++){
      if(pList->remainFlag[i][j]==FALSE) continue;
      BOOL     isRotamerTypeNew=TRUE;
      Rotamer* pRotamer=RotamerSetGet(pRotamerSet, j);
      for(int k=0; k<StringArrayGetCount(&pStringArray[i]); k++){
        if(strcmp(RotamerGetType(pRotamer),StringArrayGet(&pStringArray[i],k)) == 0){
          isRotamerTypeNew=FALSE;
          break;
        }
      }
      if(isRotamerTypeNew==TRUE){
        int length=IntArrayGetLength(&pIntArray[i]);
        StringArrayAppend(&pStringArray[i],RotamerGetType(pRotamer));
        length++;
        IntArrayResize(&pIntArray[i],length);
        IntArraySet(&pIntArray[i],length-1,1);
      }else{
        int pos, rotamerCount;
        StringArrayFind(&pStringArray[i],RotamerGetType(pRotamer),&pos);
        rotamerCount=IntArrayGet(&pIntArray[i],pos);
        rotamerCount++;
        IntArraySet(&pIntArray[i],pos,rotamerCount);
      }
    }
    for(int j=0; j<IntArrayGetLength(&pIntArray[i]); j++){
      remainRotamerCount += IntArrayGet(&pIntArray[i],j);
    }
    totalRemainRotamerCount += remainRotamerCount;
    printf("design site %3d : %3s %s %4d, %6d rotamers:  ",i,
      ResidueGetName(pDesignSite->pResidue),
      ResidueGetChainName(pDesignSite->pResidue),
      ResidueGetPosInChain(pDesignSite->pResidue),
      remainRotamerCount);
    for(int j=0; j<StringArrayGetCount(&pStringArray[i]); j++){
      printf("%4d %s ",IntArrayGet(&pIntArray[i],j),StringArrayGet(&pStringArray[i],j));
    }
    printf("\n");
  }

  return Success;
}



int SequenceRandomSingleSiteIndex(int *mutSiteIndex, int designSiteCount){
  *mutSiteIndex=rand()%designSiteCount;
  return Success;
}

int SequenceRandomSingleRotamerIndex(Structure* pStructure, RotamerList* pList, int mutSiteIndex, int *mutRotIndex){
  int reducedRotIndex =rand()%pList->remainRotamerCount[mutSiteIndex];
  int oriRotIndex=-1;
  RotamerOriginalIndexGet(pList,mutSiteIndex,reducedRotIndex,&oriRotIndex);
  *mutRotIndex=oriRotIndex;
  return Success;
}

int SequenceRandomSingleRotamerIndexNew(Structure* pStructure,RotamerList* pList,StringArray** ppRotamerType,IntArray** ppRotamerCount,int siteIndex,int* rotIndex){
  StringArray* pStringArray = *ppRotamerType;
  IntArray*    pIntArray    = *ppRotamerCount;
  int typeIndex=rand()%IntArrayGetLength(&pIntArray[siteIndex]);
  int indexInType=rand()%IntArrayGet(&pIntArray[siteIndex],typeIndex);
  int typeflag = 0;
  for(int i = 0; i < pList->rotamerCount[siteIndex]; i++){
    Rotamer* pRotamer=RotamerSetGet(DesignSiteGetRotamers(StructureGetDesignSite(pStructure,siteIndex)),i);
    if(pList->remainFlag[siteIndex][i] == FALSE) continue;
    if(strcmp(StringArrayGet(&pStringArray[siteIndex],typeIndex),RotamerGetType(pRotamer)) != 0) continue;
    if(typeflag == indexInType){
      *rotIndex=i;
      break;
    }
    typeflag++;
  }
  return Success;
}



int SequenceChangeSingleSite(Sequence* pThis, int mutationSiteIndex, int mutationRotamerIndex){
  IntArraySet(&pThis->rotamerIndices, mutationSiteIndex, mutationRotamerIndex);
  return Success;
}



//////////////////////////////////////////////////////
//calculate energy of a sequence
//////////////////////////////////////////////////////
int StructureCalcSequenceTemplateEnergy(Structure* pStructure,Sequence* pSequence,double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM],double energyTermsBind[MAX_EVOEF_ENERGY_TERM_NUM]){
  for(int i = 0; i < StructureGetChainCount(pStructure); i++){
    Chain *pChainI = StructureGetChain(pStructure,i);
    for(int ir = 0; ir < ChainGetResidueCount(pChainI); ir++){
      Residue *pResIR = ChainGetResidue(pChainI,ir);
      if(pResIR->designSiteType!=Type_ResidueDesignType_Fixed){
        //the residue is a design site, set the coordinates of non-backbone atoms to be FALSE
        for(int atom1=0; atom1<ResidueGetAtomCount(pResIR); atom1++){
          Atom* pAtom1 = ResidueGetAtom(pResIR, atom1);
          if(!pAtom1->isBBAtom) pAtom1->isXyzValid=FALSE;
        }
      }
    }
  }

  for(int i=0;i<StructureGetChainCount(pStructure);i++){
    Chain *pChainI=StructureGetChain(pStructure,i);
    for(int ir=0;ir<ChainGetResidueCount(pChainI);ir++){
      Residue *pResIR=ChainGetResidue(pChainI,ir);
      if(pResIR->designSiteType==Type_ResidueDesignType_Fixed){
        energyTerms[91]+=pResIR->aapp;
        energyTerms[92]+=pResIR->ramachandran;
        energyTerms[93]+=pResIR->dunbrack;
        EVOEF_AminoAcidReferenceEnergy(pResIR->name,energyTerms);
        EVOEF_EnergyResidueIntraEnergy(pResIR,energyTerms);
      }
      for(int is=ir+1;is<ChainGetResidueCount(pChainI);is++){
        Residue *pResIS=ChainGetResidue(pChainI,is);
        if(ResidueGetPosInChain(pResIR)+1==ResidueGetPosInChain(pResIS)) EVOEF_EnergyResidueAndNextResidue(pResIR,pResIS,energyTerms);
        else EVOEF_EnergyResidueAndOtherResidueSameChain(pResIR,pResIS,energyTerms);
        //EVOEF_EnergyResidueAndResidueSameChain(pResIR,pResIS,energyTerms);
      }
      for(int k=i+1; k<StructureGetChainCount(pStructure); k++){
        Chain *pChainK = StructureGetChain(pStructure,k);
        for(int ks=0; ks<ChainGetResidueCount(pChainK); ks++){
          Residue *pResKS = ChainGetResidue(pChainK,ks);
          if(ChainGetType(pChainI)==Type_Chain_SmallMol){
            EVOEF_EnergyResidueAndLigandResidue(pResKS,pResIR,energyTerms);
            if(FLAG_PROT_LIG || FLAG_ENZYME){
              EVOEF_EnergyResidueAndLigandResidue(pResKS,pResIR,energyTermsBind);
            }
          }
          else if(ChainGetType(pChainK)==Type_Chain_SmallMol){
            EVOEF_EnergyResidueAndLigandResidue(pResIR,pResKS,energyTerms);
            if(FLAG_PROT_LIG || FLAG_ENZYME){
              EVOEF_EnergyResidueAndLigandResidue(pResIR,pResKS,energyTermsBind);
            }
          }
          else{
            EVOEF_EnergyResidueAndOtherResidueDifferentChain(pResIR,pResKS,energyTerms);
            //EVOEF_EnergyResidueAndOtherResidueDifferentChain(pResIR,pResKS,energyTermsBind);
            //note that this code can easily cause a bug!!! (2/2/2020)
            if((strstr(DES_CHAINS,ChainGetName(pChainI))!=NULL && strstr(DES_CHAINS,ChainGetName(pChainK))==NULL) ||
              (strstr(DES_CHAINS,ChainGetName(pChainI))==NULL && strstr(DES_CHAINS,ChainGetName(pChainK))!=NULL)){
                if(FLAG_PPI){
                  EVOEF_EnergyResidueAndOtherResidueDifferentChain(pResIR,pResKS,energyTermsBind);
                }
            }
          }
        }
      }
    }
  }

  //restore the coordinates
  for(int i=0; i<StructureGetChainCount(pStructure); i++){
    Chain *pChainI=StructureGetChain(pStructure,i);
    for(int ir=0; ir<ChainGetResidueCount(pChainI); ir++){
      Residue *pResIR = ChainGetResidue(pChainI,ir);
      if(pResIR->designSiteType!=Type_ResidueDesignType_Fixed){
        for(int atom1=0; atom1<ResidueGetAtomCount(pResIR); atom1++){
          Atom* pAtom1=ResidueGetAtom(pResIR,atom1);
          if(pAtom1->isBBAtom==FALSE) pAtom1->isXyzValid=TRUE;
        }
      }
    }
  }
  return Success;
}


int StructureCalcSequenceEnergy(Structure* pStructure,Sequence* pSequence,char* programpath,char* prf,char* ss,char* sa,char* phipsi){
  pSequence->etot = 0;
  pSequence->eevo = 0;
  pSequence->ephy = 0;
  pSequence->ebin = 0;
  if(FLAG_EVOLUTION){
    char fasseq[MAX_SEQ_LEN+1]="";
    StructureGetWholeSequence(pStructure,pSequence,fasseq);
    if(FLAG_EVOPHIPSI){
      pSequence->eevo += EvolutionScoreFull(programpath,prf,ss,sa,phipsi,fasseq);
    }
    else{
      pSequence->eevo += EvolutionScoreProfileFromSequence(programpath,fasseq);
      //pSequence->eevo += EvolutionScoreProfileSSSA(programpath,prf,ss,sa,phipsi,fasseq);
    }
  }

  if(FLAG_PHYSICS){
    double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
    double energyTermsBind[MAX_EVOEF_ENERGY_TERM_NUM]={0};
    //calculate the template energy for the conformation-fixed residues
    //get template energy
    StructureCalcSequenceTemplateEnergy(pStructure,pSequence,energyTerms,energyTermsBind);
    //get pairwise rotamer energy
    for(int i=0;i<StructureGetDesignSiteCount(pStructure);i++){
      DesignSite* pSiteI = StructureGetDesignSite(pStructure,i);
      RotamerSet* pSetI = DesignSiteGetRotamers(pSiteI);
      Rotamer* pRotamerI = RotamerSetGet(pSetI,IntArrayGet(&pSequence->rotamerIndices,i));
      RotamerRestore(pRotamerI,pSetI);
      //RotamerShowAtomParameter(pRotamerI);
      //RotamerShowBondInformation(pRotamerI);
      for(int k=i;k<StructureGetDesignSiteCount(pStructure);k++){
        if(k==i){
          pSequence->ephy+=pRotamerI->selfenergy;
          pSequence->ebin+=pRotamerI->selfenergyBind;
        }
        else{
          DesignSite* pSiteK = StructureGetDesignSite(pStructure,k);
          RotamerSet* pSetK = DesignSiteGetRotamers(pSiteK);
          Rotamer* pRotamerK = RotamerSetGet(pSetK,IntArrayGet(&pSequence->rotamerIndices,k));
          RotamerRestore(pRotamerK,pSetK);
          if(pSiteI->chainIndex==pSiteK->chainIndex){
            if(ChainGetType(StructureGetChain(pStructure,pSiteI->chainIndex))==Type_Chain_Protein){
              EVOEF_EnergyProteinRotamerAndRotamerSameChain(pRotamerI,pRotamerK,energyTerms);
            }
          }
          else{//different chains
            if(ChainGetType(StructureGetChain(pStructure,pSiteI->chainIndex))==Type_Chain_SmallMol ||
              ChainGetType(StructureGetChain(pStructure,pSiteK->chainIndex))==Type_Chain_SmallMol){
                double energyTermsPL[MAX_EVOEF_ENERGY_TERM_NUM]={0};
                EVOEF_EnergyLigandRotamerAndRotamer(pRotamerI,pRotamerK,energyTermsPL);
                for(int index=0;index<MAX_EVOEF_ENERGY_TERM_NUM;index++){
                  energyTerms[index]+=energyTermsPL[index];
                }
                if((strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pSiteI->chainIndex)))!=NULL && 
                  strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pSiteK->chainIndex)))==NULL) ||
                  (strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pSiteI->chainIndex)))==NULL && 
                  strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pSiteK->chainIndex)))!=NULL)){
                    if(FLAG_PROT_LIG || FLAG_ENZYME){
                      for(int index=0;index<MAX_EVOEF_ENERGY_TERM_NUM;index++){
                        energyTermsBind[index]+=energyTermsPL[index];
                      }
                    }
                }
            }
            else if(ChainGetType(StructureGetChain(pStructure,pSiteI->chainIndex))==Type_Chain_Protein &&
              ChainGetType(StructureGetChain(pStructure,pSiteK->chainIndex))==Type_Chain_Protein){
                double energyTermsDC[MAX_EVOEF_ENERGY_TERM_NUM]={0};
                EVOEF_EnergyProteinRotamerAndRotamerDifferentChain(pRotamerI,pRotamerK,energyTermsDC);
                for(int index=0;index<MAX_EVOEF_ENERGY_TERM_NUM;index++){
                  energyTerms[index]+=energyTermsDC[index];
                  //energyTermsBind[index]+=energyTermsDC[index];
                }
                //note that this code can easily cause a bug!!! (2/2/2020)
                if((strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pSiteI->chainIndex)))!=NULL && 
                  strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pSiteK->chainIndex)))==NULL) ||
                  (strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pSiteI->chainIndex)))==NULL && 
                  strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pSiteK->chainIndex)))!=NULL)){
                    if(FLAG_PPI){
                      for(int index=0;index<MAX_EVOEF_ENERGY_TERM_NUM;index++){
                        energyTermsBind[index]+=energyTermsDC[index];
                      }
                    }
                }
            }
          }
          RotamerExtract(pRotamerK);
        }
      }
      RotamerExtract(pRotamerI);
    }
    
    EnergyTermWeighting(energyTerms);
    pSequence->ephy += energyTerms[0];
    EnergyTermWeighting(energyTermsBind);
    pSequence->ebin += energyTermsBind[0];
  }

  //normalize by length
  //pSequence->eevo /= sqrt((double)TOT_SEQ_LEN);
  //pSequence->ephy /= sqrt((double)TOT_SEQ_LEN);
  //do not normalize binding energy
  //pSequence->ebin /= sqrt((double)TOT_SEQ_LEN);
  pSequence->etot = WGT_PROFILE*pSequence->eevo + 1.0*pSequence->ephy;
  return Success;
}


/////////////////////////////////////////////////////////////////
//calculate energy change upon a single mutations
/////////////////////////////////////////////////////////////////
int StructureCalcEnergyChangeUponSingleMutation(Structure* pStructure,Sequence* pSequence,int mutSiteIndex,int mutRotIndex,double* dtot,double* dphy,double *devo,char* programpath,char* prf,char* ss, char* sa,char* phipsi){
  double phyBefore = 0;
  double phyAfter = 0;
  double evoBefore=0;
  double evoAfter=0;

  if(FLAG_PHYSICS == TRUE){
    DesignSite* pCurSite=StructureGetDesignSite(pStructure,mutSiteIndex);
    RotamerSet* pCurSet = DesignSiteGetRotamers(pCurSite);
    Rotamer* pCurRotamer = RotamerSetGet(pCurSet,IntArrayGet(&pSequence->rotamerIndices,mutSiteIndex));
    Rotamer* pNewRotamer = RotamerSetGet(pCurSet,mutRotIndex);
    RotamerRestore(pCurRotamer,pCurSet);
    RotamerRestore(pNewRotamer,pCurSet);
    //get self energy (template)
    phyBefore += pCurRotamer->selfenergy;
    phyAfter  += pNewRotamer->selfenergy;
    //get pairwise rotamer energy
    double energyTerms1[MAX_EVOEF_ENERGY_TERM_NUM]={0};
    double energyTerms2[MAX_EVOEF_ENERGY_TERM_NUM]={0};
    for(int i=0;i<pStructure->designSiteCount;i++){
      if(i==mutSiteIndex) continue;
      DesignSite* pSiteI = StructureGetDesignSite(pStructure,i);
      RotamerSet* pSetI = DesignSiteGetRotamers(pSiteI);
      Rotamer* pRotamerI = RotamerSetGet(pSetI,IntArrayGet(&pSequence->rotamerIndices,i));
      RotamerRestore(pRotamerI,pSetI);
      //get energy change
      if(pCurSite->chainIndex == pSiteI->chainIndex){
        if(ChainGetType(StructureGetChain(pStructure,pSiteI->chainIndex))==Type_Chain_Protein){
          EVOEF_EnergyProteinRotamerAndRotamerSameChain(pRotamerI,pCurRotamer,energyTerms1);
          EVOEF_EnergyProteinRotamerAndRotamerSameChain(pRotamerI,pNewRotamer,energyTerms2);
        }
      }
      else{
        if(ChainGetType(StructureGetChain(pStructure,pSiteI->chainIndex))==Type_Chain_SmallMol ||
          ChainGetType(StructureGetChain(pStructure,pCurSite->chainIndex))==Type_Chain_SmallMol){
            EVOEF_EnergyLigandRotamerAndRotamer(pRotamerI,pCurRotamer,energyTerms1);
            EVOEF_EnergyLigandRotamerAndRotamer(pRotamerI,pNewRotamer,energyTerms2);
        }
        else if(ChainGetType(StructureGetChain(pStructure,pSiteI->chainIndex))==Type_Chain_Protein &&
          ChainGetType(StructureGetChain(pStructure,pCurSite->chainIndex))==Type_Chain_Protein){
            EVOEF_EnergyProteinRotamerAndRotamerDifferentChain(pRotamerI,pCurRotamer,energyTerms1);
            EVOEF_EnergyProteinRotamerAndRotamerDifferentChain(pRotamerI,pNewRotamer,energyTerms2);
        }
      }
      RotamerExtract(pRotamerI);
    }
    EnergyTermWeighting(energyTerms1);
    EnergyTermWeighting(energyTerms2);
    phyBefore += energyTerms1[0];
    phyAfter += energyTerms2[0];
    *dphy=phyAfter-phyBefore;
    RotamerExtract(pCurRotamer);
    RotamerExtract(pNewRotamer);
  }

  if(FLAG_EVOLUTION){
    DesignSite* pCurSite=StructureGetDesignSite(pStructure,mutSiteIndex);
    RotamerSet* pCurSet = DesignSiteGetRotamers(pCurSite);
    Rotamer* pCurRotamer = RotamerSetGet(pCurSet,IntArrayGet(&pSequence->rotamerIndices,mutSiteIndex));
    Rotamer* pNewRotamer = RotamerSetGet(pCurSet,mutRotIndex);
    if(RotamerAndRotamerInSameType(pCurRotamer,pNewRotamer)==FALSE){
      char fasseq1[MAX_SEQ_LEN];
      StructureGetWholeSequence(pStructure,pSequence,fasseq1);
      char fasseq2[MAX_SEQ_LEN];
      Sequence newSeq;
      SequenceCreate(&newSeq);
      SequenceCopy(&newSeq,pSequence);
      IntArraySet(&newSeq.rotamerIndices,mutSiteIndex,mutRotIndex);
      StructureGetWholeSequence(pStructure,&newSeq,fasseq2);
      if(FLAG_EVOPHIPSI){
        evoBefore += EvolutionScoreFull(programpath,prf,ss,sa,phipsi,fasseq1);
        evoAfter += EvolutionScoreFull(programpath,prf,ss,sa,phipsi,fasseq2);
      }
      else{
        evoBefore += EvolutionScoreProfileFromSequence(programpath,fasseq1);
        evoAfter += EvolutionScoreProfileFromSequence(programpath,fasseq2);
        //evoBefore += EvolutionScoreProfileSSSA(programpath,prf,ss,sa,phipsi,fasseq1);
        //evoAfter += EvolutionScoreProfileSSSA(programpath,prf,ss,sa,phipsi,fasseq2);
      }
    }
    *devo=evoAfter-evoBefore;
  }

  *dtot = WGT_PROFILE*(*devo) + 1.0*(*dphy);
  return Success;
}


int MetropolisCriteria(Sequence* pOld,Sequence* pBest,Structure* pStructure,RotamerList* pList,int *seqIndex,double temp,int stepCount, FILE* fp,FILE *fp2,char* programpath, char* prf, char* ss, char* sa, char* phipsi){
  for(int i=0; i<stepCount; i++){
    int   mutationSiteIndex;
    int   mutationRotamerIndex;
    SequenceRandomSingleSiteIndex(&mutationSiteIndex,StructureGetDesignSiteCount(pStructure));
    SequenceRandomSingleRotamerIndex(pStructure, pList, mutationSiteIndex,&mutationRotamerIndex);
    double dtot=0, dphy=0, devo=0;
    StructureCalcEnergyChangeUponSingleMutation(pStructure,pOld,mutationSiteIndex,mutationRotamerIndex,&dtot,&dphy,&devo,programpath,prf,ss,sa,phipsi);
    if(exp(-1.0*dtot/temp)>(double)(rand()+1.0)/(RAND_MAX+1.0)){
      SequenceChangeSingleSite(pOld,mutationSiteIndex,mutationRotamerIndex);
      pOld->etot += dtot;
      pOld->ephy += dphy;
      pOld->eevo += devo;
      //printf("step=%8d,temp=%12.6f,dtot=%12.6f, dphy=%12.6f, devo=%12.6f => accepted\n",*seqIndex+i,temp,dtot,dphy,devo);
      //printf("temp=%12.6f,dtot=%12.6f, dphy=%12.6f, devo=%12.6f => accepted\n",temp,pOld->tote,pOld->phye,pOld->evoe);
      //ouput all accepted sequences
      SequenceWriteDesignRotamer(pOld,pStructure,*seqIndex+i,fp);
      SequenceWriteDesignFasta(pOld,pStructure,*seqIndex+i,fp2);
      if(pOld->etot<pBest->etot){
        SequenceCopy(pBest,pOld);
      }
    }
    else{
      //printf("step=%8d,temp=%12.6f,dtot=%12.6f, dphy=%12.6f, devo=%12.6f => rejected\n",*seqIndex+i,temp,dtot,dphy,devo);
    }
  }
  *seqIndex+=stepCount;
  return Success;
}



int SequenceGenerateRandomSeed(Sequence* pThis, RotamerList* pList){
  int i;
  SequenceDestroy(pThis);
  SequenceCreate(pThis);
  pThis->designSiteCount=pList->designSiteCount;
  pThis->etot=0;
  pThis->ephy=0;
  pThis->eevo=0;
  IntArrayResize(&pThis->rotamerIndices, pThis->designSiteCount);
  for(i=0; i<pThis->designSiteCount; i++){
    int j = rand()%pList->remainRotamerCount[i];
    int trueJ;
    RotamerOriginalIndexGet(pList, i, j, &trueJ);
    IntArraySet(&pThis->rotamerIndices, i, trueJ);
  }

  return Success;
}

/******************************************************************************************************/
/* the following two functions:
/* int SequenceGenerateInitialSequenceSeed() and
/* int SequenceGenerateCrystalSequenceSeed()
/* are to be used for debugging
/******************************************************************************************************/

int SequenceGenerateInitialSequenceSeed(Structure* pStructure, RotamerList* pList,Sequence* pSequence){
  SequenceDestroy(pSequence);
  SequenceCreate(pSequence);
  pSequence->designSiteCount=pList->designSiteCount;
  pSequence->etot=0.0;
  IntArrayResize(&pSequence->rotamerIndices, pSequence->designSiteCount);
  for(int i=0; i<pSequence->designSiteCount; i++){
    double minSelfEnergy = 1e8;
    int minSelfEnergyIndex = -1;
    DesignSite* pDesignSite = StructureGetDesignSite(pStructure,i);
    RotamerSet* pRotamerSet = DesignSiteGetRotamers(pDesignSite);
    for(int j=0;j<pList->rotamerCount[i];j++){
      if(pList->remainFlag[i][j]==FALSE) continue;
      Rotamer* pRotamer=RotamerSetGet(pRotamerSet,j);
      if(pRotamer->selfenergy<minSelfEnergy){
        minSelfEnergy = pRotamer->selfenergy;
        minSelfEnergyIndex = j;
      }
    }
    IntArraySet(&pSequence->rotamerIndices, i, minSelfEnergyIndex);
  }

  return Success;
}

int SequenceGenerateCrystalSequenceSeed(Structure* pStructure, RotamerList* pList,Sequence* pSequence){
  SequenceDestroy(pSequence);
  SequenceCreate(pSequence);
  pSequence->designSiteCount=pList->designSiteCount;
  pSequence->etot=0.0;
  IntArrayResize(&pSequence->rotamerIndices, pSequence->designSiteCount);
  for(int i=0; i<pSequence->designSiteCount; i++){
    DesignSite* pDesignSite = StructureGetDesignSite(pStructure,i);
    RotamerSet* pRotamerSet = DesignSiteGetRotamers(pDesignSite);
    if(strcmp(ResidueGetName(pDesignSite->pResidue), "HSD") == 0||strcmp(ResidueGetName(pDesignSite->pResidue), "HSE") == 0){
      IntArraySet(&pSequence->rotamerIndices, i, RotamerSetGetCount(DesignSiteGetRotamers(StructureGetDesignSite(pStructure,i)))-2);
      ;
    }
    else{
      IntArraySet(&pSequence->rotamerIndices, i, RotamerSetGetCount(DesignSiteGetRotamers(StructureGetDesignSite(pStructure,i)))-1);
    }
  }

  return Success;
}


int MonteCarloOptimization(Structure* pStructure,RotamerList* pList,char* programpath,char* prf, char* ss, char* sa, char* phipsi,char* rotindexfile,char* fastafile){
  srand((unsigned int)time(NULL));
  Sequence     oldSequence;
  Sequence     bestSequence;
  int status=Success;
  int seqIndex=0;

  FILE* fp=NULL, *fp2=NULL;
  fp=fopen(rotindexfile,"w");
  fp2=fopen(fastafile,"w");

  int remainCount=0;
  for(int i=0; i<pList->designSiteCount; i++){
    remainCount += pList->remainRotamerCount[i];
  }
  printf("total remained rotamer count: %d\n");
  int maxStep= remainCount*STEP_SCALE > MC_STEP ? remainCount*STEP_SCALE : MC_STEP;
  printf("total monte carlo step: %d\n",maxStep);


  StringArray *pRotTypes=(StringArray*)malloc(sizeof(StringArray)*pList->designSiteCount);
  IntArray *pRotCounts=(IntArray*)malloc(sizeof(IntArray)*pList->designSiteCount);
  DesignSiteGetRotamerTypeAndCount(pList,pStructure,&pRotTypes,&pRotCounts);

  printf("Monte Carlo optimization at fixed temperature T = %f\n",MC_TMIN);
  SequenceCreate(&oldSequence);
  SequenceCreate(&bestSequence);
  SequenceGenerateRandomSeed(&oldSequence,pList);
  //SequenceGenerateInitialSequenceSeed(pStructure,pList,&oldSequence);
  //SequenceGenerateCrystalSequenceSeed(pStructure,pList,&oldSequence);
  StructureCalcSequenceEnergy(pStructure,&oldSequence,programpath,prf,ss,sa,phipsi);
  SequenceCopy(&bestSequence,&oldSequence);
  SequenceWriteDesignRotamer(&bestSequence,pStructure,seqIndex,fp);
  SequenceWriteDesignFasta(&bestSequence,pStructure,seqIndex,fp2);
  seqIndex++;

  //MetropolisCriteria(&oldSequence,&bestSequence,pStructure,pList,&seqIndex,MC_TMIN,maxStep,fp,fp2,deschn,programpath,prf,ss,sa,phipsi);
  MetropolisCriteriaNew(&oldSequence,&bestSequence,pStructure,pList,&pRotTypes,&pRotCounts,&seqIndex,MC_TMIN,maxStep,fp,fp2,programpath,prf,ss,sa,phipsi);
  StructureShowMinEnergyStructure(pStructure,&bestSequence,BEST_STRUCT);
  FILE* fp3=fopen(BEST_SEQ,"w");
  SequenceWriteDesignFasta(&bestSequence,pStructure,0,fp3);
  fclose(fp3);

  fclose(fp);
  fclose(fp2);
  SequenceDestroy(&oldSequence);
  SequenceDestroy(&bestSequence);

  //release memory
  for(int i=0; i<pList->designSiteCount; i++){
    StringArrayDestroy(&pRotTypes[i]);
    IntArrayDestroy(&pRotCounts[i]);
  }
  free(pRotTypes);
  free(pRotCounts);

  return Success;
}

int FlexibleBackboneOptimization(Structure* pStructure,StructureFlex* pStructureFlex,RotamerList* pList,char* DES_CHNS,AtomParamsSet* atomParams,ResiTopoSet *resiTopos,BBdepRotamerLib *bbrotlib,int ***ramaduke, double *accdis,double ***phipsiprob,float ****cancbins,sssegment *sse,int numsse,paircont *bbptb,int numbbptb,int *alphasind,int numsalpha,pairaa *natDistRestr,double **contactCA,double **contactCB){
  const int TOTAL_STEP_NUM=200;
  double newenergy;
  double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
//  double energyTermsBind[MAX_EVOEF_ENERGY_TERM_NUM]={0};

//  EnergyTermInitialize(energyTerms);
  //TODO: expand to multiple chains
  int numchains=1;
  int chains[1];
  StructureFindChainIndex(pStructure,DES_CHNS,&chains[0]);

  for(int i=0;i<StructureGetChainCount(pStructure);i++){
    Chain *pChain = StructureGetChain(pStructure,i);
    ChainFlex *pChainFlex = StructureFlexGetChain(pStructureFlex,i);
//    EnergyHbondnhoc(pChain,pChainFlex,1); //TODO: set to specific protype
    getRamaOutliers(pChain,pChainFlex,ramaduke);
    getAngleOutliers(pChain,pChainFlex);
    getBondLengthOutliers(pChain,pChainFlex);
    getClashes(pChain,pChainFlex);
  }

  char filename2[100] = "flex2.pdb";
  FILE *fileinter2 = fopen(filename2,"w");
  StructureShowInPDBFormat(pStructure, FALSE, fileinter2);
  fclose(fileinter2);




  StructureCalcEnergyFlex(pStructure,bbrotlib,energyTerms);
  EnergyTermWeighting(energyTerms);
//  EnergyTermWeighting(energyTermsBind);
  ChainFlex *pChainFlex = StructureFlexGetChain(pStructureFlex,0);

  double oldenergy = energyTerms[0];
  //oldenergy+=calcConstraintEnergy(pChainTemp,natDistRestr)*10;
  for(int step=0;step<TOTAL_STEP_NUM;step++){
    if(step%100==0) printf("step # %d\n",step);
    makeConformationalMove(pStructure,pStructureFlex,chains,numchains,oldenergy,newenergy,
                           bbrotlib,atomParams,resiTopos,ramaduke,accdis,phipsiprob,cancbins,
                           sse,numsse,bbptb,numbbptb,alphasind,numsalpha,natDistRestr,
                           contactCA,contactCB);
  }

  char filename[100] = "flex.pdb";
  FILE *fileinter = fopen(filename,"w");
  StructureShowInPDBFormat(pStructure, FALSE, fileinter);
  fclose(fileinter);

  return Success;
}


int SimulatedAnnealingOptimization(Structure* pStructure,RotamerList* pList,char* programpath,char* prf, char* ss, char* sa, char* phipsi,char* design_rotindex_file,char* fastafile){
  srand((unsigned int)time(NULL));
  int remainCount=0;
  for(int i=0; i<pList->designSiteCount; i++){
    remainCount += pList->remainRotamerCount[i];
  }
  printf("total remained rotamer count: %d\n", remainCount);

  StringArray *pRotTypes=(StringArray*)malloc(sizeof(StringArray)*pList->designSiteCount);
  IntArray *pRotCounts=(IntArray*)malloc(sizeof(IntArray)*pList->designSiteCount);
  DesignSiteGetRotamerTypeAndCount(pList,pStructure,&pRotTypes,&pRotCounts);

  printf("searching sequences using monte-carlo simulated annealing optimization\n");
  for(int i=1;i<=NTRAJ;i++){
    printf("#search for independent trajectory %d\n",i);
    Sequence     oldSequence;
    Sequence     bestSequence;
    SequenceCreate(&oldSequence);
    SequenceCreate(&bestSequence);
    SequenceGenerateRandomSeed(&oldSequence,pList);
    //SequenceGenerateInitialSequenceSeed(pStructure,pList,&oldSequence);
    //SequenceGenerateCrystalSequenceSeed(pStructure,pList,&oldSequence);
    StructureCalcSequenceEnergy(pStructure,&oldSequence,programpath,prf,ss,sa,phipsi);
    int seqIndex=0;
    FILE* PFILE1=NULL, *PFILE2=NULL;
    PFILE1=fopen(design_rotindex_file,"w");
    PFILE2=fopen(fastafile,"w");
    SequenceCopy(&bestSequence,&oldSequence);
    SequenceWriteDesignRotamer(&bestSequence,pStructure,seqIndex,PFILE1);
    SequenceWriteDesignFasta(&bestSequence,pStructure,seqIndex,PFILE2);
    seqIndex++;
    remainCount = remainCount < SA_STEP ? remainCount:SA_STEP;
    printf("set remained rotamer count to %d for number of monte carlo moves at each temperature\n", remainCount);
    for(int cycle = 1; cycle <= SA_CYCLE; cycle++){
      printf("simulated annealing cycle %3d\n", cycle);
      double t=SA_TMAX/pow(cycle,2.0);
      while(t>SA_TMIN){
        MetropolisCriteriaNew(&oldSequence,&bestSequence,pStructure,pList,&pRotTypes,&pRotCounts,&seqIndex,t,remainCount,PFILE1,PFILE2,programpath,prf,ss,sa,phipsi);
        t *= SA_TDEC;
      }
    }
    fclose(PFILE1);
    fclose(PFILE2);
    //for debug;
    //SequenceCopy(&bestSequence,&oldSequence);
    char BEST_STRUCT_IDX[MAX_LENGTH_FILE_NAME+1];
    sprintf(BEST_STRUCT_IDX,"%s%04d.pdb",BEST_STRUCT,i);
    StructureShowMinEnergyStructure(pStructure,&bestSequence,BEST_STRUCT_IDX);
    if(FLAG_MOL2){
      char BEST_MOL2_IDX[MAX_LENGTH_FILE_NAME+1];
      sprintf(BEST_MOL2_IDX,"%s%04d.mol2",BEST_MOL2,i);
      StructureShowMinEnergyStructureLigand(pStructure,&bestSequence,BEST_MOL2_IDX);
    }
    char BEST_SEQ_IDX[MAX_LENGTH_FILE_NAME+1];
    sprintf(BEST_SEQ_IDX,"%s%04d.txt",BEST_SEQ,i);
    FILE* PFILE3=fopen(BEST_SEQ_IDX,"w");
    SequenceWriteDesignFasta(&bestSequence,pStructure,0,PFILE3);
    fclose(PFILE3);
    SequenceDestroy(&oldSequence);
    SequenceDestroy(&bestSequence);
  }


  //release memory
  for(int i=0; i<pList->designSiteCount; i++){
    StringArrayDestroy(&pRotTypes[i]);
    IntArrayDestroy(&pRotCounts[i]);
  }
  free(pRotTypes);
  free(pRotCounts);

  return Success;
}


int ReannealingOptimization(Structure* pStructure,RotamerList* pList,char* programpath,char* prf, char* ss, char* sa, char* phipsi,char* design_rotindex_file,char* fastafile){
  srand((unsigned int)time(NULL));
  Sequence     oldSequence;
  Sequence     bestSequence;
  int status=Success;
  int seqIndex=0;

  FILE* fp=NULL, *fp2=NULL;
  fp=fopen(design_rotindex_file,"w");
  fp2=fopen(fastafile,"w");

  int remainCount=0;
  for(int i=0; i<pList->designSiteCount; i++){
    remainCount += pList->remainRotamerCount[i];
  }
  printf("total remained rotamer count: %d\n", remainCount);

  StringArray *pRotTypes=(StringArray*)malloc(sizeof(StringArray)*pList->designSiteCount);
  IntArray *pRotCounts=(IntArray*)malloc(sizeof(IntArray)*pList->designSiteCount);
  DesignSiteGetRotamerTypeAndCount(pList,pStructure,&pRotTypes,&pRotCounts);

  printf("searching sequences using re-annealing optimization\n");
  SequenceCreate(&oldSequence);
  SequenceCreate(&bestSequence);
  SequenceGenerateRandomSeed(&oldSequence,pList);
  //SequenceGenerateInitialSequenceSeed(pStructure,pList,&oldSequence);
  //SequenceGenerateCrystalSequenceSeed(pStructure,pList,&oldSequence);
  StructureCalcSequenceEnergy(pStructure,&oldSequence,programpath,prf,ss,sa,phipsi);
  SequenceCopy(&bestSequence,&oldSequence);
  SequenceWriteDesignRotamer(&bestSequence,pStructure,seqIndex,fp);
  SequenceWriteDesignFasta(&bestSequence,pStructure,seqIndex,fp2);
  seqIndex++;

  int totCount = RA_CYCLE*remainCount;
  int naccept = 0;
  for(int i = 0; i < totCount; i++){
    double x = (double)i/totCount;
    double n = 1.0+(double)i/totCount;
    double m = 1.0;
    double miu = 1.0;
    double pai = 3.1415926;
    double t = RA_TMIN + (RA_TMAX - RA_TMIN)*(1-pow(x, miu))*pow(cos(n*pai*x), 2.0*m);
    //MetropolisCriteria(&oldSequence,&bestSequence,pStructure,pList,&seqIndex,t,1,fp,fp2,deschn,programpath,prf,ss,sa,phipsi);
    MetropolisCriteriaNew(&oldSequence,&bestSequence,pStructure,pList,&pRotTypes,&pRotCounts,&seqIndex,t,1,fp,fp2,programpath,prf,ss,sa,phipsi);
  }
  StructureShowMinEnergyStructure(pStructure,&bestSequence,BEST_STRUCT);
  FILE* fp3=fopen(BEST_SEQ,"w");
  SequenceWriteDesignFasta(&bestSequence,pStructure,0,fp3);
  fclose(fp3);

  fclose(fp);
  fclose(fp2);
  SequenceDestroy(&oldSequence);
  SequenceDestroy(&bestSequence);

  //release memory
  for(int i=0; i<pList->designSiteCount; i++){
    StringArrayDestroy(&pRotTypes[i]);
    IntArrayDestroy(&pRotCounts[i]);
  }
  free(pRotTypes);
  free(pRotCounts);

  return Success;
}


int MetropolisCriteriaForSCP(Sequence* pOld,Sequence* pBest,Structure* pStructure,RotamerList* pList,StringArray** ppRotamerType,IntArray** ppRotamerCount,int *seqIndex,double temp,int stepCount, FILE* fp,FILE *fp2,char* programpath, char* prf, char* ss, char* sa, char* phipsi){
  for(int i=0; i<stepCount; i++){
    int   mutationSiteIndex;
    int   mutationRotamerIndex;
    SequenceRandomSingleSiteIndex(&mutationSiteIndex,StructureGetDesignSiteCount(pStructure));
    SequenceRandomSingleRotamerIndexNew(pStructure,pList,ppRotamerType,ppRotamerCount,mutationSiteIndex,&mutationRotamerIndex);
    double dtot=0,dphy=0,devo=0,dbin=0;
    //StructureCalcEnergyChangeUponSingleMutation(pStructure,pOld,mutationSiteIndex,mutationRotamerIndex,&dtot,&dphy,&devo,programpath,prf,ss,sa,phipsi);
    StructureCalcEnergyChangeUponSingleMutationNew(pStructure,pOld,mutationSiteIndex,mutationRotamerIndex,&dtot,&dphy,&dbin,&devo,programpath,prf,ss,sa,phipsi);
    if(exp(-1.0*dtot/temp)>(double)(rand()+1.0)/(RAND_MAX+1.0)){
      SequenceChangeSingleSite(pOld,mutationSiteIndex,mutationRotamerIndex);
      pOld->etot += dtot;
      pOld->ephy += dphy;
      pOld->ebin += dbin;
      pOld->eevo += devo;
      //printf("step=%8d,temp=%12.6f,dtot=%12.6f, dphy=%12.6f, devo=%12.6f => accepted\n",*seqIndex+i,temp,dtot,dphy,devo);
      //printf("temp=%12.6f,dtot=%12.6f, dphy=%12.6f, devo=%12.6f => accepted\n",temp,pOld->tote,pOld->phye,pOld->evoe);
      //ouput all accepted sequences
      SequenceWriteDesignRotamer(pOld,pStructure,*seqIndex+i,fp);
      SequenceWriteDesignFastaForSCP(pOld,pStructure,*seqIndex+i,fp2);
      if(pOld->etot<pBest->etot){
        SequenceCopy(pBest,pOld);
      }
    }
    else{
      //printf("step=%8d,temp=%12.6f,dtot=%12.6f, dphy=%12.6f, devo=%12.6f => rejected\n",*seqIndex+i,temp,dtot,dphy,devo);
    }
  }
  *seqIndex+=stepCount;
  return Success;
}


int SimulatedAnnealingOptimizationForSCP(Structure* pStructure,RotamerList* pList,char* programpath,char* prf,char* ss,char* sa,char* phipsi,char* design_rotindex_file,char* fastafile){
  srand((unsigned int)time(NULL));
  Sequence     oldSequence;
  Sequence     bestSequence;
  int status=Success;
  int seqIndex=0;

  FILE* fp=NULL, *fp2=NULL;
  fp=fopen(design_rotindex_file,"w");
  fp2=fopen(fastafile,"w");

  int remainCount=0;
  for(int i=0; i<pList->designSiteCount; i++){
    remainCount += pList->remainRotamerCount[i];
  }
  remainCount = remainCount < SA_STEP ? remainCount:SA_STEP;
  printf("total remained rotamer count: %d\n", remainCount);

  StringArray *pRotTypes=(StringArray*)malloc(sizeof(StringArray)*pList->designSiteCount);
  IntArray *pRotCounts=(IntArray*)malloc(sizeof(IntArray)*pList->designSiteCount);
  DesignSiteGetRotamerTypeAndCount(pList,pStructure,&pRotTypes,&pRotCounts);

  printf("sidechain repacking using simulated annealing optimization\n");
  SequenceCreate(&oldSequence);
  SequenceCreate(&bestSequence);
  SequenceGenerateRandomSeed(&oldSequence,pList);
  //SequenceGenerateInitialSequenceSeed(pStructure,pList,&oldSequence);
  //SequenceGenerateCrystalSequenceSeed(pStructure,pList,&oldSequence);
  StructureCalcSequenceEnergy(pStructure,&oldSequence,programpath,prf,ss,sa,phipsi);
  SequenceCopy(&bestSequence,&oldSequence);
  SequenceWriteDesignRotamer(&bestSequence,pStructure,seqIndex,fp);
  SequenceWriteDesignFastaForSCP(&bestSequence,pStructure,seqIndex,fp2);
  seqIndex++;

  for(int cycle = 1; cycle <= SA_CYCLE; cycle++){
    printf("simulated annealing cycle %3d\n", cycle);
    double t=SA_TMAX/pow(cycle,2.0);
    while(t>SA_TMIN){
      //MetropolisCriteria(&oldSequence,&bestSequence,pStructure,pList,&seqIndex,t,remainCount,fp,fp2,deschn,programpath,prf,ss,sa,phipsi);
      MetropolisCriteriaForSCP(&oldSequence,&bestSequence,pStructure,pList,&pRotTypes,&pRotCounts,&seqIndex,t,remainCount,fp,fp2,programpath,prf,ss,sa,phipsi);
      t *= SA_TDEC;
    }
  }
  StructureShowMinEnergyStructure(pStructure,&bestSequence,BEST_STRUCT);
  FILE* fp3=fopen(BEST_SEQ,"w");
  SequenceWriteDesignFastaForSCP(&bestSequence,pStructure,0,fp3);
  fclose(fp3);

  //show performance of repacked structure
  //(1)show if each residue is recovered by the rotamer library
  //(2)show the rmsd of each repacked residue to the native
  char FILE_REPACK_TORSION[MAX_LENGTH_FILE_NAME+1];
  char FILE_REPACK_RMSD[MAX_LENGTH_FILE_NAME+1];
  sprintf(FILE_REPACK_TORSION,"%s_scptorsion.txt",PDBID);
  sprintf(FILE_REPACK_RMSD,"%s_scprmsd.txt",PDBID);
  FILE* fp4=fopen(FILE_REPACK_TORSION,"w");
  FILE* fp5=fopen(FILE_REPACK_RMSD,"w");
  SequenceWriteSCPtorsionAndRmsd(&bestSequence,pStructure,fp4,fp5);
  fclose(fp4);
  fclose(fp5);

  fclose(fp);
  fclose(fp2);
  SequenceDestroy(&oldSequence);
  SequenceDestroy(&bestSequence);

  //release memory
  for(int i=0; i<pList->designSiteCount; i++){
    StringArrayDestroy(&pRotTypes[i]);
    IntArrayDestroy(&pRotCounts[i]);
  }
  free(pRotTypes);
  free(pRotCounts);

  return Success;
}


int StructureCalcEnergyChangeUponSingleMutationNew(Structure* pStructure,Sequence* pSequence,int mutSiteIndex,int mutRotIndex,double* dtot,double* dphy,double* dbin,double *devo,char* programpath,char* prf,char* ss, char* sa,char* phipsi){
  double phyBefore = 0;
  double phyAfter = 0;
  double binBefore = 0;
  double binAfter = 0;
  double evoBefore=0;
  double evoAfter=0;

  if(FLAG_PHYSICS == TRUE){
    DesignSite* pCurSite=StructureGetDesignSite(pStructure,mutSiteIndex);
    RotamerSet* pCurSet = DesignSiteGetRotamers(pCurSite);
    Rotamer* pCurRotamer = RotamerSetGet(pCurSet,IntArrayGet(&pSequence->rotamerIndices,mutSiteIndex));
    Rotamer* pNewRotamer = RotamerSetGet(pCurSet,mutRotIndex);
    RotamerRestore(pCurRotamer,pCurSet);
    RotamerRestore(pNewRotamer,pCurSet);
    //get self energy (template)
    phyBefore += pCurRotamer->selfenergy;
    phyAfter  += pNewRotamer->selfenergy;
    binBefore += pCurRotamer->selfenergyBind;
    binAfter  += pNewRotamer->selfenergyBind;
    //get pairwise rotamer energy
    double energyTerms1[MAX_EVOEF_ENERGY_TERM_NUM]={0};
    double energyTerms2[MAX_EVOEF_ENERGY_TERM_NUM]={0};
    double energyTermsBind1[MAX_EVOEF_ENERGY_TERM_NUM]={0};
    double energyTermsBind2[MAX_EVOEF_ENERGY_TERM_NUM]={0};
    for(int i=0;i<pStructure->designSiteCount;i++){
      if(i==mutSiteIndex) continue;
      DesignSite* pSiteI = StructureGetDesignSite(pStructure,i);
      RotamerSet* pSetI = DesignSiteGetRotamers(pSiteI);
      Rotamer* pRotamerI = RotamerSetGet(pSetI,IntArrayGet(&pSequence->rotamerIndices,i));
      RotamerRestore(pRotamerI,pSetI);
      //get energy change
      if(pCurSite->chainIndex==pSiteI->chainIndex){
        if(ChainGetType(StructureGetChain(pStructure,pSiteI->chainIndex))==Type_Chain_Protein){
          EVOEF_EnergyProteinRotamerAndRotamerSameChain(pRotamerI,pCurRotamer,energyTerms1);
          EVOEF_EnergyProteinRotamerAndRotamerSameChain(pRotamerI,pNewRotamer,energyTerms2);
        }
      }
      else{
        if(ChainGetType(StructureGetChain(pStructure,pSiteI->chainIndex))==Type_Chain_SmallMol ||
          ChainGetType(StructureGetChain(pStructure,pCurSite->chainIndex))==Type_Chain_SmallMol){
            //EVOEF_EnergyLigandRotamerAndRotamer(pRotamerI,pCurRotamer,energyTerms1);
            //EVOEF_EnergyLigandRotamerAndRotamer(pRotamerI,pNewRotamer,energyTerms2);
            double energyTermsPL1[MAX_EVOEF_ENERGY_TERM_NUM]={0};
            double energyTermsPL2[MAX_EVOEF_ENERGY_TERM_NUM]={0};
            EVOEF_EnergyLigandRotamerAndRotamer(pRotamerI,pCurRotamer,energyTermsPL1);
            EVOEF_EnergyLigandRotamerAndRotamer(pRotamerI,pNewRotamer,energyTermsPL2);
            for(int index=0;index<MAX_EVOEF_ENERGY_TERM_NUM;index++){
              energyTerms1[index]+=energyTermsPL1[index];
              energyTerms2[index]+=energyTermsPL2[index];
            }
            if((strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pSiteI->chainIndex)))!=NULL && 
              strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pCurSite->chainIndex)))==NULL) ||
              (strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pSiteI->chainIndex)))==NULL && 
              strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pCurSite->chainIndex)))!=NULL)){
                if(FLAG_PROT_LIG || FLAG_ENZYME){
                  for(int index=0;index<MAX_EVOEF_ENERGY_TERM_NUM;index++){
                    energyTermsBind1[index]+=energyTermsPL1[index];
                    energyTermsBind2[index]+=energyTermsPL2[index];
                  }
                }
            }
        }
        else if(ChainGetType(StructureGetChain(pStructure,pSiteI->chainIndex))==Type_Chain_Protein &&
          ChainGetType(StructureGetChain(pStructure,pCurSite->chainIndex))==Type_Chain_Protein){
            //EVOEF_EnergyProteinRotamerAndRotamerDifferentChain(pRotamerI,pCurRotamer,energyTerms1);
            //EVOEF_EnergyProteinRotamerAndRotamerDifferentChain(pRotamerI,pNewRotamer,energyTerms2);
            double energyTermsDC1[MAX_EVOEF_ENERGY_TERM_NUM]={0};
            double energyTermsDC2[MAX_EVOEF_ENERGY_TERM_NUM]={0};
            EVOEF_EnergyProteinRotamerAndRotamerDifferentChain(pRotamerI,pCurRotamer,energyTermsDC1);
            EVOEF_EnergyProteinRotamerAndRotamerDifferentChain(pRotamerI,pNewRotamer,energyTermsDC2);
            for(int index=0;index<MAX_EVOEF_ENERGY_TERM_NUM;index++){
              energyTerms1[index]+=energyTermsDC1[index];
              energyTerms2[index]+=energyTermsDC2[index];
            }
            if((strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pSiteI->chainIndex)))!=NULL && 
              strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pCurSite->chainIndex)))==NULL) ||
              (strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pSiteI->chainIndex)))==NULL && 
              strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pCurSite->chainIndex)))!=NULL)){
                if(FLAG_PPI){
                  for(int index=0;index<MAX_EVOEF_ENERGY_TERM_NUM;index++){
                    energyTermsBind1[index]+=energyTermsDC1[index];
                    energyTermsBind2[index]+=energyTermsDC2[index];
                  }
                }
            }
        }
      }
      RotamerExtract(pRotamerI);
    }
    EnergyTermWeighting(energyTerms1);
    EnergyTermWeighting(energyTerms2);
    EnergyTermWeighting(energyTermsBind1);
    EnergyTermWeighting(energyTermsBind2);
    phyBefore += energyTerms1[0];
    phyAfter  += energyTerms2[0];
    *dphy=phyAfter-phyBefore;
    binBefore += energyTermsBind1[0];
    binAfter  += energyTermsBind2[0];
    *dbin=binAfter-binBefore;
    RotamerExtract(pCurRotamer);
    RotamerExtract(pNewRotamer);
  }

  if(FLAG_EVOLUTION){
    DesignSite* pCurSite=StructureGetDesignSite(pStructure,mutSiteIndex);
    RotamerSet* pCurSet = DesignSiteGetRotamers(pCurSite);
    Rotamer* pCurRotamer = RotamerSetGet(pCurSet,IntArrayGet(&pSequence->rotamerIndices,mutSiteIndex));
    Rotamer* pNewRotamer = RotamerSetGet(pCurSet,mutRotIndex);
    if(RotamerAndRotamerInSameType(pCurRotamer,pNewRotamer)==FALSE){
      char fasseq1[MAX_SEQ_LEN];
      StructureGetWholeSequence(pStructure,pSequence,fasseq1);
      char fasseq2[MAX_SEQ_LEN];
      Sequence newSeq;
      SequenceCreate(&newSeq);
      SequenceCopy(&newSeq,pSequence);
      IntArraySet(&newSeq.rotamerIndices,mutSiteIndex,mutRotIndex);
      StructureGetWholeSequence(pStructure,&newSeq,fasseq2);
      if(FLAG_EVOPHIPSI){
        evoBefore += EvolutionScoreFull(programpath,prf,ss,sa,phipsi,fasseq1);
        evoAfter += EvolutionScoreFull(programpath,prf,ss,sa,phipsi,fasseq2);
      }
      else{
        evoBefore += EvolutionScoreProfileFromSequence(programpath,fasseq1);
        evoAfter += EvolutionScoreProfileFromSequence(programpath,fasseq2);
        //evoBefore += EvolutionScoreProfileSSSA(programpath,prf,ss,sa,phipsi,fasseq1);
        //evoAfter += EvolutionScoreProfileSSSA(programpath,prf,ss,sa,phipsi,fasseq2);
      }
      SequenceDestroy(&newSeq);
    }
    *devo=evoAfter-evoBefore;
  }

  *dtot = WGT_PROFILE*(*devo) + 1.0*(*dphy);
  return Success;
}


int MetropolisCriteriaNew(Sequence* pOld,Sequence* pBest,Structure* pStructure,RotamerList* pList,StringArray** ppRotamerType,IntArray** ppRotamerCount,int *seqIndex,double temp,int stepCount, FILE* fp,FILE *fp2,char* programpath, char* prf, char* ss, char* sa, char* phipsi){
  for(int i=0; i<stepCount; i++){
    int   mutationSiteIndex;
    int   mutationRotamerIndex;
    SequenceRandomSingleSiteIndex(&mutationSiteIndex,StructureGetDesignSiteCount(pStructure));
    SequenceRandomSingleRotamerIndexNew(pStructure,pList,ppRotamerType,ppRotamerCount,mutationSiteIndex,&mutationRotamerIndex);
    double dtot=0,dphy=0,devo=0,dbin=0;
    //StructureCalcEnergyChangeUponSingleMutation(pStructure,pOld,mutationSiteIndex,mutationRotamerIndex,&dtot,&dphy,&devo,deschn,programpath,prf,ss,sa,phipsi);
    StructureCalcEnergyChangeUponSingleMutationNew(pStructure,pOld,mutationSiteIndex,mutationRotamerIndex,&dtot,&dphy,&dbin,&devo,programpath,prf,ss,sa,phipsi);
    if(exp(-1.0*dtot/temp)>(double)(rand()+1.0)/(RAND_MAX+1.0)){
      SequenceChangeSingleSite(pOld,mutationSiteIndex,mutationRotamerIndex);
      pOld->etot += dtot;
      pOld->ephy += dphy;
      pOld->ebin += dbin;
      pOld->eevo += devo;
      //printf("step=%8d,temp=%12.6f,dtot=%12.6f, dphy=%12.6f, devo=%12.6f => accepted\n",*seqIndex+i,temp,dtot,dphy,devo);
      //printf("temp=%12.6f,dtot=%12.6f, dphy=%12.6f, devo=%12.6f => accepted\n",temp,pOld->tote,pOld->phye,pOld->evoe);
      //ouput all accepted sequences
      SequenceWriteDesignRotamer(pOld,pStructure,*seqIndex+i,fp);
      SequenceWriteDesignFasta(pOld,pStructure,*seqIndex+i,fp2);
      if(pOld->etot<pBest->etot){
        SequenceCopy(pBest,pOld);
      }
    }
    /*else{
      printf("step=%8d,temp=%12.6f,dtot=%12.6f, dphy=%12.6f, devo=%12.6f => rejected\n",*seqIndex+i,temp,dtot,dphy,devo);
    }*/
  }
  *seqIndex+=stepCount;
  return Success;
}



