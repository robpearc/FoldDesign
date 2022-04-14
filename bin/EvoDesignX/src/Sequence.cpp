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

#include "Sequence.h"
#include "Structure.h"
#include "EnergyMatrix.h"
#include "EnergyFunction.h"
#include <time.h>
#include <string.h>
#include <ctype.h>

extern double TORSION_DEVIATION_CUTOFF;
extern char MOL2[MAX_LENGTH_FILE_NAME+1];
extern char DES_CHAINS[10];

int SequenceCreate(Sequence* pThis){
  pThis->designSiteCount = 0;
  pThis->etot = 0;
  pThis->eevo = 0;
  pThis->ephy = 0;
  pThis->ebin = 0;
  pThis->ebpf = 0;
  IntArrayCreate(&pThis->rotamerIndices, 0);
  return Success;
}

int SequenceDestroy(Sequence* pThis){
  pThis->designSiteCount = 0;
  pThis->etot = 0;
  pThis->eevo = 0;
  pThis->ephy = 0;
  pThis->ebin = 0;
  pThis->ebpf = 0;
  IntArrayDestroy(&pThis->rotamerIndices);
  return Success;
}


int SequenceCopy(Sequence* pThis, Sequence* pOther){
  pThis->designSiteCount = pOther->designSiteCount;
  pThis->etot = pOther->etot;
  pThis->ephy = pOther->ephy;
  pThis->ebin = pOther->ebin;
  pThis->eevo = pOther->eevo;
  pThis->ebpf = pOther->ebpf;
  IntArrayCopy(&pThis->rotamerIndices, &pOther->rotamerIndices);
  return Success;
}


int SequenceWriteDesignRotamer(Sequence* pThis,Structure*pStructure,int index,FILE* pFile){
  BOOL firstchain=TRUE;
  int flagIndex=-1;
  for(int i=0; i<StructureGetChainCount(pStructure);i++){
    Chain* pChain=StructureGetChain(pStructure,i);
    for(int j=0;j<ChainGetResidueCount(pChain);j++){
      Residue* pResi=ChainGetResidue(pChain,j);
      char res=ThreeLetterAAToOneLetterAA(ResidueGetName(pResi));
      DesignSite* pDesignSite=StructureFindDesignSite(pStructure,i,j);
      if(pDesignSite==NULL)continue;
      int designsiteIndex=DesignSiteIndexGet(pStructure,i,j);
      RotamerSet* pRotamerSet = DesignSiteGetRotamers(pDesignSite);
      int rotamerIndex = IntArrayGet(&pThis->rotamerIndices,designsiteIndex);
      Rotamer* pRotamer = RotamerSetGet(pRotamerSet,rotamerIndex);
      char rot=ThreeLetterAAToOneLetterAA(RotamerGetType(pRotamer));
      if(i!=flagIndex){
        if(firstchain==TRUE){
          fprintf(pFile, "%d",rotamerIndex);
          firstchain=FALSE;
        }
        else{
          fprintf(pFile, ";%d",rotamerIndex);
        }
        flagIndex=i;
      }
      else{
        fprintf(pFile, ",%d",rotamerIndex);
      }
    }
  }
  fprintf(pFile, " %d\n",index);

  return Success;
}


int SequenceWriteDesignFasta(Sequence* pThis,Structure*pStructure,int index,FILE* pFile2){
  BOOL firstchain=TRUE;
  int chainFlag=-1;
  int totResCount=0;
  int samResCount=0;
  for(int i=0; i<StructureGetChainCount(pStructure);i++){
    Chain* pChain=StructureGetChain(pStructure,i);
    if(strstr(DES_CHAINS,ChainGetName(pChain))==NULL) continue;
    for(int j=0;j<ChainGetResidueCount(pChain);j++){
      Residue* pResi=ChainGetResidue(pChain,j);
      char res=ThreeLetterAAToOneLetterAA(ResidueGetName(pResi));
      DesignSite* pDesignSite=StructureFindDesignSite(pStructure,i,j);
      if(pDesignSite==NULL){
        if(i!=chainFlag){
          if(firstchain==TRUE){
            if(ChainGetType(pChain)==Type_Chain_Protein){
              fprintf(pFile2, "%c",res);
            }
            firstchain=FALSE;
          }
          else{
            if(ChainGetType(pChain)==Type_Chain_Protein){
              fprintf(pFile2, ";%c",res);
            }
          }
          chainFlag=i;
        }
        else{
          if(ChainGetType(pChain)==Type_Chain_Protein){
            fprintf(pFile2, "%c",res);
          }
        }
        continue;
      }

      int designsiteIndex=DesignSiteIndexGet(pStructure,i,j);
      RotamerSet* pRotamerSet = DesignSiteGetRotamers(pDesignSite);
      int rotamerIndex = IntArrayGet(&pThis->rotamerIndices,designsiteIndex);
      Rotamer* pRotamer = RotamerSetGet(pRotamerSet,rotamerIndex);
      char rot=ThreeLetterAAToOneLetterAA(RotamerGetType(pRotamer));

      if(i!=chainFlag){
        if(firstchain==TRUE){
          if(ChainGetType(pChain)==Type_Chain_Protein){
            fprintf(pFile2, "%c",rot);
          }
          firstchain=FALSE;
        }
        else{
          if(ChainGetType(pChain)==Type_Chain_Protein){
            fprintf(pFile2, ";%c",rot);
          }
        }
        chainFlag=i;
      }
      else{
        if(ChainGetType(pChain)==Type_Chain_Protein){
          fprintf(pFile2, "%c",rot);
        }
      }

      if(ChainGetType(pChain)==Type_Chain_Protein && pResi->designSiteType==Type_ResidueDesignType_Mutated){
        if(rot==res){
          samResCount++;
        }
        totResCount++;
      }
    }
  }
  double seqid=(totResCount>0 ? (double)samResCount/totResCount : 1.0);
  fprintf(pFile2, " %12d %13.6f %13.6f %13.6f %13.6f %13.6f\n",index,seqid,pThis->etot,pThis->eevo,pThis->ephy, pThis->ebin);

  return Success;
}


int SequenceWriteDesignFastaForSCP(Sequence* pThis,Structure*pStructure,int index,FILE* pFile2){
  BOOL firstchain=TRUE;
  int chainFlag=-1;
  int nonAG=0;
  int simSidechain=0;
  double rmsd=0;
  for(int i=0; i<StructureGetChainCount(pStructure);i++){
    Chain* pChain=StructureGetChain(pStructure,i);
    if(strstr(DES_CHAINS,ChainGetName(pChain))==NULL) continue;
    for(int j=0;j<ChainGetResidueCount(pChain);j++){
      Residue* pResi=ChainGetResidue(pChain,j);
      char res=ThreeLetterAAToOneLetterAA(ResidueGetName(pResi));
      DesignSite* pDesignSite=StructureFindDesignSite(pStructure,i,j);
      if(pDesignSite==NULL){
        if(i!=chainFlag){
          if(firstchain==TRUE){
            if(ChainGetType(pChain)==Type_Chain_Protein){
              fprintf(pFile2, "%c",res);
            }
            firstchain=FALSE;
          }
          else{
            if(ChainGetType(pChain)==Type_Chain_Protein){
              fprintf(pFile2, ";%c",res);
            }
          }
          chainFlag=i;
        }
        else{
          if(ChainGetType(pChain)==Type_Chain_Protein){
            fprintf(pFile2, "%c",res);
          }
        }
        continue;
      }

      int designsiteIndex=DesignSiteIndexGet(pStructure,i,j);
      RotamerSet* pRotamerSet = DesignSiteGetRotamers(pDesignSite);
      int rotamerIndex = IntArrayGet(&pThis->rotamerIndices,designsiteIndex);
      Rotamer* pRotamer = RotamerSetGet(pRotamerSet,rotamerIndex);
      char rot=ThreeLetterAAToOneLetterAA(RotamerGetType(pRotamer));
      RotamerRestore(pRotamer,pRotamerSet);
      if(i!=chainFlag){
        if(firstchain==TRUE){
          if(ChainGetType(pChain)==Type_Chain_Protein){
            fprintf(pFile2, "%c",rot);
          }
          firstchain=FALSE;
        }
        else{
          if(ChainGetType(pChain)==Type_Chain_Protein){
            fprintf(pFile2, ";%c",rot);
          }
        }
        chainFlag=i;
      }
      else{
        if(ChainGetType(pChain)==Type_Chain_Protein){
          fprintf(pFile2, "%c",rot);
        }
      }

      if(ChainGetType(pChain)==Type_Chain_Protein){
        if(res != 'A' || res != 'G'){
          rmsd += RotamerAndResidueSidechainRMSD(pRotamer,pResi);
          nonAG++;
          if(RotamerAndResidueWithSimilarTorsions(pRotamer,pResi,TORSION_DEVIATION_CUTOFF)==TRUE){
            simSidechain++;
          }
        }
      }
      RotamerExtract(pRotamer);
    }
  }
  double torid=(double)simSidechain/nonAG;
  fprintf(pFile2, " %12d %13.6f %13.6f %13.6f %13.6f %13.6f %13.6f\n",index,torid,rmsd/nonAG,pThis->etot,pThis->eevo,pThis->ephy, pThis->ebin);

  return Success;
}




int StructureShowMinEnergyStructure(Structure* pStructure,Sequence* sequence,char* pdbfile){
  int status=Success;
  char errMsg[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  FILE* pf = fopen(pdbfile,"w");
  if(pf==NULL){
    status = IOError;
    sprintf(errMsg,"in file %s function %s() line %d, cannot open file %s for writing", __FILE__, __FUNCTION__, __LINE__, pdbfile);
    TraceError(errMsg,status);
    return status;
  }
  //output by chain and residue position instead of design site
  for(int i=0; i<StructureGetChainCount(pStructure);i++){
    Chain* pChain=StructureGetChain(pStructure,i);
    for(int j=0;j<ChainGetResidueCount(pChain);j++){
      Residue* pResi=ChainGetResidue(pChain,j);
      if(pResi->designSiteType==Type_ResidueDesignType_Fixed){
        if(pChain->type==Type_Chain_Protein){
          ResidueShowInPDBFormat(pResi,"ATOM",ResidueGetChainName(pResi),1,ResidueGetPosInChain(pResi),TRUE,pf);
        }
        /*else{
          ResidueShowInPDBFormat(pResi,"HETATM",ResidueGetChainName(pResi),1,ResidueGetPosInChain(pResi),FALSE,pf);
        }*/
      }
      else{
        for(int k=0;k<sequence->designSiteCount;k++){
          DesignSite* pDesignSite = StructureGetDesignSite(pStructure,k);
          RotamerSet* pRotamerSet = DesignSiteGetRotamers(pDesignSite);
          int rotamerIndex = IntArrayGet(&sequence->rotamerIndices,k);
          Rotamer* pOptRotamer = RotamerSetGet(pRotamerSet,rotamerIndex);
          if(pDesignSite->chainIndex==i && pDesignSite->resiIndex==j){
            RotamerRestore(pOptRotamer,pRotamerSet);
            if(pChain->type==Type_Chain_Protein){
              RotamerShowInPDBFormat(pOptRotamer,"ATOM",RotamerGetChainName(pOptRotamer),1,RotamerGetPosInChain(pOptRotamer),TRUE,pf);
            }
            /*else{
              RotamerShowInPDBFormat(pOptRotamer,"HETATM",RotamerGetChainName(pOptRotamer),1,RotamerGetPosInChain(pOptRotamer),FALSE,pf);
            }*/
            RotamerExtract(pOptRotamer);
            break;
          }
        }
      }
    }
  }
  fclose(pf);
  return status;
}


int StructureShowMinEnergyStructureLigand(Structure* pStructure, Sequence* sequence, char* mol2file){
  int status=Success;
  char errMsg[MAX_LENGTH_ONE_LINE_IN_FILE+1];

  StringArray lines;
  StringArrayCreate(&lines);
  FILE* pfile=fopen(MOL2,"r");
  char buffer[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  int start=-1,end=-1;
  int counter=0;
  while(fgets(buffer,MAX_LENGTH_ONE_LINE_IN_FILE,pfile)){
    StringArrayAppend(&lines,buffer);
    if(strstr(buffer,"@<TRIPOS>ATOM")!=NULL){
      start=counter+1;
    }
    else if(strstr(buffer,"@<TRIPOS>BOND")!=NULL){
      end=counter-1;
    }
    counter++;
  }
  fclose(pfile);

  FILE* pf = fopen(mol2file,"w");
  if(pf==NULL){
    status = IOError;
    sprintf(errMsg,"in file %s function %s() line %d, cannot open file %s for writing", __FILE__, __FUNCTION__, __LINE__, mol2file);
    TraceError(errMsg,status);
    return status;
  }
  //output by chain and residue position instead of design site
  for(int i=0; i<StructureGetChainCount(pStructure);i++){
    Chain* pChain=StructureGetChain(pStructure,i);
    for(int j=0;j<ChainGetResidueCount(pChain);j++){
      Residue* pResi=ChainGetResidue(pChain,j);
      if(pResi->designSiteType==Type_ResidueDesignType_Fixed){
        if(pChain->type==Type_Chain_SmallMol){
          for(int l=0;l<StringArrayGetCount(&lines);l++){
            if(l<start || l>end){
              fprintf(pf,"%s",StringArrayGet(&lines,l));
            }
            else{
              char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
              strcpy(line,StringArrayGet(&lines,l));
              StringArray words;
              StringArrayCreate(&words);
              StringArraySplitString(&words,line,' ');
              XYZ xyz=ResidueGetAtom(pResi,l-start)->xyz;
              fprintf(pf,"%7s %-10s %10.3f %10.3f %10.3f %-10s %s %s %10s\n",
                StringArrayGet(&words,0),
                StringArrayGet(&words,1),
                xyz.X,
                xyz.Y,
                xyz.Z,
                StringArrayGet(&words,5),
                StringArrayGet(&words,6),
                StringArrayGet(&words,7),
                StringArrayGet(&words,8));
              StringArrayDestroy(&words);
            }
          }
        }
      }
      else{
        for(int k=0;k<sequence->designSiteCount;k++){
          DesignSite* pDesignSite = StructureGetDesignSite(pStructure,k);
          RotamerSet* pRotamerSet = DesignSiteGetRotamers(pDesignSite);
          int rotamerIndex = IntArrayGet(&sequence->rotamerIndices,k);
          Rotamer* pOptRotamer = RotamerSetGet(pRotamerSet,rotamerIndex);
          if(pDesignSite->chainIndex==i && pDesignSite->resiIndex==j){
            RotamerRestore(pOptRotamer,pRotamerSet);
            if(pChain->type==Type_Chain_SmallMol){
              for(int l=0;l<StringArrayGetCount(&lines);l++){
                if(l<start || l>end){
                  fprintf(pf,"%s",StringArrayGet(&lines,l));
                }
                else{
                  char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
                  strcpy(line,StringArrayGet(&lines,l));
                  StringArray words;
                  StringArrayCreate(&words);
                  StringArraySplitString(&words,line,' ');
                  XYZ xyz=RotamerGetAtom(pOptRotamer,l-start)->xyz;
                  fprintf(pf,"%7s %-10s %10.3f %10.3f %10.3f %-10s %s %s %10s\n",
                    StringArrayGet(&words,0),
                    StringArrayGet(&words,1),
                    xyz.X,
                    xyz.Y,
                    xyz.Z,
                    StringArrayGet(&words,5),
                    StringArrayGet(&words,6),
                    StringArrayGet(&words,7),
                    StringArrayGet(&words,8));
                  StringArrayDestroy(&words);
                }
              }
            }
            RotamerExtract(pOptRotamer);
            break;
          }
        }
      }
    }
  }
  fclose(pf);
  StringArrayDestroy(&lines);

  return status;
}



int StructureGetWholeSequence(Structure* pStructure,Sequence* sequence,char* fastaseq){
  //output by chain and residue position instead of design site
  int index=0;
  for(int i=0; i<StructureGetChainCount(pStructure);i++){
    Chain* pChain=StructureGetChain(pStructure,i);
    if(strcmp(DES_CHAINS,ChainGetName(pChain))!=0)continue;
    for(int j=0;j<ChainGetResidueCount(pChain);j++){
      Residue* pResi=ChainGetResidue(pChain,j);
      if(pResi->designSiteType==Type_ResidueDesignType_Fixed){
        fastaseq[index++]=ThreeLetterAAToOneLetterAA(ResidueGetName(pResi));
      }
      else{
        for(int k=0;k<sequence->designSiteCount;k++){
          DesignSite* pDesignSite = StructureGetDesignSite(pStructure,k);
          RotamerSet* pRotamerSet = DesignSiteGetRotamers(pDesignSite);
          int rotamerIndex = IntArrayGet(&sequence->rotamerIndices,k);
          Rotamer* pOptRotamer = RotamerSetGet(pRotamerSet,rotamerIndex);
          if(pDesignSite->chainIndex==i && pDesignSite->resiIndex==j){
            fastaseq[index++]=ThreeLetterAAToOneLetterAA(RotamerGetType(pOptRotamer));
            break;
          }
        }
      }
    }
  }
  fastaseq[index]='\0';
  return Success;
}


int StructureShowDesignedResiduesOnly(Structure* pStructure, Sequence* sequence, char* pdbfile){
  int status=Success;
  char errMsg[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  FILE* pf = fopen(pdbfile,"w");
  if(pf==NULL){
    status = IOError;
    sprintf(errMsg,"in file %s function %s() line %d, cannot open file %s for writing", __FILE__, __FUNCTION__, __LINE__, pdbfile);
    TraceError(errMsg,status);
    return status;
  }
  //output by chain and residue position instead of design site
  for(int i=0; i<StructureGetChainCount(pStructure);i++){
    Chain* pChain=StructureGetChain(pStructure,i);
    for(int j=0;j<ChainGetResidueCount(pChain);j++){
      Residue* pResi=ChainGetResidue(pChain,j);
      if(pResi->designSiteType==Type_ResidueDesignType_Fixed) continue;
      else{
        for(int k=0;k<sequence->designSiteCount;k++){
          DesignSite* pDesignSite = StructureGetDesignSite(pStructure,i);
          RotamerSet* pRotamerSet = DesignSiteGetRotamers(pDesignSite);
          int rotamerIndex = IntArrayGet(&sequence->rotamerIndices,i);
          Rotamer* pOptRotamer = RotamerSetGet(pRotamerSet,rotamerIndex);
          if(pDesignSite->chainIndex==i && pDesignSite->resiIndex==j){
            RotamerRestore(pOptRotamer,pRotamerSet);
            RotamerShowInPDBFormat(pOptRotamer,"ATOM",RotamerGetChainName(pOptRotamer),1,RotamerGetPosInChain(pOptRotamer),FALSE,pf);
            RotamerExtract(pOptRotamer);
            break;
          }
        }
      }
    }
  }
  fclose(pf);
  return status;
}

int SequenceWriteSCPtorsionAndRmsd(Sequence* pThis,Structure*pStructure,FILE* pTorsion,FILE* pRmsd){
  for(int i=0; i<StructureGetChainCount(pStructure);i++){
    Chain* pChain=StructureGetChain(pStructure,i);
    if(strstr(DES_CHAINS,ChainGetName(pChain))==NULL) continue;
    for(int j=0;j<ChainGetResidueCount(pChain);j++){
      Residue* pResi=ChainGetResidue(pChain,j);
      char res=ThreeLetterAAToOneLetterAA(ResidueGetName(pResi));
      DesignSite* pDesignSite=StructureFindDesignSite(pStructure,i,j);
      if(pDesignSite==NULL){
        continue;
      }

      int designsiteIndex=DesignSiteIndexGet(pStructure,i,j);
      RotamerSet* pRotamerSet = DesignSiteGetRotamers(pDesignSite);
      int rotamerIndex = IntArrayGet(&pThis->rotamerIndices,designsiteIndex);
      Rotamer* pRotamer = RotamerSetGet(pRotamerSet,rotamerIndex);
      char rot=ThreeLetterAAToOneLetterAA(RotamerGetType(pRotamer));
      RotamerRestore(pRotamer,pRotamerSet);

      if(ChainGetType(pChain)==Type_Chain_Protein){
        if(res != 'A' && res != 'G'){
          double rmsd = RotamerAndResidueSidechainRMSD(pRotamer,pResi);
          if(RotamerIsSymmetricalCheck(pRotamer)==TRUE){
            Rotamer tempRot;
            RotamerCreate(&tempRot);
            SymmetricalRotamerGenerate(&tempRot,pRotamer);
            double rmsd2=RotamerAndResidueSidechainRMSD(&tempRot,pResi);
            if(rmsd2<rmsd) rmsd=rmsd2;
            RotamerDestroy(&tempRot);
          }
          fprintf(pRmsd, "%s %s %f\n",ChainGetName(pChain),ResidueGetName(pResi),rmsd);
          if(RotamerAndResidueWithSimilarTorsions(pRotamer,pResi,TORSION_DEVIATION_CUTOFF)==TRUE){
            fprintf(pTorsion, "%s %s 1\n",ChainGetName(pChain),ResidueGetName(pResi));
          }
          else{
            fprintf(pTorsion, "%s %s 0\n",ChainGetName(pChain),ResidueGetName(pResi));
          }
        }
        else{
          double rmsd=0.0;
          fprintf(pRmsd, "%s %s %f\n",ChainGetName(pChain),ResidueGetName(pResi),rmsd);
          fprintf(pTorsion, "%s %s 1\n",ChainGetName(pChain),ResidueGetName(pResi));
        }
      }
      RotamerExtract(pRotamer);
    }
  }
  return Success;
}
