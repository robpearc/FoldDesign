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

#include "ProgramFunction.h"
#include "RotamerBuilder.h"
#include "RotamerOptimizer.h"
#include "Evolution.h"
#include <string.h>
#include <ctype.h>

extern BOOL FLAG_MONOMER;
extern BOOL FLAG_PPI;
extern BOOL FLAG_ENZYME;
extern BOOL FLAG_PROT_LIG;

extern BOOL FLAG_ADD_CRYSTAL_ROT;
extern BOOL FLAG_EXPAND_HYDROXYL_ROT;
extern double TORSION_DEVIATION_CUTOFF;

extern double ROT_PROB_CUT_MIN;
extern char PROGRAM_NAME[MAX_LENGTH_FILE_NAME+1];
extern char PROGRAM_VERSION[MAX_LENGTH_FILE_NAME+1];

int EVOEF_help(){
  printf(  
    "program options:\n\n"
    "Basic OPTIONS:\n"
    "Usage: EvoEF [OPTIONS] --pdb=pdbfile\n"  
    "   --version       print version\n"
    "   --help          print help message\n"  
    "   --command=arg   choose your computation type:\n"
    "                   RepairStructure"
    "                   ComputeStability\n"
    "                   ComputeBinding\n"
    "                   BuildMutant\n"
    );
  return Success;
}

int EVOEF_version(){
  printf("program version %s\n",PROGRAM_VERSION);
  return Success;
}


int EVOEF_interface(){
  printf(
    "############################################################################################\n"
    "                                       EvoEF %s\n"
    "  A framework for macromolecular modeling, e.g.,protein design, protein side-chain packing,\n"
    "protein structure energy minimization, add and optimize hydrogen bonds, build mutant model,\n"
    "calculate protein folding stability, calculate protein-protein binding free energy, etc\n"
    "\n\n"
    "  Copyright (c) 2020 Xiaoqiang Huang\n"
    "  tommyhuangthu@foxmail.com; xiaoqiah@umich.edu\n"
    "  Dept. of Computational Medicine & Bioinformatics\n"
    "  Medical School\n"
    "  University of Michigan\n"
    "############################################################################################\n",
    PROGRAM_VERSION);
  return Success;
}


///////////////////////////////////////////////////////////////////////////////////////////
//MAIN functions
//////////////////////////////////////////////////////////////////////////////////////////
int EVOEF_ComputeChainStability(Structure *pStructure, int chainIndex, double *energyTerms){
  Chain *pChainI = StructureGetChain(pStructure, chainIndex);
  // interaction between residues within one chain
  for(int ir = 0; ir < ChainGetResidueCount(pChainI); ir++){
    Residue *pResIR = ChainGetResidue(pChainI,ir);
    EVOEF_AminoAcidReferenceEnergy(pResIR->name, energyTerms);
    EVOEF_EnergyResidueIntraEnergy(pResIR,energyTerms);
    for(int is = ir+1; is < ChainGetResidueCount(pChainI); is++){
      Residue *pResIS = ChainGetResidue(pChainI,is);
      if(ResidueGetPosInChain(pResIR)+1==ResidueGetPosInChain(pResIS)) EVOEF_EnergyResidueAndNextResidue(pResIR,pResIS,energyTerms);
      else EVOEF_EnergyResidueAndOtherResidueSameChain(pResIR,pResIS,energyTerms);
      //EVOEF_EnergyResidueAndResidueSameChain(pResIR,pResIS,energyTerms);
    }
  }

  EnergyTermWeighting(energyTerms);
  printf("Chain %s energy details:\n", ChainGetName(pChainI));
  printf("reference_ALA         =            %12.6f\n", energyTerms[ 1]);
  printf("reference_CYS         =            %12.6f\n", energyTerms[ 2]);
  printf("reference_ASP         =            %12.6f\n", energyTerms[ 3]);
  printf("reference_GLU         =            %12.6f\n", energyTerms[ 4]);
  printf("reference_PHE         =            %12.6f\n", energyTerms[ 5]);
  printf("reference_GLY         =            %12.6f\n", energyTerms[ 6]);
  printf("reference_HIS         =            %12.6f\n", energyTerms[ 7]);
  printf("reference_ILE         =            %12.6f\n", energyTerms[ 8]);
  printf("reference_LYS         =            %12.6f\n", energyTerms[ 9]);
  printf("reference_LEU         =            %12.6f\n", energyTerms[10]);
  printf("reference_MET         =            %12.6f\n", energyTerms[11]);
  printf("reference_ASN         =            %12.6f\n", energyTerms[12]);
  printf("reference_PRO         =            %12.6f\n", energyTerms[13]);
  printf("reference_GLN         =            %12.6f\n", energyTerms[14]);
  printf("reference_ARG         =            %12.6f\n", energyTerms[15]);
  printf("reference_SER         =            %12.6f\n", energyTerms[16]);
  printf("reference_THR         =            %12.6f\n", energyTerms[17]);
  printf("reference_VAL         =            %12.6f\n", energyTerms[18]);
  printf("reference_TRP         =            %12.6f\n", energyTerms[19]);
  printf("reference_TYR         =            %12.6f\n", energyTerms[20]);
  printf("intraR_vdwatt         =            %12.6f\n", energyTerms[21]);
  printf("intraR_vdwrep         =            %12.6f\n", energyTerms[22]);
  printf("intraR_electr         =            %12.6f\n", energyTerms[23]);
  printf("intraR_deslvP         =            %12.6f\n", energyTerms[24]);
  printf("intraR_deslvH         =            %12.6f\n", energyTerms[25]);
  printf("intraR_hbscbb_dis     =            %12.6f\n", energyTerms[26]);
  printf("intraR_hbscbb_the     =            %12.6f\n", energyTerms[27]);
  printf("intraR_hbscbb_phi     =            %12.6f\n", energyTerms[28]);
  printf("interS_vdwatt         =            %12.6f\n", energyTerms[31]);
  printf("interS_vdwrep         =            %12.6f\n", energyTerms[32]);
  printf("interS_electr         =            %12.6f\n", energyTerms[33]);
  printf("interS_deslvP         =            %12.6f\n", energyTerms[34]);
  printf("interS_deslvH         =            %12.6f\n", energyTerms[35]);
  printf("interS_hbbbbb_dis     =            %12.6f\n", energyTerms[41]);
  printf("interS_hbbbbb_the     =            %12.6f\n", energyTerms[42]);
  printf("interS_hbbbbb_phi     =            %12.6f\n", energyTerms[43]);
  printf("interS_hbscbb_dis     =            %12.6f\n", energyTerms[44]);
  printf("interS_hbscbb_the     =            %12.6f\n", energyTerms[45]);
  printf("interS_hbscbb_phi     =            %12.6f\n", energyTerms[46]);
  printf("interS_hbscsc_dis     =            %12.6f\n", energyTerms[47]);
  printf("interS_hbscsc_the     =            %12.6f\n", energyTerms[48]);
  printf("interS_hbscsc_phi     =            %12.6f\n", energyTerms[49]);
  printf("----------------------------------------------------\n");
  printf("Total                 =            %12.6f\n\n", energyTerms[0]);

  return Success;
}




int EVOEF_ComputeStability(Structure *pStructure,AAppTable* pAAppTable,RamaTable* pRama,double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]){
  int aas[20]={0}; //ACDEFGHIKLMNPQRSTVWY, only for regular amino acid
  //StructureGetAminoAcidComposition(pStructure, aas);
  // if the structure is composed of several chains, the residue position could be different in the whole structure from that in the separate chain
  //StructureComputeResiduePosition(pStructure);
  for(int i = 0; i < StructureGetChainCount(pStructure); i++){
    Chain *pChainI = StructureGetChain(pStructure,i);
    for(int ir = 0; ir < ChainGetResidueCount(pChainI); ir++){
      Residue *pResIR = ChainGetResidue(pChainI,ir);
      double refer=0.0;
      if(ChainGetType(pChainI)==Type_Chain_Protein){
        EVOEF_AminoAcidReferenceEnergy(pResIR->name,energyTerms);
        EVOEF_EnergyResidueIntraEnergy(pResIR,energyTerms);
        AminoAcidPropensityAndRamachandranEnergy(pResIR,pAAppTable,pRama);
        energyTerms[91]+=pResIR->aapp;
        energyTerms[92]+=pResIR->ramachandran;
      }
      for(int is = ir+1; is < ChainGetResidueCount(pChainI); is++){
        Residue *pResIS = ChainGetResidue(pChainI,is);
        if(ResidueGetPosInChain(pResIR)+1==ResidueGetPosInChain(pResIS)) EVOEF_EnergyResidueAndNextResidue(pResIR,pResIS,energyTerms);
        else EVOEF_EnergyResidueAndOtherResidueSameChain(pResIR,pResIS,energyTerms);
        //EVOEF_EnergyResidueAndResidueSameChain(pResIR,pResIS,energyTerms);
      }
      for(int k = i+1; k < StructureGetChainCount(pStructure); k++){
        Chain *pChainK = StructureGetChain(pStructure,k);
        for(int ks = 0; ks < ChainGetResidueCount(pChainK); ks++){
          Residue *pResKS = ChainGetResidue(pChainK,ks);
          if(ChainGetType(pChainI)==Type_Chain_SmallMol){
            EVOEF_EnergyResidueAndLigandResidue(pResKS,pResIR,energyTerms);
          }
          else if(ChainGetType(pChainK)==Type_Chain_SmallMol){
            EVOEF_EnergyResidueAndLigandResidue(pResIR,pResKS,energyTerms);
          }
          else{
            EVOEF_EnergyResidueAndOtherResidueDifferentChain(pResIR,pResKS,energyTerms);
          }
        }
      }
    }
  }

  EnergyTermWeighting(energyTerms);
  printf("\nStructure energy details:\n");
  printf("reference_ALA         =            %12.6f\n", energyTerms[ 1]);
  printf("reference_CYS         =            %12.6f\n", energyTerms[ 2]);
  printf("reference_ASP         =            %12.6f\n", energyTerms[ 3]);
  printf("reference_GLU         =            %12.6f\n", energyTerms[ 4]);
  printf("reference_PHE         =            %12.6f\n", energyTerms[ 5]);
  printf("reference_GLY         =            %12.6f\n", energyTerms[ 6]);
  printf("reference_HIS         =            %12.6f\n", energyTerms[ 7]);
  printf("reference_ILE         =            %12.6f\n", energyTerms[ 8]);
  printf("reference_LYS         =            %12.6f\n", energyTerms[ 9]);
  printf("reference_LEU         =            %12.6f\n", energyTerms[10]);
  printf("reference_MET         =            %12.6f\n", energyTerms[11]);
  printf("reference_ASN         =            %12.6f\n", energyTerms[12]);
  printf("reference_PRO         =            %12.6f\n", energyTerms[13]);
  printf("reference_GLN         =            %12.6f\n", energyTerms[14]);
  printf("reference_ARG         =            %12.6f\n", energyTerms[15]);
  printf("reference_SER         =            %12.6f\n", energyTerms[16]);
  printf("reference_THR         =            %12.6f\n", energyTerms[17]);
  printf("reference_VAL         =            %12.6f\n", energyTerms[18]);
  printf("reference_TRP         =            %12.6f\n", energyTerms[19]);
  printf("reference_TYR         =            %12.6f\n", energyTerms[20]);

  printf("intraR_vdwatt         =            %12.6f\n", energyTerms[21]);
  printf("intraR_vdwrep         =            %12.6f\n", energyTerms[22]);
  printf("intraR_electr         =            %12.6f\n", energyTerms[23]);
  printf("intraR_deslvP         =            %12.6f\n", energyTerms[24]);
  printf("intraR_deslvH         =            %12.6f\n", energyTerms[25]);
  printf("intraR_hbscbb_dis     =            %12.6f\n", energyTerms[26]);
  printf("intraR_hbscbb_the     =            %12.6f\n", energyTerms[27]);
  printf("intraR_hbscbb_phi     =            %12.6f\n", energyTerms[28]);
  printf("aapropensity          =            %12.6f\n", energyTerms[91]);
  printf("ramachandran          =            %12.6f\n", energyTerms[92]);
  printf("dunbrack              =            %12.6f\n", energyTerms[93]);

  printf("interS_vdwatt         =            %12.6f\n", energyTerms[31]);
  printf("interS_vdwrep         =            %12.6f\n", energyTerms[32]);
  printf("interS_electr         =            %12.6f\n", energyTerms[33]);
  printf("interS_deslvP         =            %12.6f\n", energyTerms[34]);
  printf("interS_deslvH         =            %12.6f\n", energyTerms[35]);
  printf("interS_ssbond         =            %12.6f\n", energyTerms[36]);
  printf("interS_hbbbbb_dis     =            %12.6f\n", energyTerms[41]);
  printf("interS_hbbbbb_the     =            %12.6f\n", energyTerms[42]);
  printf("interS_hbbbbb_phi     =            %12.6f\n", energyTerms[43]);
  printf("interS_hbscbb_dis     =            %12.6f\n", energyTerms[44]);
  printf("interS_hbscbb_the     =            %12.6f\n", energyTerms[45]);
  printf("interS_hbscbb_phi     =            %12.6f\n", energyTerms[46]);
  printf("interS_hbscsc_dis     =            %12.6f\n", energyTerms[47]);
  printf("interS_hbscsc_the     =            %12.6f\n", energyTerms[48]);
  printf("interS_hbscsc_phi     =            %12.6f\n", energyTerms[49]);

  printf("interD_vdwatt         =            %12.6f\n", energyTerms[51]);
  printf("interD_vdwrep         =            %12.6f\n", energyTerms[52]);
  printf("interD_electr         =            %12.6f\n", energyTerms[53]);
  printf("interD_deslvP         =            %12.6f\n", energyTerms[54]);
  printf("interD_deslvH         =            %12.6f\n", energyTerms[55]);
  printf("interD_ssbond         =            %12.6f\n", energyTerms[56]);
  printf("interD_hbbbbb_dis     =            %12.6f\n", energyTerms[61]);
  printf("interD_hbbbbb_the     =            %12.6f\n", energyTerms[62]);
  printf("interD_hbbbbb_phi     =            %12.6f\n", energyTerms[63]);
  printf("interD_hbscbb_dis     =            %12.6f\n", energyTerms[64]);
  printf("interD_hbscbb_the     =            %12.6f\n", energyTerms[65]);
  printf("interD_hbscbb_phi     =            %12.6f\n", energyTerms[66]);
  printf("interD_hbscsc_dis     =            %12.6f\n", energyTerms[67]);
  printf("interD_hbscsc_the     =            %12.6f\n", energyTerms[68]);
  printf("interD_hbscsc_phi     =            %12.6f\n", energyTerms[69]);

  printf("prolig_vdwatt         =            %12.6f\n", energyTerms[71]);
  printf("prolig_vdwrep         =            %12.6f\n", energyTerms[72]);
  printf("prolig_electr         =            %12.6f\n", energyTerms[73]);
  printf("prolig_deslvP         =            %12.6f\n", energyTerms[74]);
  printf("prolig_deslvH         =            %12.6f\n", energyTerms[75]);
  printf("prolig_hbscbb_dis     =            %12.6f\n", energyTerms[81]);
  printf("prolig_hbscbb_the     =            %12.6f\n", energyTerms[82]);
  printf("prolig_hbscbb_phi     =            %12.6f\n", energyTerms[83]);
  printf("prolig_hbscsc_dis     =            %12.6f\n", energyTerms[84]);
  printf("prolig_hbscsc_the     =            %12.6f\n", energyTerms[85]);
  printf("prolig_hbscsc_phi     =            %12.6f\n", energyTerms[86]);
  printf("----------------------------------------------------\n");
  printf("Total                 =            %12.6f\n\n", energyTerms[0]);
  return Success;
}


int EVOEF_ComputeStabilityWithBBdepRotLib(Structure *pStructure,AAppTable* pAAppTable,RamaTable* pRama,BBdepRotamerLib* pRotLib,double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]){
  int aas[20]={0}; //ACDEFGHIKLMNPQRSTVWY, only for regular amino acid
  //StructureGetAminoAcidComposition(pStructure, aas);
  // if the structure is composed of several chains, the residue position could be different in the whole structure from that in the separate chain
  //StructureComputeResiduePosition(pStructure);
  for(int i = 0; i < StructureGetChainCount(pStructure); i++){
    Chain *pChainI = StructureGetChain(pStructure,i);
    for(int ir = 0; ir < ChainGetResidueCount(pChainI); ir++){
      Residue *pResIR = ChainGetResidue(pChainI,ir);
      double refer=0.0;
      if(ChainGetType(pChainI)==Type_Chain_Protein){
        EVOEF_AminoAcidReferenceEnergy(pResIR->name, energyTerms);
        EVOEF_EnergyResidueIntraEnergy(pResIR,energyTerms);
        AminoAcidPropensityAndRamachandranEnergy(pResIR,pAAppTable,pRama);
        energyTerms[91]+=pResIR->aapp;
        energyTerms[92]+=pResIR->ramachandran;
        AminoAcidDunbrackEnergy(pResIR,pRotLib);
        energyTerms[93]+=pResIR->dunbrack;
      }
      for(int is = ir+1; is < ChainGetResidueCount(pChainI); is++){
        Residue *pResIS = ChainGetResidue(pChainI,is);
        if(ResidueGetPosInChain(pResIR)+1==ResidueGetPosInChain(pResIS)) EVOEF_EnergyResidueAndNextResidue(pResIR,pResIS,energyTerms);
        else EVOEF_EnergyResidueAndOtherResidueSameChain(pResIR,pResIS,energyTerms);
        //EVOEF_EnergyResidueAndResidueSameChain(pResIR,pResIS,energyTerms);
      }
      for(int k = i+1; k < StructureGetChainCount(pStructure); k++){
        Chain *pChainK = StructureGetChain(pStructure,k);
        for(int ks = 0; ks < ChainGetResidueCount(pChainK); ks++){
          Residue *pResKS = ChainGetResidue(pChainK,ks);
          if(ChainGetType(pChainI)==Type_Chain_SmallMol){
            EVOEF_EnergyResidueAndLigandResidue(pResKS,pResIR,energyTerms);
          }
          else if(ChainGetType(pChainK)==Type_Chain_SmallMol){
            EVOEF_EnergyResidueAndLigandResidue(pResIR,pResKS,energyTerms);
          }
          else{
            EVOEF_EnergyResidueAndOtherResidueDifferentChain(pResIR,pResKS,energyTerms);
          }
        }
      }
    }
  }

  EnergyTermWeighting(energyTerms);
  printf("\nStructure energy details:\n");
  printf("reference_ALA         =            %12.6f\n", energyTerms[ 1]);
  printf("reference_CYS         =            %12.6f\n", energyTerms[ 2]);
  printf("reference_ASP         =            %12.6f\n", energyTerms[ 3]);
  printf("reference_GLU         =            %12.6f\n", energyTerms[ 4]);
  printf("reference_PHE         =            %12.6f\n", energyTerms[ 5]);
  printf("reference_GLY         =            %12.6f\n", energyTerms[ 6]);
  printf("reference_HIS         =            %12.6f\n", energyTerms[ 7]);
  printf("reference_ILE         =            %12.6f\n", energyTerms[ 8]);
  printf("reference_LYS         =            %12.6f\n", energyTerms[ 9]);
  printf("reference_LEU         =            %12.6f\n", energyTerms[10]);
  printf("reference_MET         =            %12.6f\n", energyTerms[11]);
  printf("reference_ASN         =            %12.6f\n", energyTerms[12]);
  printf("reference_PRO         =            %12.6f\n", energyTerms[13]);
  printf("reference_GLN         =            %12.6f\n", energyTerms[14]);
  printf("reference_ARG         =            %12.6f\n", energyTerms[15]);
  printf("reference_SER         =            %12.6f\n", energyTerms[16]);
  printf("reference_THR         =            %12.6f\n", energyTerms[17]);
  printf("reference_VAL         =            %12.6f\n", energyTerms[18]);
  printf("reference_TRP         =            %12.6f\n", energyTerms[19]);
  printf("reference_TYR         =            %12.6f\n", energyTerms[20]);

  printf("intraR_vdwatt         =            %12.6f\n", energyTerms[21]);
  printf("intraR_vdwrep         =            %12.6f\n", energyTerms[22]);
  printf("intraR_electr         =            %12.6f\n", energyTerms[23]);
  printf("intraR_deslvP         =            %12.6f\n", energyTerms[24]);
  printf("intraR_deslvH         =            %12.6f\n", energyTerms[25]);
  printf("intraR_hbscbb_dis     =            %12.6f\n", energyTerms[26]);
  printf("intraR_hbscbb_the     =            %12.6f\n", energyTerms[27]);
  printf("intraR_hbscbb_phi     =            %12.6f\n", energyTerms[28]);
  printf("aapropensity          =            %12.6f\n", energyTerms[91]);
  printf("ramachandran          =            %12.6f\n", energyTerms[92]);
  printf("dunbrack              =            %12.6f\n", energyTerms[93]);

  printf("interS_vdwatt         =            %12.6f\n", energyTerms[31]);
  printf("interS_vdwrep         =            %12.6f\n", energyTerms[32]);
  printf("interS_electr         =            %12.6f\n", energyTerms[33]);
  printf("interS_deslvP         =            %12.6f\n", energyTerms[34]);
  printf("interS_deslvH         =            %12.6f\n", energyTerms[35]);
  printf("interS_ssbond         =            %12.6f\n", energyTerms[36]);
  printf("interS_hbbbbb_dis     =            %12.6f\n", energyTerms[41]);
  printf("interS_hbbbbb_the     =            %12.6f\n", energyTerms[42]);
  printf("interS_hbbbbb_phi     =            %12.6f\n", energyTerms[43]);
  printf("interS_hbscbb_dis     =            %12.6f\n", energyTerms[44]);
  printf("interS_hbscbb_the     =            %12.6f\n", energyTerms[45]);
  printf("interS_hbscbb_phi     =            %12.6f\n", energyTerms[46]);
  printf("interS_hbscsc_dis     =            %12.6f\n", energyTerms[47]);
  printf("interS_hbscsc_the     =            %12.6f\n", energyTerms[48]);
  printf("interS_hbscsc_phi     =            %12.6f\n", energyTerms[49]);

  printf("interD_vdwatt         =            %12.6f\n", energyTerms[51]);
  printf("interD_vdwrep         =            %12.6f\n", energyTerms[52]);
  printf("interD_electr         =            %12.6f\n", energyTerms[53]);
  printf("interD_deslvP         =            %12.6f\n", energyTerms[54]);
  printf("interD_deslvH         =            %12.6f\n", energyTerms[55]);
  printf("interD_ssbond         =            %12.6f\n", energyTerms[56]);
  printf("interD_hbbbbb_dis     =            %12.6f\n", energyTerms[61]);
  printf("interD_hbbbbb_the     =            %12.6f\n", energyTerms[62]);
  printf("interD_hbbbbb_phi     =            %12.6f\n", energyTerms[63]);
  printf("interD_hbscbb_dis     =            %12.6f\n", energyTerms[64]);
  printf("interD_hbscbb_the     =            %12.6f\n", energyTerms[65]);
  printf("interD_hbscbb_phi     =            %12.6f\n", energyTerms[66]);
  printf("interD_hbscsc_dis     =            %12.6f\n", energyTerms[67]);
  printf("interD_hbscsc_the     =            %12.6f\n", energyTerms[68]);
  printf("interD_hbscsc_phi     =            %12.6f\n", energyTerms[69]);

  printf("prolig_vdwatt         =            %12.6f\n", energyTerms[71]);
  printf("prolig_vdwrep         =            %12.6f\n", energyTerms[72]);
  printf("prolig_electr         =            %12.6f\n", energyTerms[73]);
  printf("prolig_deslvP         =            %12.6f\n", energyTerms[74]);
  printf("prolig_deslvH         =            %12.6f\n", energyTerms[75]);
  printf("prolig_hbscbb_dis     =            %12.6f\n", energyTerms[81]);
  printf("prolig_hbscbb_the     =            %12.6f\n", energyTerms[82]);
  printf("prolig_hbscbb_phi     =            %12.6f\n", energyTerms[83]);
  printf("prolig_hbscsc_dis     =            %12.6f\n", energyTerms[84]);
  printf("prolig_hbscsc_the     =            %12.6f\n", energyTerms[85]);
  printf("prolig_hbscsc_phi     =            %12.6f\n", energyTerms[86]);
  printf("----------------------------------------------------\n");
  printf("Total                 =            %12.6f\n\n", energyTerms[0]);
  return Success;
}

int EVOEF_ComputeBinding(Structure *pStructure){
  if(StructureGetChainCount(pStructure)>2){
    printf("Your structure has more than two protein chains, and you should specify how to split chains "
      "before computing the binding energy\n");
    printf("Otherwise, EvoEF just output the interactions between any chain pair (DEFAULT)\n");
  }
  else if(StructureGetChainCount(pStructure)<=1){
    printf("Your structure has less than or equal to one chain, binding energy cannot be calculated\n");
    return Warning;
  }

  for(int i=0; i<StructureGetChainCount(pStructure);i++){
    Chain* pChainI=StructureGetChain(pStructure,i);
    for(int k=i+1; k<StructureGetChainCount(pStructure);k++){
      Chain* pChainK=StructureGetChain(pStructure,k);
      double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
      for(int j=0;j<ChainGetResidueCount(pChainI);j++){
        Residue* pResIJ=ChainGetResidue(pChainI,j);
        for(int s=0;s<ChainGetResidueCount(pChainK);s++){
          Residue* pResKS=ChainGetResidue(pChainK,s);
          if(ChainGetType(pChainI)==Type_Chain_SmallMol){
            EVOEF_EnergyResidueAndLigandResidue(pResKS,pResIJ,energyTerms);
          }
          else if(ChainGetType(pChainK)==Type_Chain_SmallMol){
            EVOEF_EnergyResidueAndLigandResidue(pResIJ,pResKS,energyTerms);
          }
          else{
            EVOEF_EnergyResidueAndOtherResidueDifferentChain(pResIJ,pResKS,energyTerms);
          }
        }
      }
      EnergyTermWeighting(energyTerms);
      // energy terms are weighted during the calculation, don't weight them for the difference
      printf("Binding energy details between chain(s) %s and chain(s) %s:\n",
        ChainGetName(pChainI),ChainGetName(pChainK),ChainGetName(pChainI),ChainGetName(pChainK));
      printf("reference_ALA         =            %12.6f\n", energyTerms[ 1]);
      printf("reference_CYS         =            %12.6f\n", energyTerms[ 2]);
      printf("reference_ASP         =            %12.6f\n", energyTerms[ 3]);
      printf("reference_GLU         =            %12.6f\n", energyTerms[ 4]);
      printf("reference_PHE         =            %12.6f\n", energyTerms[ 5]);
      printf("reference_GLY         =            %12.6f\n", energyTerms[ 6]);
      printf("reference_HIS         =            %12.6f\n", energyTerms[ 7]);
      printf("reference_ILE         =            %12.6f\n", energyTerms[ 8]);
      printf("reference_LYS         =            %12.6f\n", energyTerms[ 9]);
      printf("reference_LEU         =            %12.6f\n", energyTerms[10]);
      printf("reference_MET         =            %12.6f\n", energyTerms[11]);
      printf("reference_ASN         =            %12.6f\n", energyTerms[12]);
      printf("reference_PRO         =            %12.6f\n", energyTerms[13]);
      printf("reference_GLN         =            %12.6f\n", energyTerms[14]);
      printf("reference_ARG         =            %12.6f\n", energyTerms[15]);
      printf("reference_SER         =            %12.6f\n", energyTerms[16]);
      printf("reference_THR         =            %12.6f\n", energyTerms[17]);
      printf("reference_VAL         =            %12.6f\n", energyTerms[18]);
      printf("reference_TRP         =            %12.6f\n", energyTerms[19]);
      printf("reference_TYR         =            %12.6f\n", energyTerms[20]);

      printf("intraR_vdwatt         =            %12.6f\n", energyTerms[21]);
      printf("intraR_vdwrep         =            %12.6f\n", energyTerms[22]);
      printf("intraR_electr         =            %12.6f\n", energyTerms[23]);
      printf("intraR_deslvP         =            %12.6f\n", energyTerms[24]);
      printf("intraR_deslvH         =            %12.6f\n", energyTerms[25]);
      printf("intraR_hbscbb_dis     =            %12.6f\n", energyTerms[26]);
      printf("intraR_hbscbb_the     =            %12.6f\n", energyTerms[27]);
      printf("intraR_hbscbb_phi     =            %12.6f\n", energyTerms[28]);

      printf("interS_vdwatt         =            %12.6f\n", energyTerms[31]);
      printf("interS_vdwrep         =            %12.6f\n", energyTerms[32]);
      printf("interS_electr         =            %12.6f\n", energyTerms[33]);
      printf("interS_deslvP         =            %12.6f\n", energyTerms[34]);
      printf("interS_deslvH         =            %12.6f\n", energyTerms[35]);
      printf("interS_ssbond         =            %12.6f\n", energyTerms[36]);
      printf("interS_hbbbbb_dis     =            %12.6f\n", energyTerms[41]);
      printf("interS_hbbbbb_the     =            %12.6f\n", energyTerms[42]);
      printf("interS_hbbbbb_phi     =            %12.6f\n", energyTerms[43]);
      printf("interS_hbscbb_dis     =            %12.6f\n", energyTerms[44]);
      printf("interS_hbscbb_the     =            %12.6f\n", energyTerms[45]);
      printf("interS_hbscbb_phi     =            %12.6f\n", energyTerms[46]);
      printf("interS_hbscsc_dis     =            %12.6f\n", energyTerms[47]);
      printf("interS_hbscsc_the     =            %12.6f\n", energyTerms[48]);
      printf("interS_hbscsc_phi     =            %12.6f\n", energyTerms[49]);

      printf("interD_vdwatt         =            %12.6f\n", energyTerms[51]);
      printf("interD_vdwrep         =            %12.6f\n", energyTerms[52]);
      printf("interD_electr         =            %12.6f\n", energyTerms[53]);
      printf("interD_deslvP         =            %12.6f\n", energyTerms[54]);
      printf("interD_deslvH         =            %12.6f\n", energyTerms[55]);
      printf("interD_ssbond         =            %12.6f\n", energyTerms[56]);
      printf("interD_hbbbbb_dis     =            %12.6f\n", energyTerms[61]);
      printf("interD_hbbbbb_the     =            %12.6f\n", energyTerms[62]);
      printf("interD_hbbbbb_phi     =            %12.6f\n", energyTerms[63]);
      printf("interD_hbscbb_dis     =            %12.6f\n", energyTerms[64]);
      printf("interD_hbscbb_the     =            %12.6f\n", energyTerms[65]);
      printf("interD_hbscbb_phi     =            %12.6f\n", energyTerms[66]);
      printf("interD_hbscsc_dis     =            %12.6f\n", energyTerms[67]);
      printf("interD_hbscsc_the     =            %12.6f\n", energyTerms[68]);
      printf("interD_hbscsc_phi     =            %12.6f\n", energyTerms[69]);

      printf("prolig_vdwatt         =            %12.6f\n", energyTerms[71]);
      printf("prolig_vdwrep         =            %12.6f\n", energyTerms[72]);
      printf("prolig_electr         =            %12.6f\n", energyTerms[73]);
      printf("prolig_deslvP         =            %12.6f\n", energyTerms[74]);
      printf("prolig_deslvH         =            %12.6f\n", energyTerms[75]);
      printf("prolig_hbscbb_dis     =            %12.6f\n", energyTerms[81]);
      printf("prolig_hbscbb_the     =            %12.6f\n", energyTerms[82]);
      printf("prolig_hbscbb_phi     =            %12.6f\n", energyTerms[83]);
      printf("prolig_hbscsc_dis     =            %12.6f\n", energyTerms[84]);
      printf("prolig_hbscsc_the     =            %12.6f\n", energyTerms[85]);
      printf("prolig_hbscsc_phi     =            %12.6f\n", energyTerms[86]);
      printf("----------------------------------------------------\n");
      printf("Total                 =            %12.6f\n", energyTerms[0]);
    }
  }

  return Success;
}


int EVOEF_ComputeBindingWithSplittingNew(Structure *pStructure,char split1[], char split2[]){
  double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
  for(int i=0;i<StructureGetChainCount(pStructure);i++){
    Chain* pChainI=StructureGetChain(pStructure,i);
    for(int k=i+1;k<StructureGetChainCount(pStructure);k++){
      Chain* pChainK=StructureGetChain(pStructure,k);
      if((strstr(split1,ChainGetName(pChainI))!=NULL && strstr(split2,ChainGetName(pChainK))!=NULL)||
        ((strstr(split2,ChainGetName(pChainI))!=NULL && strstr(split1,ChainGetName(pChainK))!=NULL))){
          for(int j=0;j<ChainGetResidueCount(pChainI);j++){
            Residue* pResIJ=ChainGetResidue(pChainI,j);
            for(int s=0;s<ChainGetResidueCount(pChainK);s++){
              Residue* pResKS=ChainGetResidue(pChainK,s);
              if(ChainGetType(pChainI)==Type_Chain_SmallMol){
                EVOEF_EnergyResidueAndLigandResidue(pResKS,pResIJ,energyTerms);
              }
              else if(ChainGetType(pChainK)==Type_Chain_SmallMol){
                EVOEF_EnergyResidueAndLigandResidue(pResIJ,pResKS,energyTerms);
              }
              else{
                EVOEF_EnergyResidueAndOtherResidueDifferentChain(pResIJ,pResKS,energyTerms);
              }
            }
          }
      }
    }
  }
  EnergyTermWeighting(energyTerms);
  // energy terms are weighted during the calculation, don't weight them for the difference
  printf("Binding energy details between chain(s) %s and chain(s) %s (DG_bind = DG(stability,complex) - DG(stability,%s) - DG(stability,%s):\n",split1,split2,split1,split2);
  printf("reference_ALA         =            %12.6f\n", energyTerms[ 1]);
  printf("reference_CYS         =            %12.6f\n", energyTerms[ 2]);
  printf("reference_ASP         =            %12.6f\n", energyTerms[ 3]);
  printf("reference_GLU         =            %12.6f\n", energyTerms[ 4]);
  printf("reference_PHE         =            %12.6f\n", energyTerms[ 5]);
  printf("reference_GLY         =            %12.6f\n", energyTerms[ 6]);
  printf("reference_HIS         =            %12.6f\n", energyTerms[ 7]);
  printf("reference_ILE         =            %12.6f\n", energyTerms[ 8]);
  printf("reference_LYS         =            %12.6f\n", energyTerms[ 9]);
  printf("reference_LEU         =            %12.6f\n", energyTerms[10]);
  printf("reference_MET         =            %12.6f\n", energyTerms[11]);
  printf("reference_ASN         =            %12.6f\n", energyTerms[12]);
  printf("reference_PRO         =            %12.6f\n", energyTerms[13]);
  printf("reference_GLN         =            %12.6f\n", energyTerms[14]);
  printf("reference_ARG         =            %12.6f\n", energyTerms[15]);
  printf("reference_SER         =            %12.6f\n", energyTerms[16]);
  printf("reference_THR         =            %12.6f\n", energyTerms[17]);
  printf("reference_VAL         =            %12.6f\n", energyTerms[18]);
  printf("reference_TRP         =            %12.6f\n", energyTerms[19]);
  printf("reference_TYR         =            %12.6f\n", energyTerms[20]);

  printf("intraR_vdwatt         =            %12.6f\n", energyTerms[21]);
  printf("intraR_vdwrep         =            %12.6f\n", energyTerms[22]);
  printf("intraR_electr         =            %12.6f\n", energyTerms[23]);
  printf("intraR_deslvP         =            %12.6f\n", energyTerms[24]);
  printf("intraR_deslvH         =            %12.6f\n", energyTerms[25]);
  printf("intraR_hbscbb_dis     =            %12.6f\n", energyTerms[26]);
  printf("intraR_hbscbb_the     =            %12.6f\n", energyTerms[27]);
  printf("intraR_hbscbb_phi     =            %12.6f\n", energyTerms[28]);

  printf("interS_vdwatt         =            %12.6f\n", energyTerms[31]);
  printf("interS_vdwrep         =            %12.6f\n", energyTerms[32]);
  printf("interS_electr         =            %12.6f\n", energyTerms[33]);
  printf("interS_deslvP         =            %12.6f\n", energyTerms[34]);
  printf("interS_deslvH         =            %12.6f\n", energyTerms[35]);
  printf("interS_ssbond         =            %12.6f\n", energyTerms[36]);
  printf("interS_hbbbbb_dis     =            %12.6f\n", energyTerms[41]);
  printf("interS_hbbbbb_the     =            %12.6f\n", energyTerms[42]);
  printf("interS_hbbbbb_phi     =            %12.6f\n", energyTerms[43]);
  printf("interS_hbscbb_dis     =            %12.6f\n", energyTerms[44]);
  printf("interS_hbscbb_the     =            %12.6f\n", energyTerms[45]);
  printf("interS_hbscbb_phi     =            %12.6f\n", energyTerms[46]);
  printf("interS_hbscsc_dis     =            %12.6f\n", energyTerms[47]);
  printf("interS_hbscsc_the     =            %12.6f\n", energyTerms[48]);
  printf("interS_hbscsc_phi     =            %12.6f\n", energyTerms[49]);

  printf("interD_vdwatt         =            %12.6f\n", energyTerms[51]);
  printf("interD_vdwrep         =            %12.6f\n", energyTerms[52]);
  printf("interD_electr         =            %12.6f\n", energyTerms[53]);
  printf("interD_deslvP         =            %12.6f\n", energyTerms[54]);
  printf("interD_deslvH         =            %12.6f\n", energyTerms[55]);
  printf("interD_ssbond         =            %12.6f\n", energyTerms[56]);
  printf("interD_hbbbbb_dis     =            %12.6f\n", energyTerms[61]);
  printf("interD_hbbbbb_the     =            %12.6f\n", energyTerms[62]);
  printf("interD_hbbbbb_phi     =            %12.6f\n", energyTerms[63]);
  printf("interD_hbscbb_dis     =            %12.6f\n", energyTerms[64]);
  printf("interD_hbscbb_the     =            %12.6f\n", energyTerms[65]);
  printf("interD_hbscbb_phi     =            %12.6f\n", energyTerms[66]);
  printf("interD_hbscsc_dis     =            %12.6f\n", energyTerms[67]);
  printf("interD_hbscsc_the     =            %12.6f\n", energyTerms[68]);
  printf("interD_hbscsc_phi     =            %12.6f\n", energyTerms[69]);

  printf("prolig_vdwatt         =            %12.6f\n", energyTerms[71]);
  printf("prolig_vdwrep         =            %12.6f\n", energyTerms[72]);
  printf("prolig_electr         =            %12.6f\n", energyTerms[73]);
  printf("prolig_deslvP         =            %12.6f\n", energyTerms[74]);
  printf("prolig_deslvH         =            %12.6f\n", energyTerms[75]);
  printf("prolig_hbscbb_dis     =            %12.6f\n", energyTerms[81]);
  printf("prolig_hbscbb_the     =            %12.6f\n", energyTerms[82]);
  printf("prolig_hbscbb_phi     =            %12.6f\n", energyTerms[83]);
  printf("prolig_hbscsc_dis     =            %12.6f\n", energyTerms[84]);
  printf("prolig_hbscsc_the     =            %12.6f\n", energyTerms[85]);
  printf("prolig_hbscsc_phi     =            %12.6f\n", energyTerms[86]);
  printf("----------------------------------------------------\n");
  printf("Total                 =            %12.6f\n", energyTerms[0]);
  return Success;
}


//this function is used to build the structure model of mutations
int EVOEF_BuildMutant(Structure* pStructure, char* mutantfile, BBindRotamerLib* rotlib, AtomParamsSet* atomParams,ResiTopoSet* resiTopos, char* pdbid){
  FileReader fr;
  FileReaderCreate(&fr, mutantfile);
  int mutantcount = FileReaderGetLineCount(&fr);
  if(mutantcount<=0){
    printf("There is no mutant found in the mutant file\n");
    return DataNotExistError;
  }

  StringArray* mutants = (StringArray*)malloc(sizeof(StringArray)*mutantcount);
  char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  int mutantIndex=0;
  while(!FAILED(FileReaderGetNextLine(&fr, line))){
    StringArrayCreate(&mutants[mutantIndex]);
    StringArraySplitString(&mutants[mutantIndex], line, ',');
    char lastMutant[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    int lastmutindex = StringArrayGetCount(&mutants[mutantIndex])-1;
    strcpy(lastMutant, StringArrayGet(&mutants[mutantIndex], lastmutindex));
    //deal with the last char of the last single mutant
    if((!isdigit(lastMutant[strlen(lastMutant)-1])) && !isalpha(lastMutant[strlen(lastMutant)-1])){
      lastMutant[strlen(lastMutant)-1] = '\0';
    }
    StringArraySet(&mutants[mutantIndex], lastmutindex, lastMutant);
    mutantIndex++;
  }
  FileReaderDestroy(&fr);

  for(int mutantIndex = 0; mutantIndex < mutantcount; mutantIndex++){
    Structure tempStruct;
    StructureCreate(&tempStruct);
    StructureCopy(&tempStruct,pStructure);
    //for each mutant, build the rotamer-tree
    IntArray mutantArray,rotamersArray;
    IntArrayCreate(&mutantArray,0);
    IntArrayCreate(&rotamersArray,0);
    for(int cycle=0; cycle<StringArrayGetCount(&mutants[mutantIndex]); cycle++){
      char mutstr[10];
      char aa1, chn, aa2;
      int posInChain;
      strcpy(mutstr, StringArrayGet(&mutants[mutantIndex], cycle));
      sscanf(mutstr, "%c%c%d%c", &aa1, &chn, &posInChain, &aa2);
      int chainIndex = -1, residueIndex = -1;
      char chainname[MAX_LENGTH_CHAIN_NAME]; chainname[0] = chn; chainname[1] = '\0';
      StructureFindChainIndex(&tempStruct,chainname,&chainIndex);
      if(chainIndex==-1){
        printf("in file %s function %s() line %d, cannot find mutation %s\n", __FILE__, __FUNCTION__, __LINE__, mutstr);
        exit(ValueError);
      }
      ChainFindResidueByPosInChain(StructureGetChain(&tempStruct, chainIndex), posInChain, &residueIndex);
      if(residueIndex==-1){
        printf("in file %s function %s() line %d, cannot find mutation %s\n", __FILE__, __FUNCTION__, __LINE__, mutstr);
        exit(ValueError);
      }
      char mutaatype[MAX_LENGTH_RESIDUE_NAME];
      OneLetterAAToThreeLetterAA(aa2, mutaatype);
      StringArray designType, patchType;
      StringArrayCreate(&designType);
      StringArrayCreate(&patchType);
      // for histidine, the default mutaatype is HSD, we need to add HSE
      StringArrayAppend(&designType, mutaatype); StringArrayAppend(&patchType, "");
      if(aa2=='H'){StringArrayAppend(&designType, "HSE"); StringArrayAppend(&patchType, "");}
      ProteinSiteBuildMutatedRotamers(&tempStruct,chainIndex,residueIndex,rotlib,atomParams,resiTopos,&designType,&patchType);
      IntArrayAppend(&mutantArray, chainIndex);
      IntArrayAppend(&mutantArray, residueIndex);
      IntArrayAppend(&rotamersArray,chainIndex);
      IntArrayAppend(&rotamersArray,residueIndex);
      StringArrayDestroy(&designType);
      StringArrayDestroy(&patchType);
    }

    // for each mutant, find the surrounding residues and build the wild-type rotamer-tree
    for(int ii=0; ii<IntArrayGetLength(&mutantArray); ii+=2){
      int chainIndex = IntArrayGet(&mutantArray,ii);
      int resiIndex = IntArrayGet(&mutantArray,ii+1);
      Residue *pResi1 = ChainGetResidue(StructureGetChain(&tempStruct, chainIndex), resiIndex);
      for(int j = 0; j < StructureGetChainCount(&tempStruct); ++j){
        Chain* pChain = StructureGetChain(&tempStruct,j);
        for(int k=0; k<ChainGetResidueCount(pChain); k++){
          Residue* pResi2 = ChainGetResidue(pChain,k);
          if(AtomArrayCalcMinDistance(&pResi1->atoms,&pResi2->atoms)<ENERGY_DISTANCE_CUTOFF){
            if(pResi2->designSiteType==Type_ResidueDesignType_Fixed){
              ProteinSiteBuildWildtypeRotamers(&tempStruct,j,k,rotlib,atomParams,resiTopos);
              ProteinSiteAddCrystalRotamer(&tempStruct,j,k,resiTopos);
              IntArrayAppend(&rotamersArray,j);
              IntArrayAppend(&rotamersArray,k);
            }
          }
        }
      }
    }

    // optimization rotamers sequentially
    printf("EvoEF Building mutation model %d, the following sites will be optimized:\n",mutantIndex+1);
    //IntArrayShow(&rotamersArray);
    //printf("\n");
    printf("chnIndex resIndex (both of them starts from zero on the chain)\n");
    for(int ii=0;ii<IntArrayGetLength(&rotamersArray);ii+=2){
      printf("%8d %8d\n",IntArrayGet(&rotamersArray,ii),IntArrayGet(&rotamersArray,ii+1));
    }
    for(int cycle=0; cycle<3; cycle++){
      printf("optimization cycle %d ... \n",cycle+1);
      for(int ii=0; ii<IntArrayGetLength(&rotamersArray); ii+=2){
        int chainIndex = IntArrayGet(&rotamersArray, ii);
        int resiIndex = IntArrayGet(&rotamersArray, ii+1);
        ProteinSiteOptimizeRotamer(&tempStruct, chainIndex, resiIndex);
      }
    }
    IntArrayDestroy(&mutantArray);
    IntArrayDestroy(&rotamersArray);
    //remember to delete rotamers for previous mutant
    StructureRemoveAllDesignSites(&tempStruct);

    char modelfile[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    if(pdbid!=NULL)
      sprintf(modelfile,"%s_Model_%04d.pdb",pdbid,mutantIndex+1);
    else
      sprintf(modelfile,"EvoEF_Model_%04d.pdb",mutantIndex+1);
    FILE* pf=fopen(modelfile,"w");
    fprintf(pf,"REMARK EvoEF generated pdb file\n");
    fprintf(pf,"REMARK Output generated by EvoEF <BuildMutant>\n");
    StructureShowInPDBFormat(&tempStruct,TRUE,pf);
    fclose(pf);
    StructureDestroy(&tempStruct);
  }

  return Success;
}


int EVOEF_BuildMutantWithBBdepRotLib(Structure* pStructure, char* mutantfile, BBdepRotamerLib* pBBdepRotLib, AtomParamsSet* atomParams,ResiTopoSet* resiTopos, char* pdbid){
  FileReader fr;
  FileReaderCreate(&fr, mutantfile);
  int mutantcount = FileReaderGetLineCount(&fr);
  if(mutantcount<=0){
    printf("in %s %s %d, no mutation found in the mutant file\n",__FILE__,__FUNCTION__,__LINE__);
    exit(DataNotExistError);
  }
  StringArray* mutants = (StringArray*)malloc(sizeof(StringArray)*mutantcount);
  char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  int mutantIndex=0;
  while(!FAILED(FileReaderGetNextLine(&fr, line))){
    StringArrayCreate(&mutants[mutantIndex]);
    StringArraySplitString(&mutants[mutantIndex], line, ',');
    char lastMutant[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    int lastmutindex = StringArrayGetCount(&mutants[mutantIndex])-1;
    strcpy(lastMutant, StringArrayGet(&mutants[mutantIndex], lastmutindex));
    //deal with the last char of the last single mutant
    if((!isdigit(lastMutant[strlen(lastMutant)-1])) && !isalpha(lastMutant[strlen(lastMutant)-1])){
      lastMutant[strlen(lastMutant)-1] = '\0';
    }
    StringArraySet(&mutants[mutantIndex], lastmutindex, lastMutant);
    mutantIndex++;
  }
  FileReaderDestroy(&fr);

  for(int mutantIndex = 0;mutantIndex<mutantcount;mutantIndex++){
    Structure tempStruct;
    StructureCreate(&tempStruct);
    StructureCopy(&tempStruct,pStructure);
    //for each mutant, build the rotamer-tree
    IntArray mutatedArray,rotamericArray;
    IntArrayCreate(&mutatedArray,0);
    IntArrayCreate(&rotamericArray,0);
    for(int posIndex=0; posIndex<StringArrayGetCount(&mutants[mutantIndex]); posIndex++){
      char mutstr[10];
      char aa1,chn,aa2;
      int posInChain;
      strcpy(mutstr,StringArrayGet(&mutants[mutantIndex],posIndex));
      sscanf(mutstr,"%c%c%d%c",&aa1,&chn,&posInChain,&aa2);
      int chainIndex = -1, residueIndex = -1;
      char chainname[MAX_LENGTH_CHAIN_NAME]; chainname[0]=chn; chainname[1]='\0';
      StructureFindChainIndex(&tempStruct,chainname,&chainIndex);
      if(chainIndex==-1){
        printf("in %s %s %d, cannot find mutation %s\n", __FILE__, __FUNCTION__, __LINE__, mutstr);
        exit(ValueError);
      }
      ChainFindResidueByPosInChain(StructureGetChain(&tempStruct, chainIndex), posInChain, &residueIndex);
      if(residueIndex==-1){
        printf("in %s %s %d, cannot find mutation %s\n", __FILE__, __FUNCTION__, __LINE__, mutstr);
        exit(ValueError);
      }
      char mutaatype[MAX_LENGTH_RESIDUE_NAME];
      OneLetterAAToThreeLetterAA(aa2, mutaatype);
      StringArray designType, patchType;
      StringArrayCreate(&designType);
      StringArrayCreate(&patchType);
      //for histidine, the default mutaatype is HSD, we need to add HSE
      StringArrayAppend(&designType, mutaatype); StringArrayAppend(&patchType, "");
      if(aa2=='H'){StringArrayAppend(&designType, "HSE"); StringArrayAppend(&patchType, "");}
      ProteinSiteBuildMutatedRotamersByBBdepRotLib(&tempStruct,chainIndex,residueIndex,pBBdepRotLib,atomParams,resiTopos,&designType,&patchType);
      IntArrayAppend(&mutatedArray, chainIndex);
      IntArrayAppend(&mutatedArray, residueIndex);
      IntArrayAppend(&rotamericArray,chainIndex);
      IntArrayAppend(&rotamericArray,residueIndex);
      StringArrayDestroy(&designType);
      StringArrayDestroy(&patchType);
    }

    //build rotamers for surrounding residues
    for(int ii=0; ii<IntArrayGetLength(&mutatedArray); ii+=2){
      int chainIndex = IntArrayGet(&mutatedArray,ii);
      int resiIndex = IntArrayGet(&mutatedArray,ii+1);
      Residue *pResi1 = ChainGetResidue(StructureGetChain(&tempStruct, chainIndex), resiIndex);
      for(int j = 0; j < StructureGetChainCount(&tempStruct); ++j){
        Chain* pChain = StructureGetChain(&tempStruct,j);
        for(int k=0; k<ChainGetResidueCount(pChain); k++){
          Residue* pResi2 = ChainGetResidue(pChain,k);
          if(AtomArrayCalcMinDistance(&pResi1->atoms,&pResi2->atoms)<ENERGY_DISTANCE_CUTOFF){
            if(pResi2->designSiteType==Type_ResidueDesignType_Fixed){
              ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&tempStruct,j,k,pBBdepRotLib,atomParams,resiTopos);
              ProteinSiteAddCrystalRotamerByBBdepRotLib(&tempStruct,j,k,resiTopos,pBBdepRotLib);
              IntArrayAppend(&rotamericArray,j);
              IntArrayAppend(&rotamericArray,k);
            }
          }
        }
      }
    }

    //optimization rotamers sequentially
    printf("EvoEF Building mutation model %d, the following sites will be optimized:\n",mutantIndex+1);
    printf("chnIndex resIndex (both of them starts from zero on the chain)\n");
    for(int ii=0;ii<IntArrayGetLength(&rotamericArray);ii+=2){
      printf("%8d %8d\n",IntArrayGet(&rotamericArray,ii),IntArrayGet(&rotamericArray,ii+1));
    }
    for(int posIndex=0; posIndex<3; posIndex++){
      printf("optimization cycle %d ... \n",posIndex+1);
      for(int ii=0; ii<IntArrayGetLength(&rotamericArray); ii+=2){
        int chainIndex = IntArrayGet(&rotamericArray, ii);
        int resiIndex = IntArrayGet(&rotamericArray, ii+1);
        ProteinSiteOptimizeRotamerWithBBdepRotLib(&tempStruct,chainIndex,resiIndex,pBBdepRotLib);
      }
    }
    IntArrayDestroy(&mutatedArray);
    IntArrayDestroy(&rotamericArray);
    //remember to delete rotamers for previous mutant
    StructureRemoveAllDesignSites(&tempStruct);

    char modelfile[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    if(pdbid!=NULL) sprintf(modelfile,"%s_Model_%04d.pdb",pdbid,mutantIndex+1);
    else sprintf(modelfile,"EvoEF_Model_%04d.pdb",mutantIndex+1);
    FILE* pf=fopen(modelfile,"w");
    fprintf(pf,"REMARK EvoEF generated pdb file\n");
    fprintf(pf,"REMARK Output generated by EvoEF <BuildMutant>\n");
    StructureShowInPDBFormat(&tempStruct,TRUE,pf);
    fclose(pf);
    StructureDestroy(&tempStruct);
  }

  return Success;
}


int EVOEF_RepairStructure(Structure* pStructure, BBindRotamerLib* pBBindRotLib, AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* pdbid){
  for(int cycle=0; cycle<1; cycle++){
    printf("EvoEF Repairing PDB: optimization cycle %d ... \n",cycle+1);
    for(int i=0; i<StructureGetChainCount(pStructure); ++i){
      Chain* pChain = StructureGetChain(pStructure, i);
      for(int j=0; j<ChainGetResidueCount(pChain); j++){
        Residue* pResi = ChainGetResidue(pChain, j);
        if(strcmp(ResidueGetName(pResi),"ALA")==0 || strcmp(ResidueGetName(pResi),"GLY")==0) continue;
        //skip CYS which may form disulfide bonds
        if(strcmp(ResidueGetName(pResi),"CYS")==0) continue;
        if(strcmp(ResidueGetName(pResi),"ASN")==0||strcmp(ResidueGetName(pResi),"GLN")==0||strcmp(ResidueGetName(pResi),"HSD")==0||strcmp(ResidueGetName(pResi),"HSE")==0){
          printf("We will flip residue %s%d%c to optimize hbond\n", ResidueGetChainName(pResi),ResidueGetPosInChain(pResi),ThreeLetterAAToOneLetterAA(ResidueGetName(pResi)));
          ProteinSiteAddCrystalRotamer(pStructure,i,j,resiTopos);
          ProteinSiteBuildFlippedCrystalRotamer(pStructure,i,j,resiTopos);
          ProteinSiteOptimizeRotamerHBondEnergy(pStructure,i,j);
        }
        else if(strcmp(ResidueGetName(pResi),"SER")==0 || strcmp(ResidueGetName(pResi),"THR")==0 || strcmp(ResidueGetName(pResi),"TYR")==0){
          printf("We will rotate hydroxyl group of residue %s%d%c to optimize hbond\n", ResidueGetChainName(pResi),ResidueGetPosInChain(pResi),ThreeLetterAAToOneLetterAA(ResidueGetName(pResi)));
          ProteinSiteAddCrystalRotamer(pStructure,i,j,resiTopos);
          ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
          ProteinSiteOptimizeRotamerHBondEnergy(pStructure,i,j);
        }
        if(TRUE){
          printf("We optimize side chain of residue %s%d%c\n", ResidueGetChainName(pResi),ResidueGetPosInChain(pResi),ThreeLetterAAToOneLetterAA(ResidueGetName(pResi)));
          ProteinSiteBuildWildtypeRotamers(pStructure,i,j,pBBindRotLib,atomParams,resiTopos);
          ProteinSiteAddCrystalRotamer(pStructure,i,j,resiTopos);
          ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
          ProteinSiteOptimizeRotamer(pStructure,i,j);
          //ProteinSiteOptimizeRotamerLocally(pStructure,i,j, 1.0);
        }
        ProteinSiteRemoveDesignSite(pStructure,i,j);
      }
    }
  }

  //output the repaired structure
  char modelfile[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  if(pdbid!=NULL){sprintf(modelfile,"%s_Repair.pdb",pdbid);}
  else{strcpy(modelfile,"EvoEF_Repair.pdb");}
  FILE* pf=fopen(modelfile,"w");
  fprintf(pf,"REMARK EvoEF generated pdb file\n");
  fprintf(pf,"REMARK Output generated by EvoEF <RepairStructure>\n");
  StructureShowInPDBFormat(pStructure,TRUE,pf);
  fclose(pf);

  return Success;
}


int EVOEF_RepairStructureWithBBdepRotLib(Structure* pStructure,BBdepRotamerLib* pBBdepRotLib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* pdbid){
  for(int cycle=0; cycle<1; cycle++){
    printf("EvoEF Repairing PDB: optimization cycle %d ... \n",cycle+1);
    for(int i=0; i<StructureGetChainCount(pStructure); ++i){
      Chain* pChain = StructureGetChain(pStructure, i);
      for(int j=0; j<ChainGetResidueCount(pChain); j++){
        Residue* pResi = ChainGetResidue(pChain, j);
        if(strcmp(ResidueGetName(pResi),"ALA")==0 || strcmp(ResidueGetName(pResi),"GLY")==0) continue;
        //skip CYS which may form disulfide bonds
        if(strcmp(ResidueGetName(pResi),"CYS")==0) continue;
        if(strcmp(ResidueGetName(pResi),"ASN")==0||strcmp(ResidueGetName(pResi),"GLN")==0||strcmp(ResidueGetName(pResi),"HSD")==0||strcmp(ResidueGetName(pResi),"HSE")==0){
          printf("We will flip residue %s%d%c to optimize hbond\n", ResidueGetChainName(pResi),ResidueGetPosInChain(pResi),ThreeLetterAAToOneLetterAA(ResidueGetName(pResi)));
          ProteinSiteAddCrystalRotamer(pStructure,i,j,resiTopos);
          ProteinSiteBuildFlippedCrystalRotamer(pStructure,i,j,resiTopos);
          ProteinSiteOptimizeRotamerHBondEnergy(pStructure,i,j);
        }
        else if(strcmp(ResidueGetName(pResi),"SER")==0 || strcmp(ResidueGetName(pResi),"THR")==0 || strcmp(ResidueGetName(pResi),"TYR")==0){
          printf("We will rotate hydroxyl group of residue %s%d%c to optimize hbond\n", ResidueGetChainName(pResi),ResidueGetPosInChain(pResi),ThreeLetterAAToOneLetterAA(ResidueGetName(pResi)));
          ProteinSiteAddCrystalRotamer(pStructure,i,j,resiTopos);
          ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
          ProteinSiteOptimizeRotamerHBondEnergy(pStructure,i,j);
        }
        if(TRUE){
          printf("We optimize side chain of residue %s%d%c\n", ResidueGetChainName(pResi),ResidueGetPosInChain(pResi),ThreeLetterAAToOneLetterAA(ResidueGetName(pResi)));
          ProteinSiteBuildWildtypeRotamersByBBdepRotLib(pStructure,i,j,pBBdepRotLib,atomParams,resiTopos);
          ProteinSiteAddCrystalRotamerByBBdepRotLib(pStructure,i,j,resiTopos,pBBdepRotLib);
          ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
          ProteinSiteOptimizeRotamerWithBBdepRotLib(pStructure,i,j,pBBdepRotLib);
        }
        ProteinSiteRemoveDesignSite(pStructure,i,j);
      }
    }
  }

  //output the repaired structure
  char modelfile[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  if(pdbid!=NULL){sprintf(modelfile,"%s_Repair.pdb",pdbid);}
  else{strcpy(modelfile,"EvoEF_Repair.pdb");}
  FILE* pf=fopen(modelfile,"w");
  fprintf(pf,"REMARK EvoEF generated pdb file\n");
  fprintf(pf,"REMARK Output generated by EvoEF <RepairStructure>\n");
  StructureShowInPDBFormat(pStructure,TRUE,pf);
  fclose(pf);

  return Success;
}


int EVOEF_WriteStructureToFile(Structure* pStructure, char* pdbfile){
  FILE* pf=fopen(pdbfile,"w");
  if(pf!=NULL){
    StructureShowInPDBFormat(pStructure,TRUE,pf);
    fclose(pf);
  }
  else{
    printf("failed to open file for writing structure coordinates\n");
    return IOError;
  }
  return Success;
}

int EVOEF_AddHydrogen(Structure* pStructure, char* pdbid){
  //polar hydrogens are automatically added, so we just output the repaired structure
  char modelfile[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  if(pdbid!=NULL){sprintf(modelfile,"%s_PolH.pdb",pdbid);}
  else{strcpy(modelfile,"EvoEF_PolH.pdb");}
  FILE* pf=fopen(modelfile,"w");
  fprintf(pf,"REMARK EvoEF generated pdb file\n");
  fprintf(pf,"REMARK Output generated by EvoEF <AddHydrogen>\n");
  StructureShowInPDBFormat(pStructure,TRUE,pf);
  fclose(pf);
  return Success;
}


int EVOEF_OptimizeHydrogen(Structure* pStructure, AtomParamsSet* atomParams,ResiTopoSet* resiTopos, char* pdbid){
  for(int cycle=0; cycle<1; cycle++){
    printf("EvoEF Repairing PDB: optimization cycle %d ... \n",cycle+1);
    for(int i=0; i<StructureGetChainCount(pStructure); ++i){
      Chain* pChain = StructureGetChain(pStructure, i);
      for(int j=0; j<ChainGetResidueCount(pChain); j++){
        Residue* pResi = ChainGetResidue(pChain, j);
        if(strcmp(ResidueGetName(pResi),"SER")==0 || strcmp(ResidueGetName(pResi),"THR")==0 || strcmp(ResidueGetName(pResi),"TYR")==0){
          printf("We will rotate hydroxyl group of residue %s%d%c to optimize hbond\n", ResidueGetChainName(pResi),ResidueGetPosInChain(pResi),ThreeLetterAAToOneLetterAA(ResidueGetName(pResi)));
          ProteinSiteAddCrystalRotamer(pStructure,i,j,resiTopos);
          ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
          ProteinSiteOptimizeRotamerHBondEnergy(pStructure,i,j);
          ProteinSiteRemoveDesignSite(pStructure,i,j);
        }
      }
    }
  }

  //output the repaired structure
  char modelfile[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  if(pdbid!=NULL){sprintf(modelfile,"%s_OptH.pdb",pdbid);}
  else{strcpy(modelfile,"EvoEF_OptH.pdb");}
  FILE* pf=fopen(modelfile,"w");
  fprintf(pf,"REMARK EvoEF generated pdb file\n");
  fprintf(pf,"REMARK Output generated by EvoEF <OptimizeHydrogen>\n");
  StructureShowInPDBFormat(pStructure,TRUE,pf);
  fclose(pf);

  return Success;
}


int EVOEF_StructureComputeResidueInteractionWithFixedSurroundingResidues(Structure *pStructure, int chainIndex, int residueIndex){
  double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
  // if the structure is composed of several chains, the residue position could be different in the whole structure from that in the separate chain
  StructureComputeResiduePosition(pStructure);
  Chain *pChainI=StructureGetChain(pStructure, chainIndex);
  Residue *pResIR= ChainGetResidue(pChainI, residueIndex);
  //step 1: find out residues within 5 angstroms to the design site of interest;
  int surroundingResiNum = 0;
  Residue **ppSurroundingResidues = NULL;
  for(int i = 0; i < pStructure->chainNum; i++){
    Chain *pChainI = pStructure->chains + i;
    for(int j = 0; j < pChainI->residueNum; j++){
      Residue *pResi2 = pChainI->residues + j;
      if(strcmp(pResIR->name, pResi2->name) == 0 && pResIR->posInChain == pResi2->posInChain){
        continue;
      }
      else if(AtomArrayCalcMinDistance(&pResIR->atoms, &pResi2->atoms)<ENERGY_DISTANCE_CUTOFF){
        surroundingResiNum++;
        ppSurroundingResidues = (Residue **)realloc(ppSurroundingResidues, sizeof(Residue*)*surroundingResiNum);
        ppSurroundingResidues[surroundingResiNum-1] = pResi2;
      }
    }
  }
  // calculate energy between residue IR and other residues
  EVOEF_AminoAcidReferenceEnergy(pResIR->name,energyTerms);
  EVOEF_EnergyResidueIntraEnergy(pResIR, energyTerms);
  for(int is = 0; is < surroundingResiNum; is++){
    Residue *pResIS = ppSurroundingResidues[is];
    if(strcmp(ResidueGetChainName(pResIR),ResidueGetChainName(pResIS))==0){
      if(pResIR->posInChain == pResIS->posInChain-1){
        EVOEF_EnergyResidueAndNextResidue(pResIR,pResIS,energyTerms);
      }
      else if(pResIR->posInChain == pResIS->posInChain+1){
        EVOEF_EnergyResidueAndNextResidue(pResIS,pResIR,energyTerms);
      }
      else{
        EVOEF_EnergyResidueAndOtherResidueSameChain(pResIR,pResIS,energyTerms);
      }
    }
    else{
      if(ChainGetType(StructureFindChainByName(pStructure,ResidueGetName(pResIS)))==Type_Chain_SmallMol){
        EVOEF_EnergyResidueAndLigandResidue(pResIR,pResIS,energyTerms);
      }
      else{
        EVOEF_EnergyResidueAndOtherResidueDifferentChain(pResIR,pResIS,energyTerms);
      }
    }
  }

  if(TRUE){
    //for debug: energy details, not weighted
    printf("Energy details between residue %s%d%c and fixed surrounding residues:\n", ResidueGetChainName(pResIR),ResidueGetPosInChain(pResIR),ThreeLetterAAToOneLetterAA(ResidueGetName(pResIR)));
    printf("reference_ALA         =            %12.6f\n", energyTerms[ 1]);
    printf("reference_CYS         =            %12.6f\n", energyTerms[ 2]);
    printf("reference_ASP         =            %12.6f\n", energyTerms[ 3]);
    printf("reference_GLU         =            %12.6f\n", energyTerms[ 4]);
    printf("reference_PHE         =            %12.6f\n", energyTerms[ 5]);
    printf("reference_GLY         =            %12.6f\n", energyTerms[ 6]);
    printf("reference_HIS         =            %12.6f\n", energyTerms[ 7]);
    printf("reference_ILE         =            %12.6f\n", energyTerms[ 8]);
    printf("reference_LYS         =            %12.6f\n", energyTerms[ 9]);
    printf("reference_LEU         =            %12.6f\n", energyTerms[10]);
    printf("reference_MET         =            %12.6f\n", energyTerms[11]);
    printf("reference_ASN         =            %12.6f\n", energyTerms[12]);
    printf("reference_PRO         =            %12.6f\n", energyTerms[13]);
    printf("reference_GLN         =            %12.6f\n", energyTerms[14]);
    printf("reference_ARG         =            %12.6f\n", energyTerms[15]);
    printf("reference_SER         =            %12.6f\n", energyTerms[16]);
    printf("reference_THR         =            %12.6f\n", energyTerms[17]);
    printf("reference_VAL         =            %12.6f\n", energyTerms[18]);
    printf("reference_TRP         =            %12.6f\n", energyTerms[19]);
    printf("reference_TYR         =            %12.6f\n", energyTerms[20]);

    printf("intraR_vdwatt         =            %12.6f\n", energyTerms[21]);
    printf("intraR_vdwrep         =            %12.6f\n", energyTerms[22]);
    printf("intraR_electr         =            %12.6f\n", energyTerms[23]);
    printf("intraR_deslvP         =            %12.6f\n", energyTerms[24]);
    printf("intraR_deslvH         =            %12.6f\n", energyTerms[25]);
    printf("intraR_hbscbb_dis     =            %12.6f\n", energyTerms[26]);
    printf("intraR_hbscbb_the     =            %12.6f\n", energyTerms[27]);
    printf("intraR_hbscbb_phi     =            %12.6f\n", energyTerms[28]);

    printf("interS_vdwatt         =            %12.6f\n", energyTerms[31]);
    printf("interS_vdwrep         =            %12.6f\n", energyTerms[32]);
    printf("interS_electr         =            %12.6f\n", energyTerms[33]);
    printf("interS_deslvP         =            %12.6f\n", energyTerms[34]);
    printf("interS_deslvH         =            %12.6f\n", energyTerms[35]);
    printf("interS_ssbond         =            %12.6f\n", energyTerms[36]);
    printf("interS_hbbbbb_dis     =            %12.6f\n", energyTerms[41]);
    printf("interS_hbbbbb_the     =            %12.6f\n", energyTerms[42]);
    printf("interS_hbbbbb_phi     =            %12.6f\n", energyTerms[43]);
    printf("interS_hbscbb_dis     =            %12.6f\n", energyTerms[44]);
    printf("interS_hbscbb_the     =            %12.6f\n", energyTerms[45]);
    printf("interS_hbscbb_phi     =            %12.6f\n", energyTerms[46]);
    printf("interS_hbscsc_dis     =            %12.6f\n", energyTerms[47]);
    printf("interS_hbscsc_the     =            %12.6f\n", energyTerms[48]);
    printf("interS_hbscsc_phi     =            %12.6f\n", energyTerms[49]);

    printf("interD_vdwatt         =            %12.6f\n", energyTerms[51]);
    printf("interD_vdwrep         =            %12.6f\n", energyTerms[52]);
    printf("interD_electr         =            %12.6f\n", energyTerms[53]);
    printf("interD_deslvP         =            %12.6f\n", energyTerms[54]);
    printf("interD_deslvH         =            %12.6f\n", energyTerms[55]);
    printf("interD_ssbond         =            %12.6f\n", energyTerms[56]);
    printf("interD_hbbbbb_dis     =            %12.6f\n", energyTerms[61]);
    printf("interD_hbbbbb_the     =            %12.6f\n", energyTerms[62]);
    printf("interD_hbbbbb_phi     =            %12.6f\n", energyTerms[63]);
    printf("interD_hbscbb_dis     =            %12.6f\n", energyTerms[64]);
    printf("interD_hbscbb_the     =            %12.6f\n", energyTerms[65]);
    printf("interD_hbscbb_phi     =            %12.6f\n", energyTerms[66]);
    printf("interD_hbscsc_dis     =            %12.6f\n", energyTerms[67]);
    printf("interD_hbscsc_the     =            %12.6f\n", energyTerms[68]);
    printf("interD_hbscsc_phi     =            %12.6f\n", energyTerms[69]);
    //total energy: weighted
    EnergyTermWeighting(energyTerms);
    printf("----------------------------------------------------\n");
    printf("Total                 =            %12.6f\n\n", energyTerms[0]);
  }

  return Success;
}

//////////////////////////////////////////////////////////////////////////////
//The following are new program functions
//////////////////////////////////////////////////////////////////////////////
int EVOEF_ComputeRotamersEnergy(Structure* pStructure,BBindRotamerLib* pBBindRotLib,AAppTable* pAAppTable,RamaTable* pRama,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* pdbid){
  char energyfile[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  sprintf(energyfile,"%s_rotenergy.txt",pdbid);
  FILE* fp=fopen(energyfile,"w");
  for(int i=0; i<StructureGetChainCount(pStructure); ++i){
    Chain* pChain = StructureGetChain(pStructure, i);
    if(ChainGetType(pChain)!=Type_Chain_Protein) continue;
    for(int j=0; j<ChainGetResidueCount(pChain); j++){
      Residue* pResi = ChainGetResidue(pChain, j);
      if(FLAG_PPI==TRUE){
        //check if residue is an interface residue
        for(int k=0;k<StructureGetChainCount(pStructure);k++){
          if(k==i) continue;
          Chain* pChainK=StructureGetChain(pStructure,k);
          for(int s=0;s<ChainGetResidueCount(pChainK);s++){
            Residue* pResiKS=ChainGetResidue(pChainK,s);
            if(AtomArrayCalcMinDistance(&pResi->atoms,&pResiKS->atoms)<ENERGY_DISTANCE_CUTOFF){
              ProteinSiteBuildAllRotamers(pStructure,i,j,pBBindRotLib,atomParams,resiTopos);
              if(FLAG_ADD_CRYSTAL_ROT) ProteinSiteAddCrystalRotamer(pStructure,i,j,resiTopos);
              if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
              ProteinSiteCalcRotamersEnergy(pStructure,pAAppTable,pRama,i,j,fp);
              ProteinSiteRemoveDesignSite(pStructure,i,j);
            }
          }
        }
      }
      else if(FLAG_ENZYME==TRUE || FLAG_PROT_LIG==TRUE){
        Residue* pSmallMol=NULL;
        StructureFindSmallMol(pStructure,&pSmallMol);
        if(AtomArrayCalcMinDistance(&pResi->atoms,&pSmallMol->atoms)<ENERGY_DISTANCE_CUTOFF){
          ProteinSiteBuildAllRotamers(pStructure,i,j,pBBindRotLib,atomParams,resiTopos);
          if(FLAG_ADD_CRYSTAL_ROT) ProteinSiteAddCrystalRotamer(pStructure,i,j,resiTopos);
          if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
          ProteinSiteCalcRotamersEnergy(pStructure,pAAppTable,pRama,i,j,fp);
          ProteinSiteRemoveDesignSite(pStructure,i,j);
        }
      }
      else{
        ProteinSiteBuildAllRotamers(pStructure,i,j,pBBindRotLib,atomParams,resiTopos);
        if(FLAG_ADD_CRYSTAL_ROT) ProteinSiteAddCrystalRotamer(pStructure,i,j,resiTopos);
        if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
        ProteinSiteCalcRotamersEnergy(pStructure,pAAppTable,pRama,i,j,fp);
        ProteinSiteRemoveDesignSite(pStructure,i,j);
      }
    }
  }
  fclose(fp);

  return Success;
}


int EVOEF_ComputeWildtypeRotamersEnergy(Structure* pStructure,BBindRotamerLib* pBBindRotLib,AAppTable* pAAppTable,RamaTable* pRama,AtomParamsSet* atomParams,ResiTopoSet* resiTopos, char* pdbid){
  char energyfile[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  sprintf(energyfile,"%s_rotenergy.txt",pdbid);
  FILE* fp=fopen(energyfile,"w");
  for(int i=0; i<StructureGetChainCount(pStructure); ++i){
    Chain* pChain = StructureGetChain(pStructure, i);
    if(ChainGetType(pChain)!=Type_Chain_Protein) continue;
    for(int j=0; j<ChainGetResidueCount(pChain); j++){
      Residue* pResi = ChainGetResidue(pChain, j);
      if(FLAG_PPI==TRUE){
        //check if residue is an interface residue
        for(int k=0;k<StructureGetChainCount(pStructure);k++){
          if(k==i) continue;
          Chain* pChainK=StructureGetChain(pStructure,k);
          for(int s=0;s<ChainGetResidueCount(pChainK);s++){
            Residue* pResiKS=ChainGetResidue(pChainK,s);
            if(AtomArrayCalcMinDistance(&pResi->atoms,&pResiKS->atoms)<ENERGY_DISTANCE_CUTOFF){
              ProteinSiteBuildWildtypeRotamers(pStructure,i,j,pBBindRotLib,atomParams,resiTopos);
              if(FLAG_ADD_CRYSTAL_ROT) ProteinSiteAddCrystalRotamer(pStructure,i,j,resiTopos);
              if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
              ProteinSiteCalcRotamersEnergy(pStructure,pAAppTable,pRama,i,j,fp);
              ProteinSiteRemoveDesignSite(pStructure,i,j);
            }
          }
        }
      }
      else if(FLAG_ENZYME==TRUE || FLAG_PROT_LIG==TRUE){
        Residue* pSmallMol=NULL;
        StructureFindSmallMol(pStructure,&pSmallMol);
        if(AtomArrayCalcMinDistance(&pResi->atoms,&pSmallMol->atoms)<ENERGY_DISTANCE_CUTOFF){
          ProteinSiteBuildWildtypeRotamers(pStructure,i,j,pBBindRotLib,atomParams,resiTopos);
          if(FLAG_ADD_CRYSTAL_ROT) ProteinSiteAddCrystalRotamer(pStructure,i,j,resiTopos);
          if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
          ProteinSiteCalcRotamersEnergy(pStructure,pAAppTable,pRama,i,j,fp);
          ProteinSiteRemoveDesignSite(pStructure,i,j);
        }
      }
      else{
        ProteinSiteBuildWildtypeRotamers(pStructure,i,j,pBBindRotLib,atomParams,resiTopos);
        if(FLAG_ADD_CRYSTAL_ROT) ProteinSiteAddCrystalRotamer(pStructure,i,j,resiTopos);
        if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
        ProteinSiteCalcRotamersEnergy(pStructure,pAAppTable,pRama,i,j,fp);
        ProteinSiteRemoveDesignSite(pStructure,i,j);
      }
    }
  }
  fclose(fp);

  return Success;
}


int EVOEF_StructureFindInterfaceResidues(Structure *pStructure, double cutoff,char* outputfile){
  if(pStructure->chainNum < 2){
    printf("there is only one chain in the whole structure, no protein-protein interface found\n");
    exit(ValueError);
  }

  IntArray* chainArrays=(IntArray*)malloc(sizeof(IntArray)*StructureGetChainCount(pStructure));
  for(int i=0;i<StructureGetChainCount(pStructure);i++){
    Chain* pChainI=StructureGetChain(pStructure,i);
    IntArrayCreate(&chainArrays[i],pChainI->residueNum);
    for(int j=0;j<pChainI->residueNum;j++){
      IntArraySet(&chainArrays[i],j,0);
    }
  }

  for(int i = 0; i < pStructure->chainNum; i++){
    Chain *pChainI = StructureGetChain(pStructure, i);
    for(int k = i+1; k < pStructure->chainNum; k++){
      Chain *pChainK = StructureGetChain(pStructure, k);
      for(int j = 0; j < pChainI->residueNum; j++){
        Residue* pResiIJ = ChainGetResidue(pChainI, j);
        //if(IntArrayGet(&chainArrays[i],j)==1)continue;
        for(int s = 0; s < pChainK->residueNum; s++){
          Residue* pResiKS = ChainGetResidue(pChainK, s);
          //if(IntArrayGet(&chainArrays[k],s)==1)continue;
          if(AtomArrayCalcMinDistance(&pResiIJ->atoms, &pResiKS->atoms) < cutoff){
            IntArraySet(&chainArrays[i], j, 1);
            IntArraySet(&chainArrays[k], s, 1);
          }
        }
      }
    }
  }

  FILE* fout=fopen(outputfile,"w");
  printf("Interface residues: \n");
  for(int j=0;j<StructureGetChainCount(pStructure);j++){
    Chain* pChain=StructureGetChain(pStructure,j);
    if(ChainGetType(pChain)==Type_Chain_Protein || ChainGetType(pChain)==Type_Chain_Nucleotide){
      for(int i=0; i<ChainGetResidueCount(pChain); i++){
        Residue* pResidue=ChainGetResidue(pChain,i);
        if(IntArrayGet(&chainArrays[j],i)==1){
          //printf("%c",ThreeLetterAAToOneLetterAA(ResidueGetName(pResidue)));
          printf("%s%d%c\n",ResidueGetChainName(pResidue),ResidueGetPosInChain(pResidue),ThreeLetterAAToOneLetterAA(ResidueGetName(pResidue)));
          fprintf(fout,"%c",ThreeLetterAAToOneLetterAA(ResidueGetName(pResidue)));
        }
        else{
          //printf("-");
          fprintf(fout,"-");
        }
      }
      //printf("\n");
      fprintf(fout,"\n");
    }
  }
  fclose(fout);

  for(int i=0;i<StructureGetChainCount(pStructure);i++){
    IntArrayDestroy(&chainArrays[i]);
  }
  free(chainArrays);
  chainArrays=NULL;

  return Success;
}


int EVOEF_StructureShowPhiPsi(Structure* pStructure,char* phipsifile){
  FILE* fout=NULL;
  if(phipsifile==NULL) fout=stdout;
  else fout=fopen(phipsifile,"w");
  fprintf(fout,"#phi-psi angles of protein residues: \n");
  for(int i=0;i<StructureGetChainCount(pStructure);i++){
    Chain* pChain=StructureGetChain(pStructure,i);
    if(ChainGetType(pChain)!=Type_Chain_Protein) continue;
    for(int j=0;j<ChainGetResidueCount(pChain);j++){
      Residue* pResi1=ChainGetResidue(pChain,j);
      fprintf(fout,"%s %s %4d %8.3f %8.3f\n",ResidueGetChainName(pResi1),ResidueGetName(pResi1),ResidueGetPosInChain(pResi1),pResi1->phipsi[0],pResi1->phipsi[1]);
    }
  }
  if(fout != stdout) fclose(fout);

  return Success;
}


int EVOEF_ComputeRotamersEnergyByBBdepRotLib(Structure* pStructure,BBdepRotamerLib* pBBdepRotLib,AAppTable* pAAppTable,RamaTable* pRama,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* pdbid){
  char energyfile[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  sprintf(energyfile,"%s_rotenergy.txt",pdbid);
  FILE* fp=fopen(energyfile,"w");
  for(int i=0; i<StructureGetChainCount(pStructure); ++i){
    Chain* pChain = StructureGetChain(pStructure, i);
    if(ChainGetType(pChain)!=Type_Chain_Protein) continue;
    for(int j=0; j<ChainGetResidueCount(pChain); j++){
      Residue* pResi = ChainGetResidue(pChain, j);
      if(FLAG_PPI==TRUE){
        //check if residue is an interface residue
        for(int k=0;k<StructureGetChainCount(pStructure);k++){
          if(k==i) continue;
          Chain* pChainK=StructureGetChain(pStructure,k);
          for(int s=0;s<ChainGetResidueCount(pChainK);s++){
            Residue* pResiKS=ChainGetResidue(pChainK,s);
            if(AtomArrayCalcMinDistance(&pResi->atoms,&pResiKS->atoms)<ENERGY_DISTANCE_CUTOFF){
              ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,i,j,pBBdepRotLib,atomParams,resiTopos);
              if(FLAG_ADD_CRYSTAL_ROT) ProteinSiteAddCrystalRotamerByBBdepRotLib(pStructure,i,j,resiTopos,pBBdepRotLib);
              if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
              ProteinSiteCalcRotamersEnergy(pStructure,pAAppTable,pRama,i,j,fp);
              ProteinSiteRemoveDesignSite(pStructure,i,j);
            }
          }
        }
      }
      else if(FLAG_ENZYME==TRUE || FLAG_PROT_LIG==TRUE){
        Residue* pSmallMol=NULL;
        StructureFindSmallMol(pStructure,&pSmallMol);
        if(AtomArrayCalcMinDistance(&pResi->atoms,&pSmallMol->atoms)<ENERGY_DISTANCE_CUTOFF){
          ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,i,j,pBBdepRotLib,atomParams,resiTopos);
          if(FLAG_ADD_CRYSTAL_ROT) ProteinSiteAddCrystalRotamerByBBdepRotLib(pStructure,i,j,resiTopos,pBBdepRotLib);
          if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
          ProteinSiteCalcRotamersEnergy(pStructure,pAAppTable,pRama,i,j,fp);
          ProteinSiteRemoveDesignSite(pStructure,i,j);
        }
      }
      else{
        ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,i,j,pBBdepRotLib,atomParams,resiTopos);
        if(FLAG_ADD_CRYSTAL_ROT) ProteinSiteAddCrystalRotamerByBBdepRotLib(pStructure,i,j,resiTopos,pBBdepRotLib);
        if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
        ProteinSiteCalcRotamersEnergy(pStructure,pAAppTable,pRama,i,j,fp);
        ProteinSiteRemoveDesignSite(pStructure,i,j);
      }
    }
  }
  fclose(fp);

  return Success;
}


int EVOEF_ComputeWildtypeRotamersEnergyByBBdepRotLib(Structure* pStructure,BBdepRotamerLib* pBBdepRotLib,AAppTable* pAAppTable,RamaTable* pRama,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* pdbid){
  char energyfile[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  sprintf(energyfile,"%s_rotenergy.txt",pdbid);
  FILE* fp=fopen(energyfile,"w");
  for(int i=0; i<StructureGetChainCount(pStructure); ++i){
    Chain* pChain = StructureGetChain(pStructure, i);
    if(ChainGetType(pChain)!=Type_Chain_Protein) continue;
    for(int j=0; j<ChainGetResidueCount(pChain); j++){
      Residue* pResi = ChainGetResidue(pChain, j);
      if(FLAG_PPI==TRUE){
        //check if residue is an interface residue
        for(int k=0;k<StructureGetChainCount(pStructure);k++){
          if(k==i) continue;
          Chain* pChainK=StructureGetChain(pStructure,k);
          for(int s=0;s<ChainGetResidueCount(pChainK);s++){
            Residue* pResiKS=ChainGetResidue(pChainK,s);
            if(AtomArrayCalcMinDistance(&pResi->atoms,&pResiKS->atoms)<ENERGY_DISTANCE_CUTOFF){
              ProteinSiteBuildWildtypeRotamersByBBdepRotLib(pStructure,i,j,pBBdepRotLib,atomParams,resiTopos);
              if(FLAG_ADD_CRYSTAL_ROT) ProteinSiteAddCrystalRotamerByBBdepRotLib(pStructure,i,j,resiTopos,pBBdepRotLib);
              if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
              ProteinSiteCalcRotamersEnergy(pStructure,pAAppTable,pRama,i,j,fp);
              ProteinSiteRemoveDesignSite(pStructure,i,j);
            }
          }
        }
      }
      else if(FLAG_ENZYME==TRUE || FLAG_PROT_LIG==TRUE){
        Residue* pSmallMol=NULL;
        StructureFindSmallMol(pStructure,&pSmallMol);
        if(AtomArrayCalcMinDistance(&pResi->atoms,&pSmallMol->atoms)<ENERGY_DISTANCE_CUTOFF){
          ProteinSiteBuildWildtypeRotamersByBBdepRotLib(pStructure,i,j,pBBdepRotLib,atomParams,resiTopos);
          if(FLAG_ADD_CRYSTAL_ROT) ProteinSiteAddCrystalRotamerByBBdepRotLib(pStructure,i,j,resiTopos,pBBdepRotLib);
          if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
          ProteinSiteCalcRotamersEnergy(pStructure,pAAppTable,pRama,i,j,fp);
          ProteinSiteRemoveDesignSite(pStructure,i,j);
        }
      }
      else{
        ProteinSiteBuildWildtypeRotamersByBBdepRotLib(pStructure,i,j,pBBdepRotLib,atomParams,resiTopos);
        if(FLAG_ADD_CRYSTAL_ROT) ProteinSiteAddCrystalRotamerByBBdepRotLib(pStructure,i,j,resiTopos,pBBdepRotLib);
        if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
        ProteinSiteCalcRotamersEnergy(pStructure,pAAppTable,pRama,i,j,fp);
        ProteinSiteRemoveDesignSite(pStructure,i,j);
      }
    }
  }
  fclose(fp);

  return Success;
}



int EVOEF_CheckRotamerInBBindRotLib(Structure* pStructure,BBindRotamerLib* pBBindRotLib,ResiTopoSet* pTopos,double cutoff,char* pdbid){
  char FILE_TORSION[MAX_LENGTH_FILE_NAME+1];
  sprintf(FILE_TORSION,"%s_torsionrecover.txt",pdbid);
  FILE* pf=fopen(FILE_TORSION,"w");  for(int i=0;i<StructureGetChainCount(pStructure);i++){
    Chain* pChain=StructureGetChain(pStructure,i);
    for(int j=0;j<ChainGetResidueCount(pChain);j++){
      Residue* pResidue=ChainGetResidue(pChain,j);
      if(strcmp(ResidueGetName(pResidue),"ALA")!=0 && strcmp(ResidueGetName(pResidue),"GLY")!=0){
        BOOL result=ProteinSiteCheckCrystalRotamerInBBindRotLib(pStructure,i,j,pTopos,pBBindRotLib,cutoff);
        if(result==TRUE){
          fprintf(pf,"%s %s 1\n",ChainGetName(pChain),ResidueGetName(pResidue));
        }
        else{
          fprintf(pf,"%s %s 0\n",ChainGetName(pChain),ResidueGetName(pResidue));
        }
      }
      else{
        fprintf(pf,"%s %s 1\n",ChainGetName(pChain),ResidueGetName(pResidue));
      }
    }
  }
  fclose(pf);
  return Success;
}


int EVOEF_CheckRotamerInBBdepRotLib(Structure* pStructure,BBdepRotamerLib* pBBdepRotLib,ResiTopoSet* pTopos,double cutoff,char* pdbid){
  char FILE_TORSION[MAX_LENGTH_FILE_NAME+1];
  sprintf(FILE_TORSION,"%s_torsionrecover.txt",pdbid);
  FILE* pf=fopen(FILE_TORSION,"w");
  for(int i=0;i<StructureGetChainCount(pStructure);i++){
    Chain* pChain=StructureGetChain(pStructure,i);
    for(int j=0;j<ChainGetResidueCount(pChain);j++){
      Residue* pResidue=ChainGetResidue(pChain,j);
      if(strcmp(ResidueGetName(pResidue),"ALA")!=0 && strcmp(ResidueGetName(pResidue),"GLY")!=0){
        BOOL result=ProteinSiteCheckCrystalRotamerInBBdepRotLib(pStructure,i,j,pTopos,pBBdepRotLib,cutoff);
        if(result==TRUE){
          fprintf(pf,"%s %s 1\n",ChainGetName(pChain),ResidueGetName(pResidue));
        }
        else{
          fprintf(pf,"%s %s 0\n",ChainGetName(pChain),ResidueGetName(pResidue));
        }
      }
      else{
        fprintf(pf,"%s %s 1\n",ChainGetName(pChain),ResidueGetName(pResidue));
      }
    }
  }
  fclose(pf);
  return Success;
}


int EVOEF_GetResiMinRmsdRotFromLab(Structure* pStructure,char* pdbid){
  char FILE_RMSD[MAX_LENGTH_FILE_NAME+1];
  sprintf(FILE_RMSD,"%s_minrmsd.txt",pdbid);
  FILE* pf=fopen(FILE_RMSD,"w");
  for(int i=0;i<StructureGetChainCount(pStructure);i++){
    Chain* pChain=StructureGetChain(pStructure,i);
    for(int j=0;j<ChainGetResidueCount(pChain);j++){
      Residue* pResidue=ChainGetResidue(pChain,j);
      DesignSite* pSite=StructureFindDesignSite(pStructure,i,j);
      if(pSite!=NULL){
        double minRMSD=1e8;
        RotamerSet* pSet=DesignSiteGetRotamers(pSite);
        for(int k=0;k<RotamerSetGetCount(pSet);k++){
          Rotamer* pRotamer=RotamerSetGet(pSet,k);
          RotamerRestore(pRotamer,pSet);
          double rmsd = RotamerAndResidueSidechainRMSD(pRotamer,pResidue);
          if(RotamerIsSymmetricalCheck(pRotamer)==TRUE){
            Rotamer tempRot;
            RotamerCreate(&tempRot);
            SymmetricalRotamerGenerate(&tempRot,pRotamer);
            double rmsd2=RotamerAndResidueSidechainRMSD(&tempRot,pResidue);
            if(rmsd2<rmsd) rmsd=rmsd2;
            RotamerDestroy(&tempRot);
          }
          if(rmsd<minRMSD) minRMSD=rmsd;
          RotamerExtract(pRotamer);
        }
        fprintf(pf,"%s %s %f\n",ChainGetName(pChain),ResidueGetName(pResidue),minRMSD);
      }
    }
  }
  fclose(pf);
  return Success;
}

int EVOEF_CompareTwoStructureSidechainTorsionAndRmsd(Structure*pStructure,Structure* pStructure2,FILE* pTorsion,FILE* pRmsd){
  for(int i=0; i<StructureGetChainCount(pStructure);i++){
    Chain* pChain=StructureGetChain(pStructure,i);
    Chain* pChain2=StructureGetChain(pStructure2,i);
    for(int j=0;j<ChainGetResidueCount(pChain);j++){
      Residue* pResi=ChainGetResidue(pChain,j);
      Residue* pResi2=ChainGetResidue(pChain2,j);
      char res=ThreeLetterAAToOneLetterAA(ResidueGetName(pResi));
      char res2=ThreeLetterAAToOneLetterAA(ResidueGetName(pResi2));
      if(ChainGetType(pChain)==Type_Chain_Protein){
        if(res != 'A' && res != 'G'){
          double rmsd = ResidueAndResidueSidechainRMSD(pResi,pResi2);
          if(ResidueIsSymmetricalCheck(pResi)==TRUE){
            Residue tempRot;
            ResidueCreate(&tempRot);
            SymmetricalResidueGenerate(&tempRot,pResi);
            double rmsd2=ResidueAndResidueSidechainRMSD(&tempRot,pResi2);
            if(rmsd2<rmsd) rmsd=rmsd2;
            ResidueDestroy(&tempRot);
          }
          fprintf(pRmsd, "%s %s %f\n",ChainGetName(pChain),ResidueGetName(pResi),rmsd);
#if 0
          if(ResidueAndResidueAllTorsionsAreSimilar(pResi,pResi2,TORSION_DEVIATION_CUTOFF)==TRUE){
            fprintf(pTorsion, "%s %s 1\n",ChainGetName(pChain),ResidueGetName(pResi));
          }
          else{
            fprintf(pTorsion, "%s %s 0\n",ChainGetName(pChain),ResidueGetName(pResi));
          }
#endif
#if 1
          IntArray simArray;
          IntArrayCreate(&simArray,0);
          //BOOL sim=TRUE;
          fprintf(pTorsion,"%s %s",ChainGetName(pChain),ResidueGetName(pResi));
          if(ResidueAndResidueCheckTorsionSimilarity(pResi,pResi2,TORSION_DEVIATION_CUTOFF,&simArray)==TRUE){
            fprintf(pTorsion," 1");
          }
          else{
            fprintf(pTorsion," 0");
          }
          for(int k=0;k<IntArrayGetLength(&simArray);k++){
            if(IntArrayGet(&simArray,k)==1 /*&& sim==TRUE*/){
              fprintf(pTorsion," 1");
            }
            else if(IntArrayGet(&simArray,k)==2){
              fprintf(pTorsion," 0");
              //sim=FALSE;
            }
            /*else{
              fprintf(pTorsion," -");
            }*/
          }
          fprintf(pTorsion,"\n");
          IntArrayDestroy(&simArray);
#endif
        }
        else{
          double rmsd=0.0;
          fprintf(pRmsd, "%s %s %f\n",ChainGetName(pChain),ResidueGetName(pResi),rmsd);
          fprintf(pTorsion, "%s %s 1\n",ChainGetName(pChain),ResidueGetName(pResi));
        }
      }
    }
  }
  return Success;
}

int EVOEF_CheckClash0(Structure *pStructure, double clashRatio){
  printf("checking structure clashes with atom pairwise distance < %.2f*(ri+rj)\n",clashRatio);
  for(int i=0;i<StructureGetChainCount(pStructure);i++){
    Chain *pChainI = StructureGetChain(pStructure, i);

    for(int j=0;j<ChainGetResidueCount(pChainI);j++){
      Residue* pResiIJ=ChainGetResidue(pChainI,j);
      for(int j2=j+1;j2<ChainGetResidueCount(pChainI);j2++){
        Residue* pResiIJ2=ChainGetResidue(pChainI,j2);
        for(int atom1=0;atom1<ResidueGetAtomCount(pResiIJ);atom1++){
          Atom* pAtom1=ResidueGetAtom(pResiIJ,atom1);
          if(pAtom1->isBBAtom==TRUE || AtomIsHydrogen(pAtom1)==TRUE) continue;
          for(int atom2=0;atom2<ResidueGetAtomCount(pResiIJ2);atom2++){
            Atom* pAtom2=ResidueGetAtom(pResiIJ2,atom2);
            if(pAtom2->isBBAtom==TRUE || AtomIsHydrogen(pAtom2)==TRUE) continue;
            double dist=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
            if(dist<clashRatio*(pAtom1->vdw_radius+pAtom2->vdw_radius)){
              printf("protein residue %s%d%s atom %s and protein residue %s%d%s atom %s has clash, dist=%f\n",
                AtomGetChainName(pAtom1),AtomGetPosInChain(pAtom1),ResidueGetName(pResiIJ),AtomGetName(pAtom1),
                AtomGetChainName(pAtom2),AtomGetPosInChain(pAtom2),ResidueGetName(pResiIJ2),AtomGetName(pAtom2),
                dist);
            }
          }
        }
      }
    }

    for(int k=i+1;k<StructureGetChainCount(pStructure);k++){
      Chain* pChainK=StructureGetChain(pStructure,k);
      for(int j=0;j<ChainGetResidueCount(pChainI);j++){
        Residue* pResiIJ=ChainGetResidue(pChainI,j);
        for(int l=0;l<ChainGetResidueCount(pChainK);l++){
          Residue* pResiKL=ChainGetResidue(pChainK,l);
          for(int atom1=0;atom1<ResidueGetAtomCount(pResiIJ);atom1++){
            Atom* pAtom1=ResidueGetAtom(pResiIJ,atom1);
            if(pAtom1->isBBAtom==TRUE || AtomIsHydrogen(pAtom1)==TRUE) continue;
            for(int atom2=0;atom2<ResidueGetAtomCount(pResiKL);atom2++){
              Atom* pAtom2=ResidueGetAtom(pResiKL,atom2);
              if(pAtom2->isBBAtom==TRUE || AtomIsHydrogen(pAtom2)==TRUE) continue;
              double dist=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
              if(dist<clashRatio*(pAtom1->vdw_radius+pAtom2->vdw_radius)){
                if(ChainGetType(pChainI)==Type_Chain_SmallMol){
                  printf("ligand residue %s%d%s atom %s and protein residue %s%d%s atom %s has clash, dist=%f\n",
                    AtomGetChainName(pAtom1),AtomGetPosInChain(pAtom1),ResidueGetName(pResiIJ),AtomGetName(pAtom1),
                    AtomGetChainName(pAtom2),AtomGetPosInChain(pAtom2),ResidueGetName(pResiKL),AtomGetName(pAtom2),
                    dist);
                }
                else if(ChainGetType(pChainK)==Type_Chain_SmallMol){
                  printf("protein residue %s%d%s atom %s and ligand residue %s%d%s atom %s has clash, dist=%f\n",
                    AtomGetChainName(pAtom1),AtomGetPosInChain(pAtom1),ResidueGetName(pResiIJ),AtomGetName(pAtom1),
                    AtomGetChainName(pAtom2),AtomGetPosInChain(pAtom2),ResidueGetName(pResiKL),AtomGetName(pAtom2),
                    dist);
                }
                else{
                  printf("protein residue %s%d%s atom %s and protein residue %s%d%s atom %s has clash, dist=%f\n",
                    AtomGetChainName(pAtom1),AtomGetPosInChain(pAtom1),ResidueGetName(pResiIJ),AtomGetName(pAtom1),
                    AtomGetChainName(pAtom2),AtomGetPosInChain(pAtom2),ResidueGetName(pResiKL),AtomGetName(pAtom2),
                    dist);
                }
              }
            }
          }
        }
      }
    }
  }

  return Success;
}

int EVOEF_CheckClash1(Structure *pStructure, double distThreshold){
  printf("checking structure clashes with atom pairwise distance < ri+rj-%.2f\n",distThreshold);
  for(int i=0;i<StructureGetChainCount(pStructure);i++){
    Chain *pChainI = StructureGetChain(pStructure, i);

    for(int j=0;j<ChainGetResidueCount(pChainI);j++){
      Residue* pResiIJ=ChainGetResidue(pChainI,j);
      for(int j2=j+1;j2<ChainGetResidueCount(pChainI);j2++){
        Residue* pResiIJ2=ChainGetResidue(pChainI,j2);
        for(int atom1=0;atom1<ResidueGetAtomCount(pResiIJ);atom1++){
          Atom* pAtom1=ResidueGetAtom(pResiIJ,atom1);
          if(pAtom1->isBBAtom==TRUE || AtomIsHydrogen(pAtom1)==TRUE) continue;
          for(int atom2=0;atom2<ResidueGetAtomCount(pResiIJ2);atom2++){
            Atom* pAtom2=ResidueGetAtom(pResiIJ2,atom2);
            if(pAtom2->isBBAtom==TRUE || AtomIsHydrogen(pAtom2)==TRUE) continue;
            double dist=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
            if(pAtom1->vdw_radius+pAtom2->vdw_radius-distThreshold>dist){
              printf("protein residue %s%d%s atom %s and protein residue %s%d%s atom %s has clash, dist=%f\n",
                AtomGetChainName(pAtom1),AtomGetPosInChain(pAtom1),ResidueGetName(pResiIJ),AtomGetName(pAtom1),
                AtomGetChainName(pAtom2),AtomGetPosInChain(pAtom2),ResidueGetName(pResiIJ2),AtomGetName(pAtom2),
                dist);
            }
          }
        }
      }
    }

    for(int k=i+1;k<StructureGetChainCount(pStructure);k++){
      Chain* pChainK=StructureGetChain(pStructure,k);
      for(int j=0;j<ChainGetResidueCount(pChainI);j++){
        Residue* pResiIJ=ChainGetResidue(pChainI,j);
        for(int l=0;l<ChainGetResidueCount(pChainK);l++){
          Residue* pResiKL=ChainGetResidue(pChainK,l);
          for(int atom1=0;atom1<ResidueGetAtomCount(pResiIJ);atom1++){
            Atom* pAtom1=ResidueGetAtom(pResiIJ,atom1);
            if(pAtom1->isBBAtom==TRUE || AtomIsHydrogen(pAtom1)==TRUE) continue;
            for(int atom2=0;atom2<ResidueGetAtomCount(pResiKL);atom2++){
              Atom* pAtom2=ResidueGetAtom(pResiKL,atom2);
              if(pAtom2->isBBAtom==TRUE || AtomIsHydrogen(pAtom2)==TRUE) continue;
              double dist=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
              if(pAtom1->vdw_radius+pAtom2->vdw_radius-distThreshold>dist){
                if(ChainGetType(pChainI)==Type_Chain_SmallMol){
                  printf("ligand residue %s%d%s atom %s and protein residue %s%d%s atom %s has clash, dist=%f\n",
                    AtomGetChainName(pAtom1),AtomGetPosInChain(pAtom1),ResidueGetName(pResiIJ),AtomGetName(pAtom1),
                    AtomGetChainName(pAtom2),AtomGetPosInChain(pAtom2),ResidueGetName(pResiKL),AtomGetName(pAtom2),
                    dist);
                }
                else if(ChainGetType(pChainK)==Type_Chain_SmallMol){
                  printf("protein residue %s%d%s atom %s and ligand residue %s%d%s atom %s has clash, dist=%f\n",
                    AtomGetChainName(pAtom1),AtomGetPosInChain(pAtom1),ResidueGetName(pResiIJ),AtomGetName(pAtom1),
                    AtomGetChainName(pAtom2),AtomGetPosInChain(pAtom2),ResidueGetName(pResiKL),AtomGetName(pAtom2),
                    dist);
                }
                else{
                  printf("protein residue %s%d%s atom %s and protein residue %s%d%s atom %s has clash, dist=%f\n",
                    AtomGetChainName(pAtom1),AtomGetPosInChain(pAtom1),ResidueGetName(pResiIJ),AtomGetName(pAtom1),
                    AtomGetChainName(pAtom2),AtomGetPosInChain(pAtom2),ResidueGetName(pResiKL),AtomGetName(pAtom2),
                    dist);
                }
              }
            }
          }
        }
      }
    }
  }

  return Success;
}


int EVOEF_CheckClash2(Structure *pStructure, double clashRatio){
  printf("checking structure clashes with atom pairwise distance < %.2f*(ri+rj)\n",clashRatio);
  for(int i=0;i<StructureGetChainCount(pStructure);i++){
    Chain *pChainI = StructureGetChain(pStructure, i);

    for(int j=0;j<ChainGetResidueCount(pChainI);j++){
      Residue* pResiIJ=ChainGetResidue(pChainI,j);
      for(int j2=j+1;j2<ChainGetResidueCount(pChainI);j2++){
        Residue* pResiIJ2=ChainGetResidue(pChainI,j2);
        for(int atom1=0;atom1<ResidueGetAtomCount(pResiIJ);atom1++){
          Atom* pAtom1=ResidueGetAtom(pResiIJ,atom1);
          if(pAtom1->isBBAtom==TRUE || AtomIsHydrogen(pAtom1)==TRUE) continue;
          for(int atom2=0;atom2<ResidueGetAtomCount(pResiIJ2);atom2++){
            Atom* pAtom2=ResidueGetAtom(pResiIJ2,atom2);
            if(pAtom2->isBBAtom==TRUE || AtomIsHydrogen(pAtom2)==TRUE) continue;
            double dist=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
            if(dist<clashRatio*(pAtom1->vdw_radius+pAtom2->vdw_radius)){
              printf("protein residue %s%d%s atom %s and protein residue %s%d%s atom %s has clash, dist=%f\n",
                AtomGetChainName(pAtom1),AtomGetPosInChain(pAtom1),ResidueGetName(pResiIJ),AtomGetName(pAtom1),
                AtomGetChainName(pAtom2),AtomGetPosInChain(pAtom2),ResidueGetName(pResiIJ2),AtomGetName(pAtom2),
                dist);
            }
          }
        }
      }
    }

    for(int k=i+1;k<StructureGetChainCount(pStructure);k++){
      Chain* pChainK=StructureGetChain(pStructure,k);
      for(int j=0;j<ChainGetResidueCount(pChainI);j++){
        Residue* pResiIJ=ChainGetResidue(pChainI,j);
        for(int l=0;l<ChainGetResidueCount(pChainK);l++){
          Residue* pResiKL=ChainGetResidue(pChainK,l);
          for(int atom1=0;atom1<ResidueGetAtomCount(pResiIJ);atom1++){
            Atom* pAtom1=ResidueGetAtom(pResiIJ,atom1);
            if(pAtom1->isBBAtom==TRUE || AtomIsHydrogen(pAtom1)==TRUE) continue;
            for(int atom2=0;atom2<ResidueGetAtomCount(pResiKL);atom2++){
              Atom* pAtom2=ResidueGetAtom(pResiKL,atom2);
              if(pAtom2->isBBAtom==TRUE || AtomIsHydrogen(pAtom2)==TRUE) continue;
              double dist=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
              if(dist<clashRatio*(pAtom1->vdw_radius+pAtom2->vdw_radius)){
                if(ChainGetType(pChainI)==Type_Chain_SmallMol){
                  printf("ligand residue %s%d%s atom %s and protein residue %s%d%s atom %s has clash, dist=%f\n",
                    AtomGetChainName(pAtom1),AtomGetPosInChain(pAtom1),ResidueGetName(pResiIJ),AtomGetName(pAtom1),
                    AtomGetChainName(pAtom2),AtomGetPosInChain(pAtom2),ResidueGetName(pResiKL),AtomGetName(pAtom2),
                    dist);
                }
                else if(ChainGetType(pChainK)==Type_Chain_SmallMol){
                  printf("protein residue %s%d%s atom %s and ligand residue %s%d%s atom %s has clash, dist=%f\n",
                    AtomGetChainName(pAtom1),AtomGetPosInChain(pAtom1),ResidueGetName(pResiIJ),AtomGetName(pAtom1),
                    AtomGetChainName(pAtom2),AtomGetPosInChain(pAtom2),ResidueGetName(pResiKL),AtomGetName(pAtom2),
                    dist);
                }
                else{
                  printf("protein residue %s%d%s atom %s and protein residue %s%d%s atom %s has clash, dist=%f\n",
                    AtomGetChainName(pAtom1),AtomGetPosInChain(pAtom1),ResidueGetName(pResiIJ),AtomGetName(pAtom1),
                    AtomGetChainName(pAtom2),AtomGetPosInChain(pAtom2),ResidueGetName(pResiKL),AtomGetName(pAtom2),
                    dist);
                }
              }
            }
          }
        }
      }
    }
  }

  return Success;
}



int EVOEF_ComputeStabilityWithBBdepRotLibNew(Structure *pStructure,AAppTable* pAAppTable,RamaTable* pRama,char* dunlibfile,double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]){
  /***********************************************************************************/
  /* calculate the dunbrack energy for the given protein
  /***********************************************************************************/
  /*ACDEFGHIKLMNPQRSTVWY, only for regular amino acid*/
  int xcount[20]={0,1,2,3,2,0,2,2,4,2,3,2,2,3,4,1,1,1,2,2};
  int nrot[20]={0,3,18,54,18,0,36,9,73,9,27,36,2,108,75,3,3,3,36,18};
  int lrot[20]={0,129,111,240,448,0,294,330,348,339,421,75,466,132,0,468,471,528,474,510};
  FILE* FIN=fopen(dunlibfile,"rb");
  for(int i = 0; i < StructureGetChainCount(pStructure); i++){
    Chain *pChainI = StructureGetChain(pStructure,i);
    if(ChainGetType(pChainI)==Type_Chain_Protein){
      for(int ir=0;ir<ChainGetResidueCount(pChainI);ir++){
        Residue* pResiIR=ChainGetResidue(pChainI,ir);
        int aaIndex=ThreeLetterAAGetIndex(ResidueGetName(pResiIR));
        if(aaIndex==0 || aaIndex==5) continue;
        int nchi=xcount[aaIndex];
        //int phiIdx=(int)floor((pResiIR->phipsi[0]+180.)/10.0);
        //if(phiIdx>=36) phiIdx-=36;
        //int psiIdx=(int)floor((pResiIR->phipsi[1]+180.)/10.0);
        //if(psiIdx>=36) psiIdx-=36;
        //int binindex=36*phiIdx+psiIdx;
        int binindex=((int)(pResiIR->phipsi[0]+180)/10)*36+(int)(pResiIR->phipsi[1]+180)/10;
        //printf("binindex = %d\n",binindex);
        //the last '36' means 36 bytes
        fseek(FIN,(1296*lrot[aaIndex]+binindex*nrot[aaIndex])*36,SEEK_SET);
        int matchIdx=-1;
        double pmin=10;
        double pmatch;
        for(int rotIdx=0;rotIdx<nrot[aaIndex];rotIdx++){
          float p=0.;
          fread((char*)&p,sizeof(char),4,FIN);
          if(p<ROT_PROB_CUT_MIN) break;
          if(p<pmin){pmin=p;}
          float val;
          double libtorsions[4]={0};
          double libdeviations[4]={0};
          for(int k=0;k<nchi;k++){
            fread((char*)&val,sizeof(char),4,FIN);
            libtorsions[k]=(float)val;
          }
          fseek(FIN,(4-nchi)*4,SEEK_CUR);
          for(int k=0;k<nchi;k++){
            fread((char*)&val,sizeof(char),4,FIN);
            libdeviations[k]=(float)val;
          }
          fseek(FIN,(4-nchi)*4,SEEK_CUR);
          /******************************************************************
          /* check if the residue torsion matches with one rotamer torsion
          /*****************************************************************/
          BOOL match=TRUE;
          for(int torIdx=0;torIdx<DoubleArrayGetLength(&pResiIR->xtorsions);torIdx++){
            double min=DegToRad(libtorsions[torIdx])-DegToRad(libdeviations[torIdx]);
            double max=DegToRad(libtorsions[torIdx])+DegToRad(libdeviations[torIdx]);
            //use a strict criteria
            //double min=DegToRad(libtorsions[torIdx])-DegToRad(5.0);
            //double max=DegToRad(libtorsions[torIdx])+DegToRad(5.0);
            double torsion=DoubleArrayGet(&pResiIR->xtorsions,torIdx);
            double torsionm2pi=torsion-2*PI;
            double torsionp2pi=torsion+2*PI;
            double torsion2=torsion;
            if((strcmp(ResidueGetName(pResiIR),"PHE")==0 && torIdx==1)||
              (strcmp(ResidueGetName(pResiIR),"TYR")==0 && torIdx==1)||
              (strcmp(ResidueGetName(pResiIR),"ASP")==0 && torIdx==1)||
              strcmp(ResidueGetName(pResiIR),"GLU")==0 && torIdx==2){
                torsion2=torsion+PI;
                torsion2=torsion>0?torsion-PI:torsion2;
            }
            double torsion2m2pi=torsion2-2*PI;
            double torsion2p2pi=torsion2+2*PI;
            if(!(
              (torsion    <=max && torsion>=min) ||
              (torsionm2pi<=max && torsionm2pi>=min) ||
              (torsionp2pi<=max && torsionp2pi>=min) ||
              (torsion2    <=max && torsion2>=min) ||
              (torsion2m2pi<=max && torsion2m2pi>=min) ||
              (torsion2p2pi<=max && torsion2p2pi>=min)
              )){
                match=FALSE;
                break;
            }
          }
          if(match==TRUE){
            matchIdx=rotIdx;
            pmatch=p;
            break;
          }
        }
        double pdelta=1e-7;
        if(matchIdx!=-1) energyTerms[93]+=-log(pmatch+pdelta);
        else energyTerms[93]+=-log(pmin+pdelta);
      }
    }
  }
  fclose(FIN);

  /**************************************************************/
  /* calculate the physical energy terms
  /**************************************************************/
  for(int i = 0; i < StructureGetChainCount(pStructure); i++){
    Chain *pChainI = StructureGetChain(pStructure,i);
    for(int ir = 0; ir < ChainGetResidueCount(pChainI); ir++){
      Residue *pResIR = ChainGetResidue(pChainI,ir);
      //ResidueShowAtomParameter(pResIR);
      //ResidueShowBondInformation(pResIR);
      if(ChainGetType(pChainI)==Type_Chain_Protein){
        EVOEF_AminoAcidReferenceEnergy(pResIR->name, energyTerms);
        EVOEF_EnergyResidueIntraEnergy(pResIR,energyTerms);
        AminoAcidPropensityAndRamachandranEnergy(pResIR,pAAppTable,pRama);
        energyTerms[91]+=pResIR->aapp;
        energyTerms[92]+=pResIR->ramachandran;
      }
      for(int is=ir+1;is<ChainGetResidueCount(pChainI);is++){
        Residue *pResIS = ChainGetResidue(pChainI,is);
        if(ResidueGetPosInChain(pResIR)+1==ResidueGetPosInChain(pResIS)) EVOEF_EnergyResidueAndNextResidue(pResIR,pResIS,energyTerms);
        else EVOEF_EnergyResidueAndOtherResidueSameChain(pResIR,pResIS,energyTerms);
        //EVOEF_EnergyResidueAndResidueSameChain(pResIR,pResIS,energyTerms);
      }
      for(int k = i+1; k < StructureGetChainCount(pStructure); k++){
        Chain *pChainK = StructureGetChain(pStructure,k);
        for(int ks = 0; ks < ChainGetResidueCount(pChainK); ks++){
          Residue *pResKS = ChainGetResidue(pChainK,ks);
          if(ChainGetType(pChainI)==Type_Chain_SmallMol){
            EVOEF_EnergyResidueAndLigandResidue(pResKS,pResIR,energyTerms);
          }
          else if(ChainGetType(pChainK)==Type_Chain_SmallMol){
            EVOEF_EnergyResidueAndLigandResidue(pResIR,pResKS,energyTerms);
          }
          else{
            EVOEF_EnergyResidueAndOtherResidueDifferentChain(pResIR,pResKS,energyTerms);
          }
        }
      }
    }
  }

  EnergyTermWeighting(energyTerms);
  printf("\nStructure energy details:\n");
  printf("reference_ALA         =            %12.6f\n", energyTerms[ 1]);
  printf("reference_CYS         =            %12.6f\n", energyTerms[ 2]);
  printf("reference_ASP         =            %12.6f\n", energyTerms[ 3]);
  printf("reference_GLU         =            %12.6f\n", energyTerms[ 4]);
  printf("reference_PHE         =            %12.6f\n", energyTerms[ 5]);
  printf("reference_GLY         =            %12.6f\n", energyTerms[ 6]);
  printf("reference_HIS         =            %12.6f\n", energyTerms[ 7]);
  printf("reference_ILE         =            %12.6f\n", energyTerms[ 8]);
  printf("reference_LYS         =            %12.6f\n", energyTerms[ 9]);
  printf("reference_LEU         =            %12.6f\n", energyTerms[10]);
  printf("reference_MET         =            %12.6f\n", energyTerms[11]);
  printf("reference_ASN         =            %12.6f\n", energyTerms[12]);
  printf("reference_PRO         =            %12.6f\n", energyTerms[13]);
  printf("reference_GLN         =            %12.6f\n", energyTerms[14]);
  printf("reference_ARG         =            %12.6f\n", energyTerms[15]);
  printf("reference_SER         =            %12.6f\n", energyTerms[16]);
  printf("reference_THR         =            %12.6f\n", energyTerms[17]);
  printf("reference_VAL         =            %12.6f\n", energyTerms[18]);
  printf("reference_TRP         =            %12.6f\n", energyTerms[19]);
  printf("reference_TYR         =            %12.6f\n", energyTerms[20]);

  printf("intraR_vdwatt         =            %12.6f\n", energyTerms[21]);
  printf("intraR_vdwrep         =            %12.6f\n", energyTerms[22]);
  printf("intraR_electr         =            %12.6f\n", energyTerms[23]);
  printf("intraR_deslvP         =            %12.6f\n", energyTerms[24]);
  printf("intraR_deslvH         =            %12.6f\n", energyTerms[25]);
  printf("intraR_hbscbb_dis     =            %12.6f\n", energyTerms[26]);
  printf("intraR_hbscbb_the     =            %12.6f\n", energyTerms[27]);
  printf("intraR_hbscbb_phi     =            %12.6f\n", energyTerms[28]);
  printf("aapropensity          =            %12.6f\n", energyTerms[91]);
  printf("ramachandran          =            %12.6f\n", energyTerms[92]);
  printf("dunbrack              =            %12.6f\n", energyTerms[93]);

  printf("interS_vdwatt         =            %12.6f\n", energyTerms[31]);
  printf("interS_vdwrep         =            %12.6f\n", energyTerms[32]);
  printf("interS_electr         =            %12.6f\n", energyTerms[33]);
  printf("interS_deslvP         =            %12.6f\n", energyTerms[34]);
  printf("interS_deslvH         =            %12.6f\n", energyTerms[35]);
  printf("interS_ssbond         =            %12.6f\n", energyTerms[36]);
  printf("interS_hbbbbb_dis     =            %12.6f\n", energyTerms[41]);
  printf("interS_hbbbbb_the     =            %12.6f\n", energyTerms[42]);
  printf("interS_hbbbbb_phi     =            %12.6f\n", energyTerms[43]);
  printf("interS_hbscbb_dis     =            %12.6f\n", energyTerms[44]);
  printf("interS_hbscbb_the     =            %12.6f\n", energyTerms[45]);
  printf("interS_hbscbb_phi     =            %12.6f\n", energyTerms[46]);
  printf("interS_hbscsc_dis     =            %12.6f\n", energyTerms[47]);
  printf("interS_hbscsc_the     =            %12.6f\n", energyTerms[48]);
  printf("interS_hbscsc_phi     =            %12.6f\n", energyTerms[49]);

  printf("interD_vdwatt         =            %12.6f\n", energyTerms[51]);
  printf("interD_vdwrep         =            %12.6f\n", energyTerms[52]);
  printf("interD_electr         =            %12.6f\n", energyTerms[53]);
  printf("interD_deslvP         =            %12.6f\n", energyTerms[54]);
  printf("interD_deslvH         =            %12.6f\n", energyTerms[55]);
  printf("interD_ssbond         =            %12.6f\n", energyTerms[56]);
  printf("interD_hbbbbb_dis     =            %12.6f\n", energyTerms[61]);
  printf("interD_hbbbbb_the     =            %12.6f\n", energyTerms[62]);
  printf("interD_hbbbbb_phi     =            %12.6f\n", energyTerms[63]);
  printf("interD_hbscbb_dis     =            %12.6f\n", energyTerms[64]);
  printf("interD_hbscbb_the     =            %12.6f\n", energyTerms[65]);
  printf("interD_hbscbb_phi     =            %12.6f\n", energyTerms[66]);
  printf("interD_hbscsc_dis     =            %12.6f\n", energyTerms[67]);
  printf("interD_hbscsc_the     =            %12.6f\n", energyTerms[68]);
  printf("interD_hbscsc_phi     =            %12.6f\n", energyTerms[69]);

  printf("prolig_vdwatt         =            %12.6f\n", energyTerms[71]);
  printf("prolig_vdwrep         =            %12.6f\n", energyTerms[72]);
  printf("prolig_electr         =            %12.6f\n", energyTerms[73]);
  printf("prolig_deslvP         =            %12.6f\n", energyTerms[74]);
  printf("prolig_deslvH         =            %12.6f\n", energyTerms[75]);
  printf("prolig_hbscbb_dis     =            %12.6f\n", energyTerms[81]);
  printf("prolig_hbscbb_the     =            %12.6f\n", energyTerms[82]);
  printf("prolig_hbscbb_phi     =            %12.6f\n", energyTerms[83]);
  printf("prolig_hbscsc_dis     =            %12.6f\n", energyTerms[84]);
  printf("prolig_hbscsc_the     =            %12.6f\n", energyTerms[85]);
  printf("prolig_hbscsc_phi     =            %12.6f\n", energyTerms[86]);
  printf("----------------------------------------------------\n");
  printf("Total                 =            %12.6f\n\n", energyTerms[0]);
  return Success;
}


