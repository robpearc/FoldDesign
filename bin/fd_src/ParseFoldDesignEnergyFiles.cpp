///////////////////////////////////////////////////////////////////////////////////////////////
//         Copyright (C) Regents of the University of Michigan - All Rights Reserved         //
//          Unauthorized copying of this file, via any medium is strictly prohibited         //
//                                                                                           //
//                          Author: Robin Pearce <robpearc@umich.edu>                        //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "ParseFoldDesignEnergyFiles.h"

ParseFoldDesignEnergyFiles::ParseFoldDesignEnergyFiles()
{
  //---Ramachandran energy tables--->
  H_general_rama=NULL;
  E_general_rama=NULL;
  C_general_rama=NULL;
  H_specific_rama=NULL;
  E_specific_rama=NULL;
  C_specific_rama=NULL;
  G_specific_rama=NULL;
  I_specific_rama=NULL;
  B_specific_rama=NULL;
  S_specific_rama=NULL;
  T_specific_rama=NULL;

  //-----SSE packing energy tables---->
  energyHHPackAngle2=NULL;
  energyHHPackAngle3=NULL;
  energyHHPackAngle4=NULL;
  energyHHPackDist2=NULL;
  energyHHPackDist3=NULL;
  energySSPackAngle2=NULL;
  energySSPackAngle3=NULL;
  energySSPackDist1=NULL;
  energySSPackDist2=NULL;
  energySSPackDist3=NULL;
  energySSPackDist4=NULL;
  energyHSPackAngle2=NULL;
  energyHSPackAngle3=NULL;
  energyHSPackAngle4=NULL;
  energyHSPackDist2=NULL;
  energyHSPackDist3=NULL;
  energyHSPackDist4=NULL;

  kbp=NULL;
  countptr=NULL;
  solweight=NULL;
  contactConstr=NULL;
  distanceConstr=NULL;
  distanceWeight=NULL;
  distRestrType=NULL;
  paa=NULL;
}

ParseFoldDesignEnergyFiles::~ParseFoldDesignEnergyFiles()
{
  const int NUM_PHI_BINS=360;
  const int NUM_PHI_BINS_SSE=360;
  const int NUM_PSI_BINS_SSE=181;
  const int NUM_DIST_BINS_SSE=201;
  const int NUM_ATOMS=8;
  int i,j;

  if(H_general_rama)
  {
    for(i=0;i<NUM_PHI_BINS;i++)
    {
      delete[]H_general_rama[i];
    }
    delete[]H_general_rama;
    H_general_rama=NULL;
  }
  if(E_general_rama)
  {
    for(i=0;i<NUM_PHI_BINS;i++)
    {
      delete[]E_general_rama[i];
    }
    delete[]E_general_rama;
    E_general_rama=NULL;
  }
  if(C_general_rama)
  {
    for(i=0;i<NUM_PHI_BINS;i++)
    {
      delete[]C_general_rama[i];
    }
    delete[]C_general_rama;
    C_general_rama=NULL;
  }
  if(H_specific_rama)
  {
    for(i=0;i<NUM_PHI_BINS;i++)
    {
      delete[]H_specific_rama[i];
    }
    delete[]H_specific_rama;
    H_specific_rama=NULL;
  }
  if(E_specific_rama)
  {
    for(i=0;i<NUM_PHI_BINS;i++)
    {
      delete[]E_specific_rama[i];
    }
    delete[]E_specific_rama;
    E_specific_rama=NULL;
  }
  if(C_specific_rama)
  {
    for(i=0;i<NUM_PHI_BINS;i++)
    {
      delete[]C_specific_rama[i];
    }
    delete[]C_specific_rama;
    C_specific_rama=NULL;
  }
  if(G_specific_rama)
  {
    for(i=0;i<NUM_PHI_BINS;i++)
    {
      delete[]G_specific_rama[i];
    }
    delete[]G_specific_rama;
    G_specific_rama=NULL;
  }
  if(I_specific_rama)
  {
    for(i=0;i<NUM_PHI_BINS;i++)
    {
      delete[]I_specific_rama[i];
    }
    delete[]I_specific_rama;
    I_specific_rama=NULL;
  }
  if(B_specific_rama)
  {
    for(i=0;i<NUM_PHI_BINS;i++)
    {
      delete[]B_specific_rama[i];
    }
    delete[]B_specific_rama;
    B_specific_rama=NULL;
  }
  if(S_specific_rama)
  {
    for(i=0;i<NUM_PHI_BINS;i++)
    {
      delete[]S_specific_rama[i];
    }
    delete[]S_specific_rama;
    S_specific_rama=NULL;
  }
  if(T_specific_rama)
  {
    for(i=0;i<NUM_PHI_BINS;i++)
    {
      delete[]T_specific_rama[i];
    }
    delete[]T_specific_rama;
    T_specific_rama=NULL;
  }

  if(energyHSPackDist2)
  {
    for(i=0;i<NUM_DIST_BINS_SSE;i++)
    {
      delete[]energyHSPackDist2[i];
    }
    delete[]energyHSPackDist2;
    energyHSPackDist2=NULL;
  }

  if(energyHSPackDist3)
  {
    for(i=0;i<NUM_DIST_BINS_SSE;i++)
    {
      delete[]energyHSPackDist3[i];
    }
    delete[]energyHSPackDist3;
    energyHSPackDist3=NULL;
  }

  if(energyHSPackDist4)
  {
    for(i=0;i<NUM_DIST_BINS_SSE;i++)
    {
      delete[]energyHSPackDist4[i];
    }
    delete[]energyHSPackDist4;
    energyHSPackDist4=NULL;
  }

  if(energyHHPackDist2)
  {
    for(i=0;i<NUM_DIST_BINS_SSE;i++)
    {
      delete[]energyHHPackDist2[i];
    }
    delete[]energyHHPackDist2;
    energyHHPackDist2=NULL;
  }
  if(energyHHPackDist3)
  {
    for(i=0;i<NUM_DIST_BINS_SSE;i++)
    {
      delete[]energyHHPackDist3[i];
    }
    delete[]energyHHPackDist3;
    energyHHPackDist3=NULL;
  }


  if(energySSPackDist1)
  {
    for(i=0;i<NUM_DIST_BINS_SSE;i++)
    {
      delete[]energySSPackDist1[i];
    }
    delete[]energySSPackDist1;
    energySSPackDist1=NULL;
  }
  if(energySSPackDist2)
  {
    for(i=0;i<NUM_DIST_BINS_SSE;i++)
    {
      delete[]energySSPackDist2[i];
    }
    delete[]energySSPackDist2;
    energySSPackDist2=NULL;
  }
  if(energySSPackDist3)
  {
    for(i=0;i<NUM_DIST_BINS_SSE;i++)
    {
      delete[]energySSPackDist3[i];
    }
    delete[]energySSPackDist3;
    energySSPackDist3=NULL;
  }
  if(energySSPackDist4)
  {
    for(i=0;i<NUM_DIST_BINS_SSE;i++)
    {
      delete[]energySSPackDist4[i];
    }
    delete[]energySSPackDist4;
    energySSPackDist4=NULL;
  }
  if(energySSPackAngle2)
  {
    for(i=0;i<NUM_PHI_BINS_SSE;i++)
    {
      for(j=0;j<NUM_PSI_BINS_SSE;j++)
      {
        delete[]energySSPackAngle2[i][j];
      }
      delete[]energySSPackAngle2[i];
    }
    delete[]energySSPackAngle2;
    energySSPackAngle2=NULL;
  }
  if(energySSPackAngle3)
  {
    for(i=0;i<NUM_PHI_BINS_SSE;i++)
    {
      for(j=0;j<NUM_PSI_BINS_SSE;j++)
      {
        delete[]energySSPackAngle3[i][j];
      }
      delete[]energySSPackAngle3[i];
    }
    delete[]energySSPackAngle3;
    energySSPackAngle3=NULL;
  }

  if(energyHHPackAngle2)
  {
    for(i=0;i<NUM_PHI_BINS_SSE;i++)
    {
      for(j=0;j<NUM_PSI_BINS_SSE;j++)
      {
        delete[]energyHHPackAngle2[i][j];
      }
      delete[]energyHHPackAngle2[i];
    }
    delete[]energyHHPackAngle2;
    energyHHPackAngle2=NULL;
  }
  if(energyHHPackAngle3)
  {
    for(i=0;i<NUM_PHI_BINS_SSE;i++)
    {
      for(j=0;j<NUM_PSI_BINS_SSE;j++)
      {
        delete[]energyHHPackAngle3[i][j];
      }
      delete[]energyHHPackAngle3[i];
    }
    delete[]energyHHPackAngle3;
    energyHHPackAngle3=NULL;
  }
  if(energyHHPackAngle4)
  {
    for(i=0;i<NUM_PHI_BINS_SSE;i++)
    {
      for(j=0;j<NUM_PSI_BINS_SSE;j++)
      {
        delete[]energyHHPackAngle4[i][j];
      }
      delete[]energyHHPackAngle4[i];
    }
    delete[]energyHHPackAngle4;
    energyHHPackAngle4=NULL;
  }


  if(energyHSPackAngle2)
  {
    for(i=0;i<NUM_PHI_BINS_SSE;i++)
    {
      for(j=0;j<NUM_PSI_BINS_SSE;j++)
      {
        delete[]energyHSPackAngle2[i][j];
      }
      delete[]energyHSPackAngle2[i];
    }
    delete[]energyHSPackAngle2;
    energyHSPackAngle2=NULL;
  }
  if(energyHSPackAngle3)
  {
    for(i=0;i<NUM_PHI_BINS_SSE;i++)
    {
      for(j=0;j<NUM_PSI_BINS_SSE;j++)
      {
        delete[]energyHSPackAngle3[i][j];
      }
      delete[]energyHSPackAngle3[i];
    }
    delete[]energyHSPackAngle3;
    energyHSPackAngle3=NULL;
  }
  if(energyHSPackAngle4)
  {
    for(i=0;i<NUM_PHI_BINS_SSE;i++)
    {
      for(j=0;j<NUM_PSI_BINS_SSE;j++)
      {
        delete[]energyHSPackAngle4[i][j];
      }
      delete[]energyHSPackAngle4[i];
    }
    delete[]energyHSPackAngle4;
    energyHSPackAngle4=NULL;
  }

  if(countptr)
  {
    for(i=0;i<160;i++)
    {
      for(j=0;j<160;j++)
      {
        delete[]countptr[i][j];
      }
      delete[]countptr[i];
    }
    delete[]countptr;
  }
  if(kbp)
  {
    for(i=0;i<160;i++)
    {
      for(j=0;j<160;j++)
      {
        delete[]kbp[i][j];
        kbp[i][j]=NULL;
      }
      delete[]kbp[i];
      kbp[i]=NULL;
    }
    delete[]kbp;
    kbp=NULL;
  }

  if(solweight)
  {
    delete[]solweight;
    solweight=NULL;
  }
  if(contactConstr)
  {
    for(i=0;i<NUM_ATOMS;i++)
    {
      for(j=0;j<NUM_ATOMS;j++)
      {
        delete[]contactConstr[i][j];
      }
      delete[]contactConstr[i];
    }
    delete[]contactConstr;
    contactConstr=NULL;
  }
  if(distanceConstr)
  {
    for(i=0;i<NUM_ATOMS;i++)
    {
      for(j=0;j<NUM_ATOMS;j++)
      {
        delete[]distanceConstr[i][j];
      }
      delete[]distanceConstr[i];
    }
    delete[]distanceConstr;
    distanceConstr=NULL;
  }
  if(distanceWeight)
  {
    for(i=0;i<NUM_ATOMS;i++)
    {
      for(j=0;j<NUM_ATOMS;j++)
      {
        delete[]distanceWeight[i][j];
      }
      delete[]distanceWeight[i];
    }
    delete[]distanceWeight;
    distanceWeight=NULL;
  }
  if(distRestrType)
  {
    for(i=0;i<NUM_ATOMS;i++)
    {
      for(j=0;j<NUM_ATOMS;j++)
      {
        delete[]distRestrType[i][j];
      }
      delete[]distRestrType[i];
    }
    delete[]distRestrType;
    distRestrType=NULL;
  }
  if(paa)
  {
    delete[]paa;
    paa=NULL;
  }
}

bool ParseFoldDesignEnergyFiles::loadFiles(
  char *libDir,
  char *dataDir,
  char *energyWeightFile,
  char *contactConstrFile,
  char *distConstrFile,
  InputData inputInfo
)
{
  int seqLength=inputInfo.getSeqLength();
  char *sequence=inputInfo.getSequence();
  //sseInformation.configureSSEFoldDesign(dataDir);
  char fileName[STD_FILE_NAME_LENGTH+1]; 

  fprintf(stdout,"Loading data files for FoldDesign energy function.\n");

  //--- All energy weights --->
  loadEnergyWeights(libDir,energyWeightFile);
  
  //--- Ramachandran energy files --->
  loadAllRamachandranFiles(libDir);
  
  //--- Secondary structure packing files --->
  loadAllHelixHelixPackingFiles(libDir);
  loadAllStrandStrandPackingFiles(libDir);
  loadAllHelixStrandPackingFiles(libDir);
  
  //--- Solvation file --->
  sprintf(fileName,"%s/20_solv.txt",dataDir);
  loadSolSeq(fileName); 
  
  //--- Dfire --->
  sprintf(fileName,"%s/data.dat",libDir);
  loadDfire(fileName);
  
  //--- Read contact/distance constrints --->
  sprintf(fileName,"%s/%s",dataDir,contactConstrFile);
  loadContactRestr(fileName,seqLength); //read contact constraints
  sprintf(fileName,"%s/%s",dataDir,distConstrFile);
  loadDistanceRestr(fileName,seqLength); //read distance constraints
  getFragDistConstr(sequence,inputInfo);
 
  return true;
}

bool ParseFoldDesignEnergyFiles::loadEnergyWeights(
  char *libDir,
  char *weightFile
)
{
  for(int i=0;i<MAX_ENERGY_TERM_WEIGHTS;i++){
    weights[i]=1.0;
  }

  char fileName[STD_FILE_NAME_LENGTH+1];
  sprintf(fileName,"%s/%s",libDir,weightFile);

  FILE* file=fopen(fileName,"r");
  if(file)
  {
    char line[STD_ARRAY_SIZE];
    while(fgets(line,STD_ARRAY_SIZE,file))
    {
      char term[STD_ARRAY_SIZE];
      double val=0.0;
      sscanf(line,"%s %lf",term,&val);
      if(!strcmp(term,"ExcludedVolume"))                           weights[ 1]=val;
      else if(!strcmp(term,"HbondRamp1"))                          weights[ 2]=val;
      else if(!strcmp(term,"HbondRamp2"))                          weights[ 3]=val;
      else if(!strcmp(term,"UserContactConstr"))                   weights[ 4]=val;
      else if(!strcmp(term,"UserDistanceConstr"))                  weights[ 5]=val;
      else if(!strcmp(term,"Ramachandran"))                        weights[ 6]=val;
      else if(!strcmp(term,"RadiusGyration"))                      weights[ 7]=val;
      else if(!strcmp(term,"DistProfileContact"))                  weights[ 8]=val;
      else if(!strcmp(term,"BABmotifPenalty"))                     weights[ 9]=val;
      else if(!strcmp(term,"HelixHelixPackingPhiPsiTheta"))        weights[10]=val;
      else if(!strcmp(term,"HelixHelixPackingDistTheta"))          weights[11]=val;
      else if(!strcmp(term,"HelixHelixPackingPhiPsiThetaAlpha"))   weights[12]=val;
      else if(!strcmp(term,"HelixHelixPackingDistThetaAlpha"))     weights[13]=val;
      else if(!strcmp(term,"HelixStrandPackingPhiPsiTheta"))       weights[14]=val;
      else if(!strcmp(term,"HelixStrandPackingDistTheta"))         weights[15]=val;
      else if(!strcmp(term,"StrandStrandPackingPhiPsiTheta"))      weights[16]=val;
      else if(!strcmp(term,"StrandStrandPackingDistTheta"))        weights[17]=val;
      else if(!strcmp(term,"SecondaryStructureConstr"))            weights[18]=val;
      else if(!strcmp(term,"CA-CABondBreak"))                      weights[19]=val;
      else if(!strcmp(term,"Solvation"))                           weights[20]=val;
    }
    fclose(file);
  }
  else
  {
    fprintf(stderr,"Error in file %s function %s line %d, cannot open file %s\n",__FILE__,
            __FUNCTION__,__LINE__,weightFile);
    return true;
  }

  return true;
}

bool ParseFoldDesignEnergyFiles::loadAllRamachandranFiles(
  char *libDir
)
{
  char fileName[STD_FILE_NAME_LENGTH+1];
  bool flagSuccess; 
  int i;
  const int NUM_PHI_BINS=360,NUM_PSI_BINS=360;

  if(!H_general_rama)
  {
    H_general_rama=new double*[NUM_PHI_BINS];
    for(i=0;i<NUM_PHI_BINS;i++)
    {
      H_general_rama[i]=new double[NUM_PSI_BINS];
    }
  }
  if(!E_general_rama)
  {
    E_general_rama=new double*[NUM_PHI_BINS];
    for(i=0;i<NUM_PHI_BINS;i++)
    {
      E_general_rama[i]=new double[NUM_PSI_BINS];
    }
  }
  if(!C_general_rama)
  {
    C_general_rama=new double*[NUM_PHI_BINS];
    for(i=0;i<NUM_PHI_BINS;i++)
    {
      C_general_rama[i]=new double[NUM_PSI_BINS];
    }
  }
  if(!H_specific_rama)
  {
    H_specific_rama=new double*[NUM_PHI_BINS];
    for(i=0;i<NUM_PHI_BINS;i++)
    {
     H_specific_rama[i]=new double[NUM_PSI_BINS];
    }
  }
  if(!E_specific_rama)
  {
    E_specific_rama=new double*[NUM_PHI_BINS];
    for(i=0;i<NUM_PHI_BINS;i++)
    {
     E_specific_rama[i]=new double[NUM_PSI_BINS];
    }
  }
  if(!C_specific_rama)
  {
    C_specific_rama=new double*[NUM_PHI_BINS];
    for(i=0;i<NUM_PHI_BINS;i++)
    {
     C_specific_rama[i]=new double[NUM_PSI_BINS];
    }
  }
  if(!G_specific_rama)
  {
    G_specific_rama=new double*[NUM_PHI_BINS];
    for(i=0;i<NUM_PHI_BINS;i++)
    {
     G_specific_rama[i]=new double[NUM_PSI_BINS];
    }
  }
  if(!I_specific_rama)
  {
    I_specific_rama=new double*[NUM_PHI_BINS];
    for(i=0;i<NUM_PHI_BINS;i++)
    {
     I_specific_rama[i]=new double[NUM_PSI_BINS];
    }
  }
  if(!B_specific_rama)
  {
    B_specific_rama=new double*[NUM_PHI_BINS];
    for(i=0;i<NUM_PHI_BINS;i++)
    {
     B_specific_rama[i]=new double[NUM_PSI_BINS];
    }
  }
  if(!S_specific_rama)
  {
    S_specific_rama=new double*[NUM_PHI_BINS];
    for(i=0;i<NUM_PHI_BINS;i++)
    {
     S_specific_rama[i]=new double[NUM_PSI_BINS];
    }
  }
  if(!T_specific_rama)
  {
    T_specific_rama=new double*[NUM_PHI_BINS];
    for(i=0;i<NUM_PHI_BINS;i++)
    {
     T_specific_rama[i]=new double[NUM_PSI_BINS];
    }
  }

  sprintf(fileName,"%s/H_general_rama.txt",libDir);
  flagSuccess=loadRamaEnergyFile(fileName,H_general_rama);
  if(flagSuccess)
  {
    sprintf(fileName,"%s/E_general_rama.txt",libDir);
    flagSuccess=loadRamaEnergyFile(fileName,E_general_rama);
  }
  if(flagSuccess)
  {
    sprintf(fileName,"%s/C_general_rama.txt",libDir);
    flagSuccess=loadRamaEnergyFile(fileName,C_general_rama);
  }
  if(flagSuccess)
  {
    sprintf(fileName,"%s/H_specific_rama.txt",libDir);
    flagSuccess=loadRamaEnergyFile(fileName,H_specific_rama);
  }
  if(flagSuccess)
  {
    sprintf(fileName,"%s/E_specific_rama.txt",libDir);
    flagSuccess=loadRamaEnergyFile(fileName,E_specific_rama);
  }
  if(flagSuccess)
  {
    sprintf(fileName,"%s/C_specific_rama.txt",libDir);
    flagSuccess=loadRamaEnergyFile(fileName,C_specific_rama);
  }
  if(flagSuccess)
  {
    sprintf(fileName,"%s/G_specific_rama.txt",libDir);
    flagSuccess=loadRamaEnergyFile(fileName,G_specific_rama);
  }
  if(flagSuccess)
  {
    sprintf(fileName,"%s/I_specific_rama.txt",libDir);
    flagSuccess=loadRamaEnergyFile(fileName,I_specific_rama);
  }
  if(flagSuccess)
  {
    sprintf(fileName,"%s/B_specific_rama.txt",libDir);
    flagSuccess=loadRamaEnergyFile(fileName,B_specific_rama);
  }
  if(flagSuccess)
  {
    sprintf(fileName,"%s/S_specific_rama.txt",libDir);
    flagSuccess=loadRamaEnergyFile(fileName,S_specific_rama);
  }
  if(flagSuccess)
  {
    sprintf(fileName,"%s/T_specific_rama.txt",libDir);
    flagSuccess=loadRamaEnergyFile(fileName,T_specific_rama);
  }

  return flagSuccess;
}

bool ParseFoldDesignEnergyFiles::loadAllHelixHelixPackingFiles(
  char *libDir
)
{
  char fileName[STD_FILE_NAME_LENGTH+1];
  bool flagSuccess;
  const int NUM_PHI_BINS=360,NUM_PSI_BINS=181,NUM_THETA_BINS=181,NUM_DIST_BINS=201;
  int i,j;

  if(!energyHHPackAngle2)
  {
    energyHHPackAngle2=new double**[NUM_PHI_BINS];
    for(i=0;i<NUM_PHI_BINS;i++)
    {
      energyHHPackAngle2[i]=new double*[NUM_PSI_BINS];
      for(j=0;j<NUM_PSI_BINS;j++)
      {
        energyHHPackAngle2[i][j]=new double[NUM_THETA_BINS];
      }
    }
  }
  if(!energyHHPackAngle3)
  {
    energyHHPackAngle3=new double**[NUM_PHI_BINS];
    for(i=0;i<NUM_PHI_BINS;i++)
    {
      energyHHPackAngle3[i]=new double*[NUM_PSI_BINS];
      for(j=0;j<NUM_PSI_BINS;j++)
      {
        energyHHPackAngle3[i][j]=new double[NUM_THETA_BINS];
      }
    }
  }
  if(!energyHHPackAngle4)
  {
    energyHHPackAngle4=new double**[NUM_PHI_BINS];
    for(i=0;i<NUM_PHI_BINS;i++)
    {
      energyHHPackAngle4[i]=new double*[NUM_PSI_BINS];
      for(j=0;j<NUM_PSI_BINS;j++)
      {
        energyHHPackAngle4[i][j]=new double[NUM_THETA_BINS];
      }
    }
  }
  if(!energyHHPackDist2)
  {
    energyHHPackDist2=new double*[NUM_DIST_BINS];
    for(i=0;i<NUM_DIST_BINS;i++)
    {
      energyHHPackDist2[i]=new double[NUM_THETA_BINS];
    }
  }
  if(!energyHHPackDist3)
  {
    energyHHPackDist3=new double*[NUM_DIST_BINS];
    for(i=0;i<NUM_DIST_BINS;i++)
    {
      energyHHPackDist3[i]=new double[NUM_THETA_BINS];
    }
  }
 
  sprintf(fileName,"%s/energy_helix_helix_phithetasigma_seqsep2.txt",libDir);
  flagSuccess=loadSSEPackingFilePhiPsiTheta(fileName,energyHHPackAngle2);
  if(flagSuccess)
  {
    sprintf(fileName,"%s/energy_helix_helix_phithetasigma_seqsep3.txt",libDir);
    flagSuccess=loadSSEPackingFilePhiPsiTheta(fileName,energyHHPackAngle3);
  }
  if(flagSuccess)
  {
    sprintf(fileName,"%s/energy_helix_helix_phithetasigma_seqsep4.txt",libDir);
    flagSuccess=loadSSEPackingFilePhiPsiTheta(fileName,energyHHPackAngle4);
  }
  if(flagSuccess)
  {
    sprintf(fileName,"%s/hh_distsigma_seqsep2_energy.txt",libDir);
    flagSuccess=loadSSEPackingFileDistTheta(fileName,energyHHPackDist2);
  }
  if(flagSuccess)
  {
    sprintf(fileName,"%s/hh_distsigma_seqsep3_energy.txt",libDir);
    flagSuccess=loadSSEPackingFileDistTheta(fileName,energyHHPackDist3);
  }

  return flagSuccess;
}

bool ParseFoldDesignEnergyFiles::loadAllStrandStrandPackingFiles(
  char *libDir
)
{
  char fileName[STD_FILE_NAME_LENGTH+1];
  bool flagSuccess;
  const int NUM_PHI_BINS=360,NUM_PSI_BINS=181,NUM_THETA_BINS=181,NUM_DIST_BINS=201;
  int i,j;

  if(!energySSPackAngle2)
  {
    energySSPackAngle2=new double**[NUM_PHI_BINS];
    for(i=0;i<NUM_PHI_BINS;i++)
    {
      energySSPackAngle2[i]=new double*[NUM_PSI_BINS];
      for(j=0;j<NUM_PSI_BINS;j++)
      {
        energySSPackAngle2[i][j]=new double[NUM_THETA_BINS];
      }
    }
  }
  if(!energySSPackAngle3)
  {
    energySSPackAngle3=new double**[NUM_PHI_BINS];
    for(i=0;i<NUM_PHI_BINS;i++)
    {
      energySSPackAngle3[i]=new double*[NUM_PSI_BINS];
      for(j=0;j<NUM_PSI_BINS;j++)
      {
        energySSPackAngle3[i][j]=new double[NUM_THETA_BINS];
      }
    }
  }
  if(!energySSPackDist1)
  {
    energySSPackDist1=new double*[NUM_DIST_BINS];
    for(i=0;i<NUM_DIST_BINS;i++)
    {
      energySSPackDist1[i]=new double[NUM_THETA_BINS];
    }
  }
  if(!energySSPackDist2)
  {
    energySSPackDist2=new double*[NUM_DIST_BINS];
    for(i=0;i<NUM_DIST_BINS;i++)
    {
      energySSPackDist2[i]=new double[NUM_THETA_BINS];
    }
  }
  if(!energySSPackDist3)
  {
    energySSPackDist3=new double*[NUM_DIST_BINS];
    for(i=0;i<NUM_DIST_BINS;i++)
    {
      energySSPackDist3[i]=new double[NUM_THETA_BINS];
    }
  }
  if(!energySSPackDist4)
  {
    energySSPackDist4=new double*[NUM_DIST_BINS];
    for(i=0;i<NUM_DIST_BINS;i++)
    {
      energySSPackDist4[i]=new double[NUM_THETA_BINS];
    }
  }

  sprintf(fileName,"%s/energy_strand_strand_phithetasigma_seqsep2.txt",libDir); 
  flagSuccess=loadSSEPackingFilePhiPsiTheta(fileName,energySSPackAngle2);
  if(flagSuccess)
  {
    sprintf(fileName,"%s/energy_strand_strand_phithetasigma_seqsep3.txt",libDir);
    flagSuccess=loadSSEPackingFilePhiPsiTheta(fileName,energySSPackAngle3);
  }
  if(flagSuccess)
  {
    sprintf(fileName,"%s/ss_distsigma_dot1_energy.txt",libDir);
    flagSuccess=loadSSEPackingFileDistTheta(fileName,energySSPackDist1);
  }
  if(flagSuccess)
  {
    sprintf(fileName,"%s/ss_distsigma_dot2_energy.txt",libDir);
    flagSuccess=loadSSEPackingFileDistTheta(fileName,energySSPackDist2);
  }
  if(flagSuccess)
  {
    sprintf(fileName,"%s/ss_distsigma_dot3_energy.txt",libDir);
    flagSuccess=loadSSEPackingFileDistTheta(fileName,energySSPackDist3);
  }
  if(flagSuccess)
  {
    sprintf(fileName,"%s/ss_distsigma_dot4_energy.txt",libDir);
    flagSuccess=loadSSEPackingFileDistTheta(fileName,energySSPackDist4);
  }

  return flagSuccess;
}

bool ParseFoldDesignEnergyFiles::loadAllHelixStrandPackingFiles(
  char *libDir
)
{
  char fileName[STD_FILE_NAME_LENGTH+1];
  bool flagSuccess;
  const int NUM_PHI_BINS=360,NUM_PSI_BINS=181,NUM_THETA_BINS=181,NUM_DIST_BINS=201;
  int i,j;

  if(!energyHSPackAngle2)
  {
    energyHSPackAngle2=new double**[NUM_PHI_BINS];
    for(i=0;i<NUM_PHI_BINS;i++)
    {
      energyHSPackAngle2[i]=new double*[NUM_PSI_BINS];
      for(j=0;j<NUM_PSI_BINS;j++)
      {
        energyHSPackAngle2[i][j]=new double[NUM_THETA_BINS];
      }
    }
  }
  if(!energyHSPackAngle3)
  {
    energyHSPackAngle3=new double**[NUM_PHI_BINS];
    for(i=0;i<NUM_PHI_BINS;i++)
    {
      energyHSPackAngle3[i]=new double*[NUM_PSI_BINS];
      for(j=0;j<NUM_PSI_BINS;j++)
      {
        energyHSPackAngle3[i][j]=new double[NUM_THETA_BINS];
      }
    }
  }
  if(!energyHSPackAngle4)
  {
    energyHSPackAngle4=new double**[NUM_PHI_BINS];
    for(i=0;i<NUM_PHI_BINS;i++)
    {
      energyHSPackAngle4[i]=new double*[NUM_PSI_BINS];
      for(j=0;j<NUM_PSI_BINS;j++)
      {
        energyHSPackAngle4[i][j]=new double[NUM_THETA_BINS];
      }
    }
  }
  if(!energyHSPackDist2)
  {
    energyHSPackDist2=new double*[NUM_DIST_BINS];
    for(i=0;i<NUM_DIST_BINS;i++)
    {
      energyHSPackDist2[i]=new double[NUM_THETA_BINS];
    }
  }
  if(!energyHSPackDist3)
  {
    energyHSPackDist3=new double*[NUM_DIST_BINS];
    for(i=0;i<NUM_DIST_BINS;i++)
    {
      energyHSPackDist3[i]=new double[NUM_THETA_BINS];
    }
  }
  if(!energyHSPackDist4)
  {
    energyHSPackDist4=new double*[NUM_DIST_BINS];
    for(i=0;i<NUM_DIST_BINS;i++)
    {
      energyHSPackDist4[i]=new double[NUM_THETA_BINS];
    }
  }

  sprintf(fileName,"%s/energy_helix_strand_phithetasigma_seqsep2.txt",libDir);
  flagSuccess=loadSSEPackingFilePhiPsiTheta(fileName,energyHSPackAngle2);
  if(flagSuccess)
  {
    sprintf(fileName,"%s/energy_helix_strand_phithetasigma_seqsep3.txt",libDir);
    flagSuccess=loadSSEPackingFilePhiPsiTheta(fileName,energyHSPackAngle3);
  }
  if(flagSuccess)
  {
    sprintf(fileName,"%s/energy_helix_strand_phithetasigma_seqsep4.txt",libDir);
    flagSuccess=loadSSEPackingFilePhiPsiTheta(fileName,energyHSPackAngle4);
  }
  if(flagSuccess)
  {
    sprintf(fileName,"%s/hs_distsigma_seqsep2_energy.txt",libDir);
    flagSuccess=loadSSEPackingFileDistTheta(fileName,energyHSPackDist2);
  }
  if(flagSuccess)
  {
    sprintf(fileName,"%s/hs_distsigma_seqsep3_energy.txt",libDir);
    flagSuccess=loadSSEPackingFileDistTheta(fileName,energyHSPackDist3);
  }
  if(flagSuccess)
  {
    sprintf(fileName,"%s/hs_distsigma_seqsep4_energy.txt",libDir);
    flagSuccess=loadSSEPackingFileDistTheta(fileName,energyHSPackDist4);
  }

  return flagSuccess;
}

bool ParseFoldDesignEnergyFiles::loadRamaEnergyFile(
  char *ramaFile,
  double **ramaEnergyArray
)
{
  FILE *file;
  file=fopen(ramaFile,"rt");
  if(!file)
  {
    fprintf(stderr,"Error in file %s function %s line %d, cannot open file %s",__FILE__,
            __FUNCTION__,__LINE__,ramaFile);
    return false;
  }

  const int NUM_PHI_BINS=360,NUM_PSI_BINS=360;
  int i,j,ii,jj;
  double value;
  char line[STD_ARRAY_SIZE],label1[STD_ARRAY_SIZE],label2[STD_ARRAY_SIZE];

  for(i=0;i<NUM_PHI_BINS;i++)
  {
    for(j=0;j<NUM_PSI_BINS;j++)
    {
      fgets(line,STD_ARRAY_SIZE,file);
      sscanf(line,"%s %d %s %d %lf",&label1,&ii,&label2,&jj,&ramaEnergyArray[i][j]);
    }
  }
  fclose(file);
  
  return true;
}

bool ParseFoldDesignEnergyFiles::loadSSEPackingFilePhiPsiTheta(
  char *sseFile,
  double ***phiPsiThetaEnergyArray
)
{
  FILE *file;
  file=fopen(sseFile,"rt");
  if(!file)
  {
    fprintf(stderr,"Error in file %s function %s line %d, cannot open file %s\n",__FILE__,
            __FUNCTION__,__LINE__,sseFile);
    return false;
  }

  const int NUM_PHI_BINS=360,NUM_PSI_BINS=181,NUM_THETA_BINS=181;
  int i,j,k;
  char line[STD_ARRAY_SIZE];

  for(i=0;i<NUM_PHI_BINS;i++)
  {
    for(j=0;j<NUM_PSI_BINS;j++)
    {
      for(k=0;k<NUM_THETA_BINS;k++)
      {
        fgets(line,STD_ARRAY_SIZE,file);
        sscanf(line,"%lf",&phiPsiThetaEnergyArray[i][j][k]);
      }
    }
  }
  fclose(file);

  return true;
}

bool ParseFoldDesignEnergyFiles::loadSSEPackingFileDistTheta(
  char *sseFile,
  double **distThetaEnergyArray
)
{
  FILE *file;
  file=fopen(sseFile,"rt");
  if(!file)
  {
    fprintf(stderr,"Error in file %s function %s line %d, cannot open file %s\n",__FILE__,
            __FUNCTION__,__LINE__,sseFile);
    return false;
  }

  const int NUM_DIST_BINS=201,NUM_THETA_BINS=181;
  int i,j,k;
  int theta;
  double dist;
  char line[STD_ARRAY_SIZE],label1[STD_ARRAY_SIZE],label2[STD_ARRAY_SIZE];
  for(i=0;i<NUM_DIST_BINS;i++)
  {
    for(j=0;j<NUM_THETA_BINS;j++)
    {
      fgets(line,STD_ARRAY_SIZE,file);
      sscanf(line,"%s %lf %s %d %lf",&label1,&dist,&label2,&theta,&distThetaEnergyArray[i][j]);
    }
  }
  fclose(file);

  return true;
}

bool ParseFoldDesignEnergyFiles::loadSolSeq(
  char *solFile
)
{
  FILE *file;
  file=fopen(solFile,"rt");
  if(!file)
  {
    fprintf(stderr,"Error in file %s function %s line %d, cannot open file %s\n",__FILE__,
            __FUNCTION__,__LINE__,solFile);
    return false;
  }

  const double EXPOSED_WEIGHT=3.0,BURIED_WEIGHT=1.0;
  const double EXPOSED_THRESH=0.1;
  int i,j,solnum;
  double tval;
  char line[STD_ARRAY_SIZE];
  fgets(line,STD_ARRAY_SIZE,file);
  sscanf(line,"%d",&solnum);
  if(solweight) delete[]solweight;
  solweight=new double[solnum];
  for(i=0;i<solnum;i++)
  {
    fgets(line,STD_ARRAY_SIZE,file);
    sscanf(line,"%d %lf",&j,&tval);
    if(tval<=EXPOSED_THRESH)
    {
      solweight[i]=EXPOSED_WEIGHT;
    }
    else
    {
      solweight[i]=BURIED_WEIGHT;
    }
  }
  fclose(file);

  return true;
}

bool ParseFoldDesignEnergyFiles::loadDfire(
  char *filename
)
{
  FILE *file;
  char bf[2000];
  int i,j,k,num;

  if((file=fopen(filename,"rt"))==NULL)
  {
    printf("Unable to open dfire data %s\n",filename);
    return false;
  }
  if(!countptr)
  {
    countptr=new int**[160];
    for(i=0;i<160;i++)
    {
      countptr[i]=new int*[160];
      for(j=0;j<160;j++)
        countptr[i][j]=new int[30];
    }
  }
  for(k=0;k<30;k++)
  {
    fgets(bf,2000,file);
    for(j=1;j<159;j++)
    {
      fgets(bf,2000,file);
      for(i=1;i<159;i++)
      {
        sscanf(bf+(i-1)*10,"%10d",&num);
        countptr[i][j][k]=num;
      }
    }
    fgets(bf,2000,file);
  }
  fclose(file);
  if(!kbp)
  {
    kbp=new double**[160];
    for(i=0;i<160;i++)
    {
      kbp[i]=new double*[160];
      for(j=0;j<160;j++)
        kbp[i][j]=new double[30];
    }
  }
  for(k=0;k<30;k++)
  {
    for(j=1;j<159;j++)
    {
      for(i=j;i<159;i++)
      {
        if(countptr[i][j][k]==0 && countptr[j][i][k]==0)
        {
          kbp[i][j][k]=0.1;//origin james is 0.1, the same in dfire
        }
        else if(i==j)
        {
          kbp[i][j][k]=(0.016784)*(-0.001)*RR*TT*log((countptr[i][j][k])/pow(((k+0.5)/29.5),ALPHA)/(countptr[i][j][29]));
        }
        else
        {
          kbp[i][j][k]=(0.016784)*(-0.001)*RR*TT*log((countptr[i][j][k]+countptr[j][i][k])/pow(((k+0.5)/29.5),ALPHA)/(countptr[i][j][29]+countptr[j][i][29]));
        }
        kbp[j][i][k]=kbp[i][j][k];
      }
    }
  }
  return true;
}



bool ParseFoldDesignEnergyFiles::loadContactRestr(
  char *contactConstrFile,
  int seqLength
)
{
  const int NUM_ATOMS=8;
  FILE *file;
  int i,j,index,k=0;
  int istart,iend,atom_type1,atom_type2;
  float ialarm;
  char line[STD_ARRAY_SIZE];

  file=fopen(contactConstrFile,"rt");
  if(contactConstr==NULL)
  {
    contactConstr=new float**[NUM_ATOMS];
    for(i=0;i<NUM_ATOMS;i++)
    {
      contactConstr[i]=new float*[NUM_ATOMS];
    }
    for(i=0;i<NUM_ATOMS;i++)
    {
      for(j=0;j<NUM_ATOMS;j++)
      {
        contactConstr[i][j]=new float[seqLength*(seqLength+1)/2];
      }
    }
  }

  for(i=0;i<NUM_ATOMS;i++)
  {
    for(j=0;j<NUM_ATOMS;j++)
    {
      for(k=0;k<seqLength*(seqLength+1)/2;k++)
      {
        contactConstr[i][j][k]=0.0;
      }
    }
  }

  if(!file)
  {
    fprintf(stdout,"No contact restraint file %s\n",contactConstrFile);
    return false;
  }

  fgets(line,STD_ARRAY_SIZE,file); // read file head
  while(!feof(file))
  {
    fgets(line,STD_ARRAY_SIZE,file);
    sscanf(line,"%d %d %f %d %d",&istart,&iend,&ialarm,&atom_type1,&atom_type2);
    if(istart<=iend)
    {
      index=(istart-1)*seqLength-((istart-1)-1)*(istart-1)/2+(iend-1)-(istart-1);
    }
    else
    {
      index=(iend-1)*seqLength-((iend-1)-1)*(iend-1)/2+(istart-1)-(iend-1);
    }
    contactConstr[atom_type1][atom_type2][index]=ialarm;
  }
  fclose(file);

  return true;
}

bool ParseFoldDesignEnergyFiles::loadDistanceRestr(
  char *distanceFile,
  int seqnum
)
{
  const int NUM_ATOMS=8;
  FILE *file;
  int i,j,index,k=0;
  int istart,iend,atom_type1,atom_type2,restrType;
  float weight,dist;
  char line[STD_ARRAY_SIZE];
  file=fopen(distanceFile,"rt");

  if(distanceConstr==NULL)
  {
    distanceConstr=new float**[NUM_ATOMS];
    for(i=0;i<NUM_ATOMS;i++)
    {
      distanceConstr[i]=new float*[NUM_ATOMS];
    }
    for(i=0;i<NUM_ATOMS;i++)
    {
      for(j=0;j<NUM_ATOMS;j++)
      {
        distanceConstr[i][j]=new float[seqnum*(seqnum+1)/2];
      }
    }
  }

  if(!distanceWeight)
  {
    distanceWeight=new float**[NUM_ATOMS];
    for(i=0;i<NUM_ATOMS;i++)
    {
      distanceWeight[i]=new float*[NUM_ATOMS];
    }
    for(i=0;i<NUM_ATOMS;i++)
    {
      for(j=0;j<NUM_ATOMS;j++)
      {
        distanceWeight[i][j]=new float[seqnum*(seqnum+1)/2];
      }
    }
  }

  if(distRestrType==NULL)
  {
    distRestrType=new int**[NUM_ATOMS];
    for(i=0;i<NUM_ATOMS;i++)
    {
      distRestrType[i]=new int*[NUM_ATOMS];
    }
    for(i=0;i<NUM_ATOMS;i++)
    {
      for(j=0;j<NUM_ATOMS;j++)
      {
        distRestrType[i][j]=new int[seqnum*(seqnum+1)/2];
      }
    }
  }


  for(i=0;i<NUM_ATOMS;i++)
  {
    for(j=0;j<NUM_ATOMS;j++)
    {
      for(k=0;k<seqnum*(seqnum+1)/2;k++)
      {
        distanceConstr[i][j][k]=0.0;
        distanceWeight[i][j][k]=0.0;
        distRestrType[i][j][k]=0;
      }
    }
  }

  if(!file)
  {
    fprintf(stdout,"No distance restraint file %s\n",distanceFile);
    return false;
  }

  fgets(line,STD_ARRAY_SIZE,file); // read file head
  while(!feof(file))
  {
    fgets(line,STD_ARRAY_SIZE,file);
    sscanf(line,"%d %d %f %f %d %d %d",&istart,&iend,&weight,&dist,&atom_type1,&atom_type2,&restrType);
    if(istart<=iend)
    {
      index=(istart-1)*seqnum-((istart-1)-1)*(istart-1)/2+(iend-1)-(istart-1);
    }
    else
    {
      index=(iend-1)*seqnum-((iend-1)-1)*(iend-1)/2+(istart-1)-(iend-1);
    }
    distanceConstr[atom_type1][atom_type2][index]=dist;
    distanceWeight[atom_type1][atom_type2][index]=weight;
    distRestrType[atom_type1][atom_type2][index]=restrType;
  }
  fclose(file);

  return true;
}

void ParseFoldDesignEnergyFiles::getFragDistConstr(
  char *seqdat,
  InputData inputInfo
) // this one is used
{
  int i,j;
  int cutdel=1;
  int maxSep=3;
  int fragleng=9;//10 better in top5; 9 better in first  5 for betapro
  double wtleng=4.0;
  point3f **fragcont[20];
  int proleng=inputInfo.getSeqLength();

  for(i=0;i<20;i++)
  {
    fragcont[i]=new point3f*[(i+1)*topno];
    for(j=0;j<(i+1)*topno;j++) fragcont[i][j]=new point3f[proleng-(i+1)+1];
    fragcont[i]=inputInfo.getFragConformation(i);
  }

  if(!fragcont[fragleng-1]) return;
  if(paa) delete[]paa;
  int ii,jj,iii,jjj,aa,bb,cc,dd,aaa,bbb;
  int **paanum;
  double **paamean;
  double **paastd;
  int **distbins[20];
  int tbin;
  for(i=0;i<20;i++)
  {
    distbins[i]=new int*[proleng];
    for(j=0;j<proleng;j++)
    {
      distbins[i][j]=new int[proleng];
      for(ii=0;ii<proleng;ii++) distbins[i][j][ii]=0;
    }
  }

  paanum=new int*[proleng];
  paamean=new double*[proleng];
  paastd=new double*[proleng];
  for(i=0;i<proleng;i++)
  {
    paanum[i]=new int[proleng];
    paamean[i]=new double[proleng];
    paastd[i]=new double[proleng];
    for(j=0;j<proleng;j++)
    {
      paanum[i][j]=0;
      paamean[i][j]=0;
      paastd[i][j]=0;
    }
  }

  point3d tp[5];
  BasicFunc bf;
  double tdist;
  int inda,indb;
  for(i=0;i<proleng-fragleng+1;i++)
  {
    for(ii=0;ii<topno;ii++)
    {
      iii=fragleng*ii;
      for(aa=0;aa<fragleng;aa++)
      {
        aaa=i+aa;
        inda=bf.aminoid(seqdat[aaa]);
        if(inda>19 || inda<0) inda=5;
        tp[0]=bf.setv(fragcont[fragleng-1][iii+aa][i].x,
                      fragcont[fragleng-1][iii+aa][i].y,
                      fragcont[fragleng-1][iii+aa][i].z);
        for(bb=aa+cutdel;bb<fragleng-maxSep;bb++)
        {
          bbb=i+bb;
          indb=bf.aminoid(seqdat[bbb]);
          if(indb>19 || indb<0) indb=5;
          tp[1]=bf.setv(fragcont[fragleng-1][iii+bb][i].x,
                        fragcont[fragleng-1][iii+bb][i].y,
                        fragcont[fragleng-1][iii+bb][i].z);
          tp[2]=bf.minu(tp[1],tp[0]);
          tdist=bf.norm(tp[2]);
          paamean[aaa][bbb]+=tdist;
          paastd[aaa][bbb]+=tdist*tdist;
          paanum[aaa][bbb]++;
        }
      }
    }
  }

  paa=new pairaa[proleng*proleng];
  for(i=0;i<proleng*proleng;i++) paa[i].ispair=false;
  for(i=0;i<proleng;i++)
  {
    inda=bf.aminoid(seqdat[i]);
    if(inda>19 || inda<0) inda=5;
    for(j=i+1;j<proleng;j++)
    {
      indb=bf.aminoid(seqdat[j]);
      if(indb>19 || indb<0) indb=5;
      if(paanum[i][j]>0)
      {
        paamean[i][j]/=double(paanum[i][j]);
        paastd[i][j]=paastd[i][j]/double(paanum[i][j])-paamean[i][j]*paamean[i][j];
        if(paastd[i][j]<0.0) paastd[i][j]=0.0;
        paastd[i][j]=sqrt(paastd[i][j]);
        if(paanum[i][j]>wtleng*fragleng)// && distbins[18][i][j]<int(cutdist*2.0)-1)
        {
          ii=i*proleng+j;
          paa[ii].ispair=true;
          paa[ii].tnum=paanum[i][j];
          paa[ii].dist=paamean[i][j];
          paa[ii].dstd=paastd[i][j];
        }//if
      }
    }
  }
  for(i=0;i<proleng;i++)
  {
    delete[]paanum[i];
    delete[]paamean[i];
    delete[]paastd[i];
  }
  for(i=0;i<20;i++)
  {
    for(j=0;j<proleng;j++) delete[]distbins[i][j];
    delete[]distbins[i];
  }
}

