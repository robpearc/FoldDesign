///////////////////////////////////////////////////////////////////////////////////////////////
//         Copyright (C) Regents of the University of Michigan - All Rights Reserved         //
//          Unauthorized copying of this file, via any medium is strictly prohibited         //
//                                                                                           //
//                          Author: Robin Pearce <robpearc@umich.edu>                        //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "CommonPara.h"
#include "ParseInput.h"
#include "ParsePDB.h"
#include "PrintFunc.h"
#include "Fragments.h"
#include "ProgramFunctions.h"
#include "FoldDesignEnergyFunction.h"
#include "FoldDesignMovement.h"
#include "ParseGeometryFiles.h"
#include "ParseFoldDesignEnergyFiles.h"
#include "ParseFoldDesignMovementFiles.h"
#include "InputData.h"

void ProgramFunctions::FoldDesignScore(
  char *libDir,
  char *dataDir,
  char *inputFile,
  char *energyWeightFile,
  char *contactConstrFile,
  char *distConstrFile,
  char *pdbFile
)
{
  double energy;
  int i;
  int seqLength;

  //------ Object instantiation -------->
  ParseInput inputParser;
  ParsePDB pdbParser; 
  ParseGeometryFiles geometryFileParser;
  ParseFoldDesignEnergyFiles energyFileParser;

  InputData inputInfo;
  GeometryCalc geometry;
  BasicFunc bf;
  FoldDesignEnergyFunction energyFunction;

  //---------- Load input ------------>
  inputParser.configureInputFoldDesign(dataDir,inputFile);
  inputInfo.configureInput(inputParser);
  seqLength=inputInfo.getSeqLength();

  energyFileParser.loadFiles(libDir,dataDir,energyWeightFile,contactConstrFile,
                             distConstrFile,inputInfo);
  geometryFileParser.loadFiles(libDir);
  geometry.configureGeometry(geometryFileParser);
  energyFunction.setEnergyFunctionParameters(energyFileParser,geometry,inputInfo);

  //---------- Load input PDB ---------->
  vector<point3f> decstr_vec;
  pdbParser.loadPDBChain(pdbFile,decstr_vec);
   for(i=0;i<seqLength;i++)
   {
     decstr_vec[i].ss3=inputInfo.get3StateSSatPos(i);
     decstr_vec[i].ss8=inputInfo.get8StateSSatPos(i);
     decstr_vec[i].ssm=decstr_vec[i].ss3;
     decstr_vec[i].aaa=inputInfo.getResAtPos(i);
     decstr_vec[i].iaa=(unsigned char)(bf.aminoid(decstr_vec[i].aaa));
     if(decstr_vec[i].iaa>19) decstr_vec[i].iaa=5;
   }

   geometry.tor2strsg2(&decstr_vec[0],seqLength);
   geometry.str2tor(&decstr_vec[0],seqLength,3);

   //---------- Calc energy ------------->
   energy=energyFunction.calcTotalEnergy(&decstr_vec[0],seqLength);
   printf("Total Energy: %lf\n",energy);
}

void ProgramFunctions::generateFragments(
  char *libDir,
  char *dataDir,
  char *inputFile,
  char *templateList,
  bool zipFlag
)
{
  Fragments fragGenerator;
  ParseInput inputParser;
  InputData inputInfo;
  FILE *file;
  int i;
  int seqLength;
  int fragSizeForClustering=10;
  int maxFragSize=20;
  bool flagFileExists=false;
  char fragList[STD_FILE_NAME_LENGTH+1]; 
  char outputFile[STD_FILE_NAME_LENGTH+1];
  char seqFile[STD_FILE_NAME_LENGTH+1];
  char sseFile[STD_FILE_NAME_LENGTH+1];
  char mtxFile[STD_FILE_NAME_LENGTH+1];
  char featFile[STD_FILE_NAME_LENGTH+1];
  char cmd[STD_FILE_NAME_LENGTH+1];

  sprintf(fragList,"%s/%s",libDir,templateList);
  sprintf(mtxFile,"%s/fraglib_mtx",libDir);
  sprintf(featFile,"%s/fraglib_feat",libDir);

  inputParser.configureInputFragments(dataDir,inputFile);
  inputInfo.configureInput(inputParser); 
  seqLength=inputInfo.getSeqLength();
  for(i=1;i<=maxFragSize;i++)
  {
    sprintf(outputFile,"%s/%dseqfra.topse.bz2",dataDir,i);
    file=fopen(outputFile,"rb");
    if(file)
    {
      fclose(file);
      printf("%s exists. Skipping generation of fragments of size %d.\n",outputFile,i);
    }
    else
    {
      sprintf(outputFile,"%s/%dseqfra.topse",dataDir,i);
      file=fopen(outputFile,"rt");

      if(!file)
      {
        printf("Generating fragments of size %d.\n",i);
        fragGenerator.generateFragments(dataDir,fragList,inputInfo,mtxFile,featFile,
                                        outputFile,i);
        
        if(i==10)
        {
          sprintf(outputFile,"%s/topdh.topdh",dataDir);
          file=fopen(outputFile,"rt");
          if(file)
          {
            fclose(file);
            printf("%s exists. Skipping clustering of fragments.\n",outputFile);
          }
          else
          {
            fragGenerator.calcfragdh(dataDir,seqLength,outputFile,fragSizeForClustering,topno);
          }
        }
        
        //if(zipFlag)
        //{
        //  sprintf(cmd,"bzip2 %s/%dseqfra.topse",dataDir,i);
        //  system(cmd);
        //}
      }
      else
      {
        fclose(file);
        printf("%s exists. Skipping generation of fragments of size %d.\n",outputFile,i);
      }
    }
  }

  sprintf(outputFile,"%s/topdh.topdh",dataDir);
  file=fopen(outputFile,"rt");
  if(file)
  {
    fclose(file);
    printf("%s exists. Skipping clustering of fragments.\n",outputFile);
  }
  else
  {
    fragGenerator.calcfragdh(dataDir,seqLength,outputFile,fragSizeForClustering,topno);
  }

}

void ProgramFunctions::FoldDesignREMC(
  char *libDir,
  char *dataDir,
  char *inputFile,
  char *energyWeightFile,
  char *contactConstrFile,
  char *distConstrFile,
  char *fixedFileName,
  int cycles,
  int randomNum,
  bool fixFlag
)
{ 
  //---Object instantiation--->
  FoldDesignMovement *movement=new FoldDesignMovement[remcno];
  ParseInput inputParser;
  ParsePDB pdbParser; //QT
  ParseGeometryFiles geometryFileParser;
  ParseFoldDesignEnergyFiles energyFileParser;
  ParseFoldDesignMovementFiles movementFileParser;

  InputData inputInfo;
  FoldDesignEnergyFunction energyFunction;
  GeometryCalc geometry;
  BasicFunc bf;
  PrintFunc pf; 

  bool flagps,flagss,flagss2,flagfixpos,flagfixpdb;
  int i,j,k,l,m;
  int n_rep=40;
  char tmpname[300],tmfasta[300];
  clock_t remcstart, remcfinish;
  double duration;
  flagps=inputParser.configureInputFoldDesign(dataDir,inputFile);
  if(!flagps) return;
  
  inputInfo.configureInput(inputParser);
  int seqLength=inputInfo.getSeqLength();
  int sweepnum=int(30*sqrt(seqLength));
  int movetype=0;

  bool flagchange;
  int numpair;
  int numcycles=cycles;
  int temi,teme;
  int cyci=0;
  
  int nrand=randomNum; //random number
  int nthd=0;  //threading alignments
  double timecutoff=30;
  double trand,poschange,enechange;
  point3f *decexchange;
  
  // newly defined ---------->
  double T1=1;
  double T2=20;
  double S1=0.001;
  double S2=0.004;
  double d8,d10,dwell,da,db,dc,dd; // contact parameter
  double Ea[300];
  double E_min=10000;
  int N_Ea[300],N_swap_acc[300],N_swap_tot[300];
  int i_move_type;
  int N_acc[100],N_tot[100];
  FILE *file30;
  char outd[300];
  char move_name[100][10]; // array of string, two dimension
  double tmp;
  time_t rawtime;
  struct tm * timeinfo;
  
  //----- record starting time ----------->
  sprintf(outd,"%s/out.d",dataDir); // for output
  file30=fopen(outd,"wt");
  time (&rawtime);
  timeinfo = localtime (&rawtime);
  fprintf(file30,"simulation starting time: %s\n",asctime(timeinfo));
  printf("------> simulation starting time: %s\n",asctime(timeinfo));
  
  //------initial random number --------->
  sgenrand(nrand*1013+13);
  srand(nrand*1019+23);
  SelectStream(nrand);
  PutSeed(nrand*1021+43);
  double timethresh;
  timethresh=timecutoff*3600;
  
  double *betaseq=new double[n_rep];
  double *oldenergy=new double[n_rep];
  double *forenergy=new double[n_rep];//record before
  double *newenergy=new double[n_rep];
  bool flagmc;
  int *totatmpt=new int[n_rep];
  int *accatmpt=new int[n_rep];
  FILE *file;
  FILE *remcfile;
  point3f **mcdecstr;

  //------Load energy function, geometry, and movement files ------>
  energyFileParser.loadFiles(libDir,dataDir,energyWeightFile,contactConstrFile,
                             distConstrFile,inputInfo);
  geometryFileParser.loadFiles(libDir);
  geometry.configureGeometry(geometryFileParser); 
  energyFunction.setEnergyFunctionParameters(energyFileParser,geometry,inputInfo);
  movementFileParser.loadFiles(libDir,dataDir,inputInfo);

  fprintf(stdout,"Energy function set up\n");
 
  //------Set up simulation movement-------->
  for(j=0;j<n_rep;j++)
  {
    Ea[j]=0;
    N_Ea[j]=0;
    if(j==0)
    {
      //Configure first replica movement 
      movement[0].configureSimulation(libDir,dataDir,movementFileParser,energyFunction,
                                      inputInfo,geometry,T1,T2,S1,S2,n_rep);
    }
    else
    {
      //Overloaded assignment operator to set up all replicas
      movement[j]=movement[0];
    }
    betaseq[j]=1.0/movement[0].temparray[n_rep-1-j];
  }

  //------- for statistics of swap rate ----------->
  int *numexchange=new int[n_rep];
  for(j=0;j<n_rep;j++)
  {
    //numexchange[j]=0;
    N_swap_acc[j]=0;
    N_swap_tot[j]=0;
  }
  
  fprintf(file30, "seed of random number generator= %d\n",nrand);
  fprintf(file30, "Lch=%d, N_cycle=%d, N_rep=%d, n_sweep=%d\n",
          seqLength,numcycles,n_rep,sweepnum);
  fprintf(file30, "First_decoy_rec=%d Time_cut=%.3f T=[%.3f, %.3f]\n\n",
          cyci,timecutoff,movement[0].temparray[0],movement[0].temparray[n_rep-1]);
 
  //-------- generate initial conformation ---------------->
  fprintf(stdout,"Generating initial decoy conformations\n");
  mcdecstr=new point3f*[n_rep];
  for(m=0;m<n_rep;m++)
  {
    mcdecstr[m]=new point3f[seqLength];
    memset(mcdecstr[m],-1,seqLength*sizeof(point3f));
    movement[m].setInitDecoy(mcdecstr[m],seqLength);//phi psi
    for(j=0;j<seqLength;j++)
    {
      mcdecstr[m][j].ss3=inputInfo.get3StateSSatPos(j); 
      mcdecstr[m][j].ss8=inputInfo.get8StateSSatPos(j); 
      mcdecstr[m][j].aaa=inputInfo.getResAtPos(j);
      mcdecstr[m][j].iaa=(unsigned char)(bf.aminoid(mcdecstr[m][j].aaa));
      if(mcdecstr[m][j].iaa<0 || mcdecstr[m][j].iaa>19)
      {     
        mcdecstr[m][j].iaa=5;
      }
      mcdecstr[m][j].resind=j;
    }
    geometry.str2tor(mcdecstr[m],seqLength,3);


    sprintf(tmfasta,"%s/%s",dataDir,fixedFileName);
    FILE *fixFile;
    fixFile=fopen(tmfasta,"r");
    if(fixFile)
    {
      fclose(fixFile);
      vector<point3f> fixedStruct;
      pdbParser.loadPDBChain(tmfasta,fixedStruct);
      vector<int> fixed_pos;
      for(j=0;j<fixedStruct.size();j++)
      {
        k=fixedStruct[j].resind-1;
        mcdecstr[m][k].ptn.x=fixedStruct[j].ptn.x;
        mcdecstr[m][k].ptn.y=fixedStruct[j].ptn.y;
        mcdecstr[m][k].ptn.z=fixedStruct[j].ptn.z;
        mcdecstr[m][k].x=fixedStruct[j].x;
        mcdecstr[m][k].y=fixedStruct[j].y;
        mcdecstr[m][k].z=fixedStruct[j].z;
        mcdecstr[m][k].ptc.x=fixedStruct[j].ptc.x;
        mcdecstr[m][k].ptc.y=fixedStruct[j].ptc.y;
        mcdecstr[m][k].ptc.z=fixedStruct[j].ptc.z;
        mcdecstr[m][k].ssm=inputInfo.get3StateSSatPos(k); 
        strcpy(mcdecstr[m][k].name,"UNKON");
        mcdecstr[m][k].stype=mcdecstr[m][k].ssm;
        mcdecstr[m][k].aaa=inputInfo.getResAtPos(k);
        mcdecstr[m][k].residueid=inputInfo.getResAtPos(k);
        mcdecstr[m][k].iaa=(unsigned char)(bf.aminoid(inputInfo.getResAtPos(k)));
        if(mcdecstr[m][k].iaa<0 || mcdecstr[m][k].iaa>19)
        {
          mcdecstr[m][k].iaa=5;
        }
        fixed_pos.push_back(k);
      }
      movement[m].fixed_pos=fixed_pos;
      movement[m].fixed_flag=true;
      geometry.str2tor(mcdecstr[m],seqLength,3);
      geometry.tor2str(mcdecstr[m],seqLength,3);
      geometry.str2tor(mcdecstr[m],seqLength,3);
      movement[m].setInitDecoyFixed(mcdecstr[m],seqLength);
    }
    //movement[m].extractDistConstraints(mcdecstr[m],ps.seqnum,init-1,term); //TODO: fix start and end pos
 
    geometry.tor2strsg2(mcdecstr[m],seqLength);
    oldenergy[m]=energyFunction.calcTotalEnergy(mcdecstr[m],seqLength);
    inputParser.removeSingletonHelices(mcdecstr[m],seqLength);
  }

  for(m=0;m<n_rep;m++)
  {
    accatmpt[m]=0;
    totatmpt[m]=0;
    forenergy[m]=oldenergy[m];
  }
  
  //-------- initize trajectory files --------------->
  temi=n_rep-10;
  teme=n_rep-1;
  for(m=temi;m<=teme;m++)
  {
    sprintf(tmpname,"%s/remc%d.pdb",dataDir,m);
    file=fopen(tmpname,"wt");
    fclose(file);
  }
  for(i=0;i<=13;i++)
  {
    N_acc[i]=0;
    N_tot[i]=0;
  }

  //----------- replica-exchange Monte Carlo (REMC) starts --------------------->
  int num1=0;
  remcstart=clock(); 
  for(l=0;l<numcycles;l++) // n_cycle, llll
  {         
    remcfinish = clock();
    printf("doing i_cycle=%3d time(h)=%8.3f\n",l,(remcfinish-remcstart)/CLOCKS_PER_SEC/60/60.0);
    fprintf(file30,"doing i_cycle=%3d time(h)=%8.3f\n",l,(remcfinish-remcstart)/CLOCKS_PER_SEC/60/60.0);

    for(m=0;m<n_rep;m++)
    {// n_rep, mmmm
      movement[m].indCyc=l;
      for(j=0;j<sweepnum;j++)
      {
        totatmpt[m]++; // total number of moves
        //--------- attempt movement, run Metropolis here:
        movement[m].setHbondRamp(numcycles);
        flagmc=movement[m].attemptConformationalMove(mcdecstr[m],seqLength,
                                                     betaseq[m],oldenergy[m],
                                                     &newenergy[m],&movetype);
        num1++;
        i_move_type=movement[m].getMoveType();
        //printf("--- itype=%d\n",i_move_type);
        if(flagmc) //flagmc=true: accepted move; false: rejected
        {  
          accatmpt[m]++;  // count number of accepted moves for m-th replica
          oldenergy[m]=newenergy[m]; // record energy
          N_acc[i_move_type]++;
        }
        N_tot[i_move_type]++;
    
        //--------- check running time ------------------->
        remcfinish = clock(); 
        duration = (double)(remcfinish - remcstart)/CLOCKS_PER_SEC;
        if(duration>timethresh)
        {
          fprintf(file30,"\n premature stop at icycle=%d, i_rep=%d, i_move=%d, run_time(h)=%f time_cutoff(h)=%f\n",
                  l,m,j,duration/3600,timethresh/3600);
          break; //goto pos1;
        }

      } // jjjj: finish all monte carlo steps between
        // a pair of adjacents swaps for one replica
      if(duration>timethresh) break;
      Ea[m]+=oldenergy[m];
      N_Ea[m]++;
      if(oldenergy[m]<E_min) E_min=oldenergy[m];
      //printf("--- m=%d, Ea=%8.3f, N_Ea=%5d\n",m,Ea[m],N_Ea[m]);
      
      //--------output trajectory for SPICKER ---------------->
      if(l>=cyci && m>=temi)
      { //cyci=0, temi=n_rep-10 (top 10 replicas)
        geometry.tor2strsg2(mcdecstr[m],seqLength);
        sprintf(tmpname,"%s/remc%d.pdb",dataDir,m);
        pf.writetra(mcdecstr[m],seqLength,tmpname,l-cyci+1,oldenergy[m],1);
      }
    } // n_rep, mmmm: finish all replicas between a pair of adjacent swaps
    if(duration>timethresh) break;
    
    //--------- make swap between neighboring replicas ------------------------->
    if(l==l/2*2) //even
    { 
      numpair=n_rep/2;
      for(m=0;m<numpair;m++)
      {
        trand=Random();
        poschange=exp((betaseq[2*m+1]-betaseq[2*m])*(oldenergy[2*m+1]-oldenergy[2*m]));
        //printf("l=%5d, swap1,2=%5d,%5d\n",l,2*m,2*m+1);
        N_swap_tot[2*m]++;
        if(trand<poschange)
        {
          N_swap_acc[2*m]++;
          decexchange=mcdecstr[2*m];
          mcdecstr[2*m]=mcdecstr[2*m+1];
          mcdecstr[2*m+1]=decexchange;
          enechange=oldenergy[2*m];
          oldenergy[2*m]=oldenergy[2*m+1];
          oldenergy[2*m+1]=enechange;
        }
      }
    }
    else //odd
    {
      numpair=(n_rep-1)/2;
      for(m=0;m<numpair;m++)
      {
        trand=Random();
        poschange=exp((betaseq[2*m+2]-betaseq[2*m+1])*(oldenergy[2*m+2]-oldenergy[2*m+1]));
        //printf("l=%5d, swap1,2=%5d,%5d\n",l,2*m+1,2*m+2);
        N_swap_tot[2*m+1]++;
        if(trand<poschange)
        {
          N_swap_acc[2*m+1]++;
          decexchange=mcdecstr[2*m+1];
          mcdecstr[2*m+1]=mcdecstr[2*m+2];
          mcdecstr[2*m+2]=decexchange;
          enechange=oldenergy[2*m+1];
          oldenergy[2*m+1]=oldenergy[2*m+2];
          oldenergy[2*m+2]=enechange;
        }
      }
    }
    //^^^^^^^^^^^^ replica-swap is completed ^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
    for(m=0;m<n_rep;m++) forenergy[m]=oldenergy[m];
    
  }//n_cycle, llll: finish all swap cycles
  if(duration<=timethresh)
  {
    fprintf(file30,"\n normal stop at icycle=%d run_time(h)=%f time_cutoff(h)=%f\n",l,
            duration/3600,timethresh/3600);
  }

  //-------- output to out.d ------------->
  sprintf(move_name[0],"sub");
  sprintf(move_name[1],"phi");
  sprintf(move_name[2],"psi");
  sprintf(move_name[3],"top");
  sprintf(move_name[4],"rot");
  sprintf(move_name[5],"LMP");
  sprintf(move_name[6],"bbp");
  sprintf(move_name[7],"ome");
  sprintf(move_name[8],"len");
  sprintf(move_name[9],"ang");
  sprintf(move_name[10],"aaa");
  sprintf(move_name[11],"btn");
  sprintf(move_name[12],"sft");
  sprintf(move_name[13],"tra");
  fprintf(file30,"\n ----------- statics of individual movements -------------->\n");
  fprintf(file30," i_move  move_name  N_acc  N_tot N_acc/N_tot \n");
  for(i=0;i<=13;i++)
  {
    if(N_tot[i]>0)
      tmp=float(N_acc[i])/N_tot[i];
    else
      tmp=0;
    fprintf(file30,"%5d %5s %8d %8d %8.3f\n",i,move_name[i],N_acc[i],N_tot[i],tmp);
  }
  fprintf(file30,"\n ----------- summary of MC movements -------------->\n");
  fprintf(file30," i_rep     T(i) E_fin(i)  <E(i)> Na_swap/Nt_swap     N_a/Nt\n");
  for(m=0;m<n_rep;m++)
  {
    if(N_swap_tot[m]>0)
      tmp=N_swap_acc[m]/double(N_swap_tot[m]);
    else
      tmp=0;
    fprintf(file30,"%5d %8.3f %8.1f %8.1f %8.5f(%4d/%4d) %8.5f(%6d/%7d)\n",
                   m,1/betaseq[m],forenergy[m],Ea[m]/N_Ea[m],
                   tmp,N_swap_acc[m],N_swap_tot[m],
                   accatmpt[m]/double(totatmpt[m]),accatmpt[m],totatmpt[m]);
    
  }
  fprintf(file30,"E_min= %10.3f\n",E_min);
  
  //----record endding time ----------->
  time (&rawtime);
  timeinfo = localtime (&rawtime);
  fprintf(file30,"ending time: %s\n",asctime(timeinfo));
  fprintf(file30,"Total number of hours: %8.3f\n",duration/3600);
  printf("E_min= %10.3f\n",E_min);
  printf("ending time: %s\n",asctime(timeinfo));
  printf("Total number of hours: %8.3f\n",duration/3600);
  fclose(file30);
  //^^^^^^^^^^^^^ out.d is completed ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  //------------Free memory--------------->
  for(m=0;m<n_rep;m++) delete[]mcdecstr[m];
  delete[]mcdecstr;
  mcdecstr=NULL; 
  delete[]forenergy;
  forenergy=NULL; 
  delete[]oldenergy;
  oldenergy=NULL;
  delete[]newenergy;
  newenergy=NULL;
  delete[]betaseq;
  betaseq=NULL;
  delete[]totatmpt;
  totatmpt=NULL;
  delete[]accatmpt;
  accatmpt=NULL;
  delete[]movement;
  movement=NULL;
  delete[]numexchange;
  numexchange=NULL;

  return;
} // end of REMC_simulation
