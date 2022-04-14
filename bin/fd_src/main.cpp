///////////////////////////////////////////////////////////////////////////////////////////////
//         Copyright (C) Regents of the University of Michigan - All Rights Reserved         //
//          Unauthorized copying of this file, via any medium is strictly prohibited         //
//                                                                                           //
//                          Author: Robin Pearce <robpearc@umich.edu>                        //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

 /********************************************************************************************
 * oooooooooooo      oooo       .o8 oooooooooo.                  o8o                         *
 * `888'     `8      `888      "888 `888'   `Y8b                 `"'                         *
 *  888      .ooooo.  888  .oooo888  888      888 .ooooo.  .oooo.ooooo  .ooooooooooo. .oo.   *
 *  888oooo d88' `88b 888 d88' `888  888      888d88' `88bd88(  "8`888 888' `88b `888P"Y88b  *
 *  888     888   888 888 888   888  888      888888ooo888`"Y88b.  888 888   888  888   888  *
 *  888     888   888 888 888   888  888     d88'888    .oo.  )88b 888 `88bod8P'  888   888  *
 * o888o    `Y8bod8P o888o`Y8bod88P"o888bood8P'  `Y8bod8P'8""888P'o888o`8oooooo. o888o o888o *
 *                                                                     d"     YD             *
 *                                                                     "Y88888P'             *
 * ******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ProgramFunctions.h"
#include "Getopt.h"

using namespace std;

bool CheckCommandName(char* queryname);
void FoldDesign_interface();
void FoldDesign_help();

const int MAX_LENGTH_ERR_MSG=1024;

const char* short_opts = "hv";
struct option long_opts[] = {
  {"help",          no_argument,       NULL,  1},
  {"version",       no_argument,       NULL,  2},
  {"command",       required_argument, NULL,  3},
  {"data_dir",      required_argument, NULL,  4},
  {"library",       required_argument, NULL,  5},
  {"fix_file",      required_argument, NULL,  6},
  {"random_num",    required_argument, NULL,  7}, 
  {"train_weight",  required_argument, NULL,  8},
  {"pdb",           required_argument, NULL,  9},
  {"input",         required_argument, NULL, 10},
  {"cont_restr_file",         required_argument, NULL, 11},
  {"dist_restr_file",         required_argument, NULL, 12},
  {"remc_cycles",         required_argument, NULL, 13},
  {NULL,            no_argument,       NULL, 0}
};

char *dataDir            = "./";
char *libDir             = "./";
char *weightFile         = "fold_design_energy_weights.txt";
char *contactConstrFile  = "./contact_restr.txt";
char *distanceConstrFile = "./distance_restr.txt";
char *templateList       = "new_cullpdb";
char *fixedFile          = "fixed_motif.pdb";
char *pdbFile            = "test.pdb";
char *inputFile          = "input.txt";
double cuttime           = -1;
double trainWeight       = 1.0;
int cycles               = 500;
int rannum               = 1032;
int randomNum            = 102;
const int success=0;
const int failure=1;
bool fixFlag=false;


int main(int argc, char** argv)
{
  ProgramFunctions func;
  FoldDesign_interface();
  char *cmdName = "FoldDesign";
  
  while(true)
  {
    int opt=getopt_long(argc, argv, short_opts, long_opts, NULL);
    if(opt == -1)
    {
      break;
    }
    switch(opt)
    {
      case 'h':
        FoldDesign_help();
        exit(success);
      case 'v':
        exit(success);
      case 1:
        FoldDesign_help();
        exit(success);
      case 2:
        exit(success);
      case 3:
        cmdName = optarg;
        if(!CheckCommandName(cmdName)){
          printf("Command %s is not supported! Exiting program.\n", cmdName);
          exit(failure);
        }
        break;
      case 4:
        dataDir=optarg;
        break;
      case 5:
        libDir=optarg;
        break;
      case 6:
	fixedFile=optarg;
        break;
      case 7:
        randomNum=atoi(optarg);
        break;
      case 8:
        trainWeight=atof(optarg);
        break;
      case 9:
        pdbFile=optarg;
        break;
      case 10:
        inputFile=optarg;
        break;
      case 11:
        contactConstrFile=optarg;
        break;
      case 12:
        distanceConstrFile=optarg;
        break;
      case 13:
        cycles=atoi(optarg);
        break;
      default:
        fprintf(stderr, "in file %s function %s() line %d, unknown option! Exiting program.", 
                __FILE__, __FUNCTION__, __LINE__);
        exit(failure);
        break;    
    }
  }

  if(!strcmp(cmdName, "FoldDesign"))
  {
      bool zipFlag=false;
      func.FoldDesignREMC(libDir,dataDir,inputFile,weightFile,contactConstrFile,
                          distanceConstrFile,fixedFile,cycles,randomNum,fixFlag);
  }
  else if(!strcmp(cmdName,"GenerateFragments"))
  {
    bool zipFlag=false;
    func.generateFragments(libDir,dataDir,inputFile,templateList,zipFlag);
  }
  else if(!strcmp(cmdName, "FoldDesignScore"))
  {
    func.FoldDesignScore(libDir,dataDir,inputFile,weightFile,contactConstrFile,
                         distanceConstrFile,pdbFile);
  }
  else
  {
    FoldDesign_help();
  }

}

bool CheckCommandName(
  char* queryname
)
{
  int MAX_CMD_NUM=100;
  char *supportedCMD[]={
    "FoldDesignScore",
    "GenerateFragments",
    "FoldDesign",
    "Refine",
    NULL
  };

  bool doesExist=false;
  for(int i=0;i<MAX_CMD_NUM;i++)
  {
    if(supportedCMD[i]==NULL) 
    {
      break;
    }
    else
    {
      if(strcmp(queryname,supportedCMD[i])==0)
      {
        doesExist=true;
        break;
      }
    }
  }
  return doesExist;
}

void FoldDesign_interface()
{
  printf(
    "#######################################################################################\n"
    "                                    FoldDesign                                         \n"
    "                   A program for de novo protein fold design                           \n"
    "\n\n"
    "  Written by Robin Pearce (robpearc@umich.edu)\n"
    "  Copyright (C) Regents of the University of Michigan\n"
    "  Prof. Yang Zhang's Research Group\n"
    "  Dept. of Computational Medicine & Bioinformatics\n"
    "  Medical School\n"
    "  University of Michigan\n"
    "  Ann Arbor, MI 48109-2218, USA\n"
    "#######################################################################################\n"
    );
}

void FoldDesign_help()
{
  printf(
    "FoldDesign program options:\n\n"
    "Basic OPTIONS:\n"
    "Usage: FoldDesign [OPTIONS] --command=arg\n"
    "   --help          print help message\n"
    "   --command=arg   choose your program function:\n"
    "                   GenerateFragments\n"
    "                   FoldDesign\n"
    "   --data_dir      Path to the directory containing the input.txt file,\n"
    "                   distance restraints (distance_restr.txt), and contact\n"
    "                   restraints (contact_restr.txt)\n"
    "   --library       Path to the FoldDesign library\n"
    "   --random_num    Random number used during the simulations. Default 102\n"
    "   --remc_cycles   Number of REMC cycles to perform. Default 500.\n"
    );
}


