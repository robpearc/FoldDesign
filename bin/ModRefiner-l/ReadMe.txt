                                                                    ModRefiner Instruction
1. main-chain energy minimization
usage: ./mcrefinement data_dir bin_dir ini_name ref_name ran_num
output: ncaco.pdb (initial main-chain model) mcini_name (refined main-chain model)
2. full-atomic energy minimization
usage: ./emrefinement data_dir bin_dir ini_name ref_name str_val ran_num
output: fulinit.pdb (initial full-atomic model) emini_name (refined full-atomic model)

data_dir: data directory which contains ini_name.
bin_dir : directory which contains all the library files.
ini_name: model to be refined.
ref_name: reference model which is also in data_dir. Only Ca atoms are used as restraints. The model doesn't require full-length and continuous backbone.
str_val : strength value in [0,100]. Larger value makes the final model closer to the reference model.
ran_num : random number in [0,1000000].

Recommandation: 
Even you already have the main-chain or full-atomic model, it is better to run main-chain refinement first if the model has bad main-chain quality.
