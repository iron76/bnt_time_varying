export BTN_FOLDER=/Users/iron/Desktop/iron/myMatlab/graphicalModels
export BTN_MODIFIED_FOLDER=/Users/iron/Desktop/tmp

mkdir -p $BTN_MODIFIED_FOLDER/bnt/BNT/CPDs/@gaussian_CPD
mkdir -p $BTN_MODIFIED_FOLDER/bnt/BNT/inference/@inf_engine
mkdir -p $BTN_MODIFIED_FOLDER/bnt/KPMstats
mkdir -p $BTN_MODIFIED_FOLDER/bnt/modifiedEM

cp $BTN_FOLDER/bnt/BNT/CPDs/@gaussian_CPD/maximize_params_modified.m $BTN_MODIFIED_FOLDER/bnt/BNT/CPDs/@gaussian_CPD/maximize_params_modified.m

cp $BTN_FOLDER/bnt/BNT/CPDs/@gaussian_CPD/update_ess_modified.m $BTN_MODIFIED_FOLDER/bnt/BNT/CPDs/@gaussian_CPD/update_ess_modified.m

cp $BTN_FOLDER/bnt/BNT/CPDs/@gaussian_CPD/updateCov_modified.m $BTN_MODIFIED_FOLDER/bnt/BNT/CPDs/@gaussian_CPD/updateCov_modified.m

cp $BTN_FOLDER/bnt/BNT/inference/@inf_engine/update_engine_modified.m $BTN_MODIFIED_FOLDER/bnt/BNT/inference/@inf_engine/update_engine_modified.m

cp $BTN_FOLDER/bnt/KPMstats/mixgauss_Mstep_modified.m $BTN_MODIFIED_FOLDER/bnt/KPMstats/mixgauss_Mstep_modified.m

cp $BTN_FOLDER/bnt/modifiedEM/clg_Mstep_modified.m $BTN_MODIFIED_FOLDER/bnt/modifiedEM/clg_Mstep_modified.m

cp $BTN_FOLDER/bnt/modifiedEM/learn_params_em_modified.m $BTN_MODIFIED_FOLDER/bnt/modifiedEM/learn_params_em_modified.m