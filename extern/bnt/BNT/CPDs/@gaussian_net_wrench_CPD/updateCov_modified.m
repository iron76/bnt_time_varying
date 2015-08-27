function CPD = updateCov_modified(newCPD, oldCPD)
% UPDATE_ENGINE Update the engine to take into account the new parameters (inf_engine).
% engine = update_engine(engine, newCPDs)
%
% This generic method is suitable for engines that do not process the parameters until 'enter_evidence'.

CPD = oldCPD;
CPD.cov = newCPD.cov;
CPD.inertial_params = newCPD.inertial_params;
