function engine = update_engine_modified(engine, newCPDs, oldCPDs)
% UPDATE_ENGINE Update the engine to take into account the new parameters (inf_engine).
% engine = update_engine(engine, newCPDs)
%
% This generic method is suitable for engines that do not process the parameters until 'enter_evidence'.


[m,n] = size(newCPDs);
for i = 1 : n
    engine.bnet.CPD{i} = updateCov_modified(newCPDs{i}, oldCPDs{i});
end