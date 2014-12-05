function engine = bnet_to_engine(bnet, engine)
% BNET_TO_ENGINE Insert the bnet structure inside the engine (inf_engine)
% engine = bnet_to_engine(bnet, engine)

engine.bnet = bnet;

% We cannot write 'engine.bnet' without writing a 'subsref' function,
% since engine is an object with private parts.
% The bnet field should be the only thing external users of the engine should need access to.
% We do not pass bnet as a separate argument, since it could get out of synch with the one
% encoded inside the engine.
       
