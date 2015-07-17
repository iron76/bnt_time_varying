iCub

nBody = iCub_dmodel.NB;

aaq   = zeros(nBody,1);

% Do the computations in the CAD configuration (all 0, but second joint of
% the shoulder set to 90)
for i = 1 : iCub_dmodel.NB
   if strcmp(iCub_dmodel.jointname{i}, 'l_shoulder_roll') 
       aaq(i) = pi/2;
   end
   if strcmp(iCub_dmodel.jointname{i}, 'r_shoulder_roll') 
       aaq(i) = pi/2;
   end
end

mapObj = containers.Map(iCub_dmodel.linkname,1:nBody);

aaXup = cell(nBody,1);
aaXa = cell(nBody,1);

%% Compute transforms with respect to the parent for all links
         % obj.Xup{i} contains {}^{i}X_{\lambda(i)}
         for i = 1 : nBody
            [ XJ, ~ ] = jcalc( iCub_dmodel.jtype{i}, aaq(i) );
            aaXup{i} = XJ * iCub_dmodel.Xtree{i};
         end

%% Compute transforms with respect to the base for all links
         % obj.Xa{i} contains {}^{i}X_{0}
         for i = 1:length(iCub_dmodel.parent)
            if iCub_dmodel.parent(i) == 0
               aaXa{i} = aaXup{i};
            else
               aaXa{i} = aaXup{i} * aaXa{iCub_dmodel.parent(i)};
            end
         end
         
%% Hardcoded transforms for sensors boards
% {}^root X_{sensor} => X_root_sensor 

% torso
X_root_0B7 = pluckerFromSE3Transform([0.1161     -0.762460   -0.636534;            ...
                   -0.0841462 -0.646115    0.758588;            ...
                   -0.989667  -0.03451 -0.139172;],             ...
                  [-0.0505674  -0.0594174   0.0793540]);

X_root_0B8 = pluckerFromSE3Transform([0.107     -0.4775    -0.872;              ...
                   -0.092    -0.8781    0.469473;            ...
                  -0.989954 0.0302051 -0.138126;],           ...
                  [-0.0626  0.04433   0.0793589]);

X_root_0B9 = pluckerFromSE3Transform([0.107     0.4775    -0.872;          ...
                   0.092    -0.8781    -0.469473;       ...
                  -0.989954 -0.0302051 -0.138126;],     ...
                  [-0.0626  -0.04433   0.0793589]);

X_root_0B10 = pluckerFromSE3Transform([0.1161    0.762460   -0.636534;         ...
                    0.0841462 -0.646115  -0.758588;         ...
                   -0.989667  0.03451 -0.139172;],          ...
                  [-0.0505674  -0.0594174   0.0793540]);

% r_arm 
X_root_2B13    = pluckerFromSE3Transform([0.760908  -0.258982 -0.594935;      ...
                       0.204107  0.965882   -0.159412;     ...
                       0.615922  -0.000132127 0.787807;],  ...
                       [-0.05328 0.1697     0.20365]); 

X_root_2B11    = pluckerFromSE3Transform([-0.259073 -0.865864 0.427974;      ...
                        0.965858  -0.232311 0.114675;      ...
                        0.000129  0.443071  0.896486;],    ...
                       [-0.0239259 0.267444 0.20398]);  

X_root_2B12    = pluckerFromSE3Transform([-0.259073 -0.865864 0.427974;        ...
                        0.965858  -0.232311 0.114675;        ...
                        0.000129  0.443071  0.896486;],      ...
                       [-0.005333 0.198797 0.20398]);  

X_root_2B10    = pluckerFromSE3Transform([-0.259073 -0.865864 0.427974;      ...
                        0.965858  -0.232311 0.114675;      ...
                        0.000129  0.443071  0.896486;],    ...
                       [-0.0136571 0.222021 0.20488]);  

X_root_2B9     = pluckerFromSE3Transform([0.258122 0.947610   0.188170;      ...
                        -0.9647  0.242282   0.103211;      ...
                        0.0522136 -0.208169 0.976698;],    ...
                       [-0.0464325 0.3172   0.2023]);

X_root_2B8     = pluckerFromSE3Transform([0.147280  -0.607326 -0.78068;      ...
                        -0.984907 -0.162609 -0.0593074;    ...
                        -0.0909268 0.777634 -0.622108;],   ...
                        [-0.101328  0.3461  0.1582]);

% l_arm (to check!) 
X_root_1B13    = pluckerFromSE3Transform([0.760908  0.258982 -0.594935;      ...
                       -0.204107  0.965882   0.159412;     ...
                       0.615922  0.000132127 0.787807;],  ...
                       [-0.05328 -0.1697     0.20365]); 
                   
X_root_1B11    = pluckerFromSE3Transform([-0.259073 0.865864 0.427974;          ...
                        -0.965858  -0.232311 -0.114675;       ...
                        0.000129  -0.443071  0.896486;],      ...
                       [-0.0239259 -0.267444 0.20398]);  

X_root_1B12    = pluckerFromSE3Transform([-0.259073 0.865864 0.427974;         ...
                        -0.965858  -0.232311 -0.114675;        ...
                        0.000129  -0.443071  0.896486;],        ...
                       [-0.005333 -0.198797 0.20398]);  

X_root_1B10    = pluckerFromSE3Transform([-0.259073 0.865864 0.427974;         ...
                        -0.965858  -0.232311 -0.114675;      ...
                        0.000129  -0.443071  0.896486;],     ...
                       [-0.0136571 -0.222021 0.20488]);  

X_root_1B9     = pluckerFromSE3Transform([0.258122 -0.947610   0.188170;      ...
                        0.9647  0.242282   -0.103211;       ...
                        0.0522136 0.208169 0.976698;],      ...
                       [-0.0464325 -0.3172   0.2023]);

X_root_1B8     = pluckerFromSE3Transform([0.147280  0.607326 -0.78068;        ...
                        0.984907 -0.162609 0.0593074;       ...
                        -0.0909268 -0.777634 -0.622108;],   ...
                        [-0.101328  -0.3461  0.1582]);

% r_mmsp this need to be checked !

X_root_2B7      = pluckerFromSE3Transform([0.172319  0.985039 -0.0019222;      ...
                        -0.984998 0.172329 0.00907318;      ...
                        0.00926869 0.0003298 0.999957;],    ...
                        [0.0 0.0 0.0]'); %[-0.104839 0.457   0.177]');

% l_mmsp this need to be checked ! 

X_root_1B7      = pluckerFromSE3Transform([0.172319  -0.985039 -0.0019222;     ...
                        0.984998 0.172329 -0.00907318;      ...
                        0.00926869 -0.0003298 0.999957;],   ...
                        [0.0 0.0 0.0]'); %[-0.104839 -0.457   0.177]');

                    

% Re-extracted from CAD just to be sure
X_root_imu_mtx = pluckerFromSE3Transform([1.0  0.0 0.0;    ...
                        0.0 1.0 0.0;      ...
                        0.0 0.0 1.0;],   ...
                        [-0.02010 -0.0095 0.375897]);


X_chest_imu_mtx = aaXa{mapObj('chest+torso+neck_1+neck_2+head+imu_frame')}*X_root_imu_mtx;

%%% Convert global coordinates to local coordinates
X_chest_9B7   = aaXa{mapObj('chest+torso+neck_1+neck_2+head+imu_frame')}*X_root_0B7;
X_chest_9B8   = aaXa{mapObj('chest+torso+neck_1+neck_2+head+imu_frame')}*X_root_0B8;
X_chest_9B9   = aaXa{mapObj('chest+torso+neck_1+neck_2+head+imu_frame')}*X_root_0B9;
X_chest_9B10  = aaXa{mapObj('chest+torso+neck_1+neck_2+head+imu_frame')}*X_root_0B10;

X_r_arm_2B11  = aaXa{mapObj('r_upper_arm+r_arm')}*X_root_2B11;
X_r_arm_2B10  = aaXa{mapObj('r_upper_arm+r_arm')}*X_root_2B10;
X_r_arm_2B12  = aaXa{mapObj('r_upper_arm+r_arm')}*X_root_2B12;
X_r_arm_2B13  = aaXa{mapObj('r_upper_arm+r_arm')}*X_root_2B13;

X_l_arm_1B11  = aaXa{mapObj('l_upper_arm+l_arm')}*X_root_1B11;
X_l_arm_1B10  = aaXa{mapObj('l_upper_arm+l_arm')}*X_root_1B10;
X_l_arm_1B12  = aaXa{mapObj('l_upper_arm+l_arm')}*X_root_1B12;
X_l_arm_1B13  = aaXa{mapObj('l_upper_arm+l_arm')}*X_root_1B13;

X_r_forearm_2B8 = aaXa{mapObj('r_forearm+r_wrist_1+r_hand+r_gripper')}*X_root_2B8;
X_r_forearm_2B9 = aaXa{mapObj('r_forearm+r_wrist_1+r_hand+r_gripper')}*X_root_2B9;

X_l_forearm_1B8 = aaXa{mapObj('l_forearm+l_wrist_1+l_hand+l_gripper')}*X_root_1B8;
X_l_forearm_1B9 = aaXa{mapObj('l_forearm+l_wrist_1+l_hand+l_gripper')}*X_root_1B9;

correction_2B7 = [ 1.0 0.0 0.0; ...
                   0.0 0.0 -1.0; ...
                   0.0 1.0 0.0];
X_r_forearm_2B7 = pluckerFromSE3Transform(correction_2B7,[0.0 0.0 0.0])*aaXa{mapObj('r_forearm+r_wrist_1+r_hand+r_gripper')}*X_root_2B7;

correction_1B7 = [ 1.0 0.0 0.0; ...
                   0.0 0.0 1.0; ...
                   0.0 -1.0 0.0];
% correction_1B7 = eye(3,3);
X_l_forearm_1B7 = pluckerFromSE3Transform(correction_1B7,[0.0 0.0 0.0])*aaXa{mapObj('l_forearm+l_wrist_1+l_hand+l_gripper')}*X_root_1B7;

% For the foot is easier to express the MTB directly in the local frame
% (check orientation)
X_r_foot_closer_to_base = pluckerFromSE3Transform([0.0     0.0 1.0;             ...
                                1.0     0.0  0.0;             ...
                                0.0     1.0 0.0;],           ...
                               [0.059  0.0118  -0.07497]);

X_r_foot_far_from_base = pluckerFromSE3Transform([0.0     0.0 -1.0;            ...
                                1.0     0.0  0.0;           ...
                                0.0     -1.0 0.0;],         ...
                               [0.059  0.0008  0.102]);

X_l_foot_closer_to_base = pluckerFromSE3Transform([0.0     0.0 -1.0;           ...
                                -1.0     0.0  0.0;           ...
                                0.0     1.0 0.0;],         ...
                               [-0.059  -0.0118  -0.077]);

X_l_foot_far_from_base = pluckerFromSE3Transform([0.0     0.0 -1.0;          ...
                                -1.0     0.0  0.0;         ...
                                0.0     1.0 0.0;],       ...
                               [-0.059  -0.0008  -0.102]);

% This are manually extracted from urdf 
X_la_fts_correction = [ -1.0 0.0 0.0; ...
                        0.0 0.0 -1.0; ...
                        0.0 1.0 0.0];
X_left_arm_fts =  pluckerFromSE3Transform(X_la_fts_correction,[0.0 0.0 0.0])*urdfpose2spatialplucker([-1.5708 0 0],[0.015 9.36797e-18 0.15228]);
X_ra_fts_correction = [ -1.0 0.0 0.0; ...
                        0.0 0.0 -1.0; ...
                        0.0 1.0 0.0];
X_right_arm_fts  = pluckerFromSE3Transform(X_ra_fts_correction,[0.0 0.0 0.0])*urdfpose2spatialplucker([1.5708 -0 3.14159],[0.015 -8.93429e-18 0.15228]);
X_ll_fts_correction = [ 0.0 1.0 0.0; ...
                        -1.0 0.0  0.0; ...
                        0.0 0.0 -1.0];
X_left_leg_fts  =  pluckerFromSE3Transform(X_ll_fts_correction,[0.0 0.0 0.0])*urdfpose2spatialplucker([3.14159 0 -1.5708],[0 0 0]);
X_rl_fts_correction = [ 0.0 1.0 0.0; ...
                        -1.0 0.0  0.0; ...
                        0.0 0.0 -1.0];
X_right_leg_fts =  pluckerFromSE3Transform(X_rl_fts_correction,[0.0 0.0 0.0])*urdfpose2spatialplucker([0 -0 1.5708],[0 0 0]);
X_lf_fts_correction = [ -1.0 0.0 0.0; ...
                        0.0 -1.0  0.0; ...
                        0.0 0.0 -1.0];
X_left_foot_fts  = pluckerFromSE3Transform(X_lf_fts_correction,[0.0 0.0 0.0])*urdfpose2spatialplucker([0 1.5708 0],[-0.0035 0 0.0685]);
X_rf_fts_correction = [ 1.0 0.0 0.0; ...
                        0.0 -1.0  0.0; ...
                        0.0 0.0 1.0];
X_right_foot_fts = pluckerFromSE3Transform(X_rf_fts_correction,[0.0 0.0 0.0])*urdfpose2spatialplucker([0 1.5708 0],[-0.0035 4.28626e-19 0.0685]);

%%% 
% Get the sensor frame actually needed
% by the main.m script

% featherstone matrices are angular-linear
% our data is linear-angular
permutator = [zeros(3,3),eye(3,3); eye(3,3), zeros(3,3)];

X_chest_imu = permutator*X_chest_imu_mtx*permutator;

% 1B13 is 1:3
X_l_upper_arm_la_acc = permutator*X_l_arm_1B13*permutator;

% 9B3 is ??
X_l_upper_foot_lf_acc = permutator*X_l_foot_closer_to_base*permutator;

% 1B7 is hand mmsp
X_l_forearm_lh_imu    = permutator*X_l_forearm_1B7*permutator*diag([1,1,1,-1,-1,1]); 

% 9B7 is 1:3 
X_chest_to_acc        = permutator*X_chest_9B7*permutator;

% 2B13 is 1:3
X_r_upper_arm_ra_acc = permutator*X_r_arm_2B13*permutator;

% 9B12 is ??
X_r_upper_foot_rf_acc = permutator*X_r_foot_closer_to_base*permutator;

% 2B7 is hand mmsp
X_r_forearm_rh_imu    = permutator*X_r_forearm_2B7*permutator*diag([1,1,1,-1,-1,1]);

% for the fts, convert motion spatial transforms to motion force transforms
% 
X_l_upper_arm_la_fts_force  = permutator*inv(X_left_arm_fts)'*permutator;
X_r_upper_arm_ra_fts_force  = permutator*inv(X_right_arm_fts)'*permutator;
X_l_thigh_ll_fts_force      = permutator*inv(X_left_leg_fts)'*permutator;
X_r_thigh_rl_fts_force      = permutator*inv(X_right_leg_fts)'*permutator;
X_l_upper_foot_lf_fts_force = permutator*inv(X_left_foot_fts)'*permutator;
X_r_upper_foot_rf_fts_force = permutator*inv(X_right_foot_fts)'*permutator;