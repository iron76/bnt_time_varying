function res = testRegressors

res = 0;

% create random inertial parameters 
m = rand(1,1);
com = rand(3,1);
ll = rand(3,3);
Icom = ll'*ll;

a = rand(6,1);
v = rand(6,1);

tol = 1e-14;

I = mcI(m,com,Icom);

% test that conversion from and to inertial parameters are consistent
Iconv = inertiaMatrixFromInertialParams(inertialParamsFromInertiaMatrix(I));

if( not(all(isalmost(I,Iconv,tol))) )
    disp('Something wrong in inertial parameters/matrix conversions')
    res = 1;
end

% test the inertia regressor 
momentumDir = I*v;
momentumRegr = inertiaRegressor(v)*inertialParamsFromInertiaMatrix(I);

if( not(all(isalmost(momentumDir,momentumRegr,tol))) ) 
    disp('Something wrong in inertiaRegressor')
    res = 1;
end

% test the net wrench regressor 
momentumDotDir =  I*a + crf(v)*I*v;
momentumDotRegr = netWrenchRegressor(a,v)*inertialParamsFromInertiaMatrix(I);

if( not(all(isalmost(momentumDotDir,momentumDotRegr,tol))) )
    disp('Something wrong in inertiaRegressor')
    res = 1;
end

end