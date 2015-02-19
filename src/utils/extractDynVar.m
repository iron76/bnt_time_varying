function [a,fB,f,tau,fx,d2q] = extractDynVar(NB, d)

%EXTRACTDYNVAR Extract variables from the vector of dynamic variables d.
%    This fuction extract dynamic variables from the vector d
%    containing them. Inputs of the function are: 

%    NB - the number of rigid links in the articulated rigid body
%    d  - the vector of dynamic variables
%
%    The function works only if the dimension n of tau and d2q is equal 
%    to NB!


a = cell(1,NB);
for i = 1 : NB
    a{1,i} =  d(((i-1)*26)+1 : ((i-1)*26)+6 ,1);
end

fB = cell(1,NB);
for i = 1 : NB
    fB{1,i} =  d(((i-1)*26)+7 : ((i-1)*26)+12 ,1);
end

f = cell(1,NB);
for i = 1 : NB
    f{1,i} =  d(((i-1)*26)+13 : ((i-1)*26)+18 ,1);
end

tau_ = zeros(NB,1);
for i = 1 : NB
    tau(i,1) = d(((i-1)*26)+19);
end

fx = cell(1,NB);
for i = 1 : NB
    fx{1,i} = d(((i-1)*26)+20 : ((i-1)*26)+25 ,1);
end

d2q = zeros(NB,1);
for i = 1 : NB
    d2q(i,1) = d(((i-1)*26)+26);
end

end