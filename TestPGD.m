%TestPGD
%We first consider homogeneous partition
clc;clear; 
dx = 2;
load('n10m10.mat');
At_sdp = full(At_sdp); b_sdp = full(b_sdp); c_sdp = full(c_sdp);
opts.Maxiter =1000;
opts.dx = dx;  %size of a block
opts.t = 0.03; %barrier function parameter
opts.a = 0.2;  %linear search parameter
opts.beta = 0.5; %linear search parameter
opt = ColumnGen_Both(At_sdp,b_sdp,c_sdp,K_sdp,dx,1,true); %True solution
[obj_inner,X0] =  InnerApproximation(At_sdp,b_sdp,c_sdp,K_sdp,dx,1); %initial point
[OBJ_Barr,Obj] = PGD(At_sdp,b_sdp,c_sdp,K_sdp,opts,vec(X0),opt); 