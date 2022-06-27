%StableSet
function [OBJ_Inner,OBJ_Outer,OBJ_sdp]=StableSet(Adj,dx)
    %filename = 'n30p2'
    %Adj = readmatrix(['SedumiData\StableSet\',filename,'.txt']);%adjacency matrix 

    [At_sdp,b_sdp,c_sdp,K_sdp]=FormStabelSetProblem_Primal_Beta(Adj);
    At_sdp = full(At_sdp); b_sdp = full(b_sdp); c_sdp = full(c_sdp);
    At_sdp_copy = At_sdp;
    c_sdp_copy = c_sdp;
    prob1 = SedumiToMosek(At_sdp,b_sdp,c_sdp,K_sdp);
    [rcode1, res1] = mosekopt('minimize info', prob1);
    OBJ_sdp =  -res1.sol.itr.pobjval;

    %%
    OBJ_Inner  = [];
    n = width(Adj);
    [hei,wid] = size(At_sdp);
    m = hei;
    %dx = 5;
    NumOfComb = (n/dx)*(n/dx-1)/2;
    U = eye(n);
    opts.bfw = 1; opts.nop= n/dx; opts.dual = 0;
    Indices =BIGPSDposition(n,dx);
    IndicesAll =BIGPSDpositionAll(n,dx);
    len = (2*dx)*(2*dx+1)/2;
    for iter = 1:11
        tstart1=tic;
        [Anew, bnew, cnew, Knew, info] = factorwidth(At_sdp,b_sdp,c_sdp,K_sdp,opts,IndicesAll);

        prob =SedumiToMosek(Anew,bnew,cnew,Knew);
        [r, res] = mosekopt('minimize info', prob);
        vecx = res.sol.itr.barx;
        OBJ_Inner = [OBJ_Inner, -res.sol.itr.pobjval];
        X = zeros(n);
        tempx = zeros(n);
        start = 1;
        for num = 1:NumOfComb
            X(Indices(num,:)) = X(Indices(num,:)) + vecx(start:start+len-1)';
            start = start + len;
            tempx(Indices(num,:)) = 0; %reset
        end  
        X = X+tril(X,-1)';
        X = U'*X*U;    
        U = chol(X);
        c_sdp = vec(U*mat(c_sdp_copy)*U');
        for i = 1:m
            At_sdp(i,:) = vec(U*(mat(At_sdp_copy(i,:)))*U')';
        end    
    end

    %%
    clearvars -except OBJ_Inner filename OBJ_sdp dx Adj n
    %%
    %% Outer

    J = ones(n,n);
    P = dx*ones(1,n/dx);
    opts.bfw = 1; opts.nop= n/dx; opts.dual = 0;
    [E,Comb,D]= ConstrOPR(P); %Construct Operator and the Combination

    NumOfComb = height(Comb);

    [At_sdp,b_sdp,c_sdp,K_sdp]= FormStabelSetProblem_DualSDP(Adj);
    [Anew, bnew, cnew, Knew, info] = factorwidth_general(At_sdp,b_sdp,c_sdp,K_sdp,opts);
    Anew_copy = Anew;
    Acon = Anew(:,Knew.f+1:end);
    Acon_copy = Acon;
    U = eye(n);
    OBJ_Outer = [];
    Indices =BIGPSDposition(n,dx);
    IndicesAll =BIGPSDpositionAll(n,dx);
    len = (2*dx)*(2*dx+1)/2;
    for iter = 1:11
        %Mosek
        prob1 = SedumiToMosek(Anew,bnew,cnew,Knew);
        [rcode1, res1] = mosekopt('minimize info', prob1);
        start = 1;
        vecx = res1.sol.itr.barx;
        OBJ_Outer = [OBJ_Outer,res1.sol.itr.pobjval];
        X = zeros(n);
        tempx = zeros(n);
        for num = 1:NumOfComb
            X(Indices(num,:)) = X(Indices(num,:)) + vecx(start:start+len-1)';
            start = start + len;
            tempx(Indices(num,:)) = 0; %reset
        end  
        X = X+tril(X,-1)';
        X = U'*X*U;    
        U = chol(X);
        Acone=UpdateCone(Comb,NumOfComb,D,E,U);
        Anew(1:n^2,1+n^2+1:end) = Acone;
    end
end