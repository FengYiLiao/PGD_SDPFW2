%Change of basis Inner approximate
%Auther: Feng-Yi, Liao

function [OBJ_Inner]=InnerApproximation(At_sdp,b_sdp,c_sdp,K_sdp,dx,MaxIter)
    addpath('Module\');
    n = sqrt(width(At_sdp)); U = eye(n);
    opts.bfw = 1; opts.nop= n/dx; opts.dual = 0;

    NumOfComb = (n/dx)*(n/dx-1)/2;
    m = height(At_sdp);
    OBJ_Inner = [];
    At_sdp_copy = At_sdp; b_sdp_copy = b_sdp; c_sdp_copy = c_sdp;
    ConstrTime = []; SolveTime = []; UpdateTime = [];
    OverallTime = tic ;
    [xIndSym] = SymmetricIndices(2*dx); %Get the indices of symmetric parts
    Indices =BIGPSDposition(n,dx); %The nonzero indices in a n x n matrix (only symmetric part)
    IndicesAll =BIGPSDpositionAll(n,dx);%The nonzero indices in a n x n matrix (including lower and upper)
    len = (2*dx)*(2*dx+1)/2; %len of symmetric part
    fprintf('----start main algorithm-------\n');
    for iter = 1:MaxIter


        tstart1=tic;
        [Anew, bnew, cnew, Knew, info] = factorwidth(At_sdp,b_sdp,c_sdp,K_sdp,opts,IndicesAll);

        prob1 = SedumiToMosek(Anew,bnew,cnew,Knew);
        
        ConstrTime = [ConstrTime, toc(tstart1)];

        tstartSolve=tic;
        [rcode1, res1] = mosekopt('minimize info', prob1);
        status = res1.sol.itr.prosta;
        if ~strcmp(status,'PRIMAL_AND_DUAL_FEASIBLE')
           break; 
        end
        
        OBJ_Inner = [OBJ_Inner, res1.sol.itr.pobjval];
        SolveTime = [SolveTime,toc(tstartSolve)];


        tstartUpdate = tic;
        start = 1;
        vecx = res1.sol.itr.barx;
        X = zeros(n);
        tempx = zeros(n);
        for num = 1:NumOfComb
            X(Indices(num,:)) = X(Indices(num,:)) + vecx(start:start+len-1)';
            start = start + len;
            tempx(Indices(num,:)) = 0; %reset
        end  
        X = X+tril(X,-1)';
        X = U'*X*U;

        %Update U
        [eigvec,eigval] = eig(X); 
        %U = chol(X);
        eigval(eigval<0) = 0;
        U = sqrt(eigval)*eigvec';

        %Update problem
        c_sdp = vec(U*mat(c_sdp_copy)*U');
        for i = 1:m
            At_sdp(i,:) = vec(U*(mat(At_sdp_copy(i,:)))*U')';
        end
        UpdateTime = [UpdateTime,toc(tstartUpdate)];
        time = toc(OverallTime);
        if time >1800 %time limit
           break; 
        end
    end
end