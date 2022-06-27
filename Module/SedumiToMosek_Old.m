function prob = SedumiToMosek(At,b,c,K)
    %Convert data in Sedumi Standard Primal form to Mosek 
    %Update: 04/16/2022
    %*******Important********
    %We only look consider the cone of free(K.f) and semidefinite(K.s)
    %The reason we seperate homogeneous case and unhomogeneous case is just
    %to precipitate the process
    
    %Consider 2 cases. 
    %1. K.s is homogeneous partition. 
    %2. K.s is not homogeneous.
    

    NumOfPSD = length(K.s);
    NumOfFree = K.f;
    [hei,wei] = size(At);
    m = hei; %number of constraints
    At_PSD = At(:,NumOfFree+1:end);
    LenPSD = 0;%total number of PSD variables (only symmetric part)
    steps = K.s.*(K.s+1)./2;
    LenPSD = sum(steps,'all');
  
    
    IntIdxs = zeros(1,LenPSD);

    PSDIDX = zeros(1,LenPSD);
    symrows = zeros(1,LenPSD);
    symcols = zeros(1,LenPSD);
    start = 1;
    OffSet = 0; %For InterestedIdx
    for i = 1:NumOfPSD
        step = K.s(i)*(K.s(i)+1)/2;
        PSDIDX(start:start+step-1) = i;
        [IndSym,IndDiag] = SymmetricIndices(K.s(i));
        [row,col]=ind2sub([K.s(i),K.s(i)],IndSym);
        symrows(start:start+step-1) = row;
        symcols(start:start+step-1) = col;
        
        IntIdx = IndSym';
        if i>1
            OffSet = OffSet + K.s(i-1)^2;
        end
        IntIdx = IntIdx+OffSet;
        IntIdxs(start:start+step-1) =IntIdx;
        
        start = start+step;
    end
    
    
    
    
    %Objective
    c_f = c(1:NumOfFree); %free part of c
    c_s = c(NumOfFree+1:end); %PSD part of c
    
    %Obj_linear part
    prob.c = c_f;
    
    %upper/lower bounds for constraints and variables
    prob.blc = b;
    prob.buc = b;
    prob.blx = -inf*ones(1,NumOfFree);
    prob.bux = [];
    
    % Dimensions of PSD variables
    prob.bardim = K.s;    
    
    [r,c,v]=find(c_s(IntIdxs));
    prob.barc.subj = PSDIDX(r);
    prob.barc.subk = symrows(r);
    prob.barc.subl = symcols(r);
    prob.barc.val = v'; 


    %Constraint
    %constraint_linear part
    if NumOfFree ~= 0
        At_Free = At(:,1:NumOfFree);
        [r,c,v]=find(At_Free);
        prob.a = sparse(r,c,v,m,NumOfFree);
    else
        prob.a = sparse([], [], [], m, 0); 
    end    
    
    %constraint_psd part
    
    [r,c,v] = find(At_PSD(:,IntIdxs));
    prob.bara.subk = symrows(c);
    prob.bara.subl = symcols(c);
    prob.bara.subi = r';
    prob.bara.subj = PSDIDX(c);
    prob.bara.val = v';    
    
    

   
end