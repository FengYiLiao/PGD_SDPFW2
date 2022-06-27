%formulate dual problem in sedumi "primal" form

function [A_new,b_new,c_new,K_new]= FormStabelSetProblem_DualSDP(Adj)
%     rng(1);    
    n = width(Adj);
    %C = Adj;
    I = eye(n);
    J = ones(n,n);
    [IndSym,IndDiag,IndOffDiag,ShrinkIndDiag,ShrinkIndOffDiag,IndOffDiagCounter] = SymmetricIndices(n);
    Adj = -(Adj-ones(n,n));
    Adj_tril = tril(Adj,1);
    [row,col] = find(Adj_tril);
    m = length(row);
    Al = zeros(length(row),n^2);
    for i = 1:length(row)
        a1 = sub2ind([n,n],row(i),col(i));
        a2 = sub2ind([n,n],col(i),row(i));
        Al(i,a1) = 1;
        Al(i,a2) = 1;
    end

    
%     [E,Comb,D]= ConstrOPR(P); %Construct Operator and the Combination
%    NumOfComb = length(Comb);
    
    rows = [];
    cols = [];
    vals = [];
    
    NumOfVars = 0;
    K.f = 1+n^2;
    K.s = [n];
    start = 1;
%     for num = 1:NumOfComb
%         fprintf("%d\n",num);
%         i = Comb(num,1); j = Comb(num,2);
%         EE = ([E{i};E{j}]');
%         w = width(EE);
%         
%         
% %         Acone(:,start:start+w^2-1) =  -kron(EE,EE);
%         cons = -kron(EE,EE);
%         [r,l,v]=find(cons);
%         rows = [rows; r];
%         cols = [cols;l+start-1];
%         vals = [vals;v];
%         
%         
%         NumOfVars = NumOfVars +w^2;
%         K.s = [K.s,w];
%         start = start+w^2;
%     end
%     Acone =sparse(rows',cols',vals',n^2,sum(D));
    Acone = -eye(n^2);

    A = [-vec(I),eye(n^2),Acone];
    A = [A;zeros(height(Al),1),Al,zeros(height(Al),width(Acone))];
    %A = [A;sparse(NumOfVars,1+n^2),eye(NumOfVars)];
    b = zeros(n^2+height(Al),1);
    %b = zeros(,)
    b(1:n^2) = vec(J);
    
    c = zeros(1+n^2+n^2,1);
    c(1) = -1;
    %[x,y,info] = sedumi(A,b,c,K);
    
    A_new = A;
    b_new = b;
    c_new = c;
    K_new = K;

end