function [A_new,b_new,c_new,K_new]= FormStabelSetProblem_Primal_Beta(Adj)
%Faster version
    n = width(Adj);
    C = Adj;
    I = eye(n);
    J = ones(n,n);
    %[IndSym,IndDiag,IndOffDiag,ShrinkIndDiag,ShrinkIndOffDiag,IndOffDiagCounter] = SymmetricIndices(n);
    Adj_tril = tril(Adj,-1);
    [row,col] = find(Adj_tril);
    b = zeros(1+length(row),1);
   % K.f = 1+ length(row);
    K.s = n;
    b(1) = 1;
    %A = [];
    %A = [A;vec(I)'];
    A = sparse(1+length(row),n^2);
    A(1,:) = vec(I)';
    for i = 1:length(row)
        a1 = sub2ind([n,n],row(i),col(i));
        a2 = sub2ind([n,n],col(i),row(i));
        A(i+1,a1) = 1;
        A(i+1,a2) = 1;
        %Cons()
%         Cons = zeros(n,n);
%         Cons(row(i),col(i)) = -1;
%         Cons(col(i),row(i)) = -1;
%         A = [A;vec(Cons)'];
    end
    %A = [A;-speye(n^2)];
    
    A_new = A;
    b_new = b;
    c_new = -vec(J);
    K_new = K;

end