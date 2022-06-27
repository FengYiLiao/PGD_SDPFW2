function Acone = UpdateCone(Comb,NumOfComb,D,E,U)
    %Acone = [];
    n = width(U);
    rows = [];
    cols = [];
    vals = [];
    start = 1;
    for num = 1:NumOfComb
        fprintf("Update: %d\n",num);
        i = Comb(num,1); j = Comb(num,2);
        EE = U'*([E{i};E{j}]');
        w = width(EE);
        cons = -kron(EE,EE);
        [r,l,v]=find(cons);
        rows = [rows; r];
        cols = [cols;l+start-1];
        vals = [vals;v];
        
%         Acone=[Acone, -kron(EE,EE)];
        
        start = start+w^2;
%         NumOfVars = NumOfVars +w^2;
%         K.s = [K.s,w];
    end
    Acone = sparse(rows',cols',vals',n^2,sum(D));
    Acone = full(Acone);
end