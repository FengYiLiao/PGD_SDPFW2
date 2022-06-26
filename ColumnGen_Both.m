function [OBJ,time] = ColumnGen_Both(At_sdp,b_sdp,c_sdp,K_sdp,dx,MaxIter,Cone)
    addpath('Module\');
    %record original number of free variable
    %it will be used during the construction of X
    if ~isfield(K_sdp,'f')
        K_sdp.f = 0;
    end
    NumOfFreeBefore = K_sdp.f; 
    n = K_sdp.s;
    %n = sqrt(width(At_sdp));
    Cut = 2*dx;
    
    if Cone == true
        opts.bfw = 1; opts.nop= n/dx; opts.dual = 1; opts.keep_split = 0;
    elseif Cone == false
        opts.bfw = 1; opts.nop= n/dx; opts.dual = 0; opts.keep_split = 0;
    end
    %NumOfComb = (n/dx)*(n/dx-1)/2;
    [xIndSym,~,xIndOffDiag,~,~,xIndOffDiagCounter] = SymmetricIndices(n);
    [xIndSym_cut,xIndDiag_cut,xIndOff_cut,xShrinkIndDiag_cut,xShrinkIndOffDiag_cut,xIndOffDiagCounter_cut] = SymmetricIndices(2*dx);

    len_cut = length(xIndSym_cut);

    ZERO = zeros(n,n);
    [Anew, bnew, cnew, Knew, info] = Copy_of_factorwidth(At_sdp,b_sdp,c_sdp,K_sdp,opts);
    [hei,wei] = size(Anew);
    NumOfFreeAfter = Knew.f;
    mnew = hei;
    PSDVarsCount = sum(Knew.s.*Knew.s,'all') ;
    TotalVarCount = wei; 
    OBJ = zeros(MaxIter,1);
    
    %generate data first
    if Cone == true
        At_br = zeros(len_cut,Cut^2);
        ind = sub2ind([len_cut,Cut^2],xShrinkIndDiag_cut',xIndDiag_cut);
        At_br(ind) = -1;
        ind = sub2ind([len_cut,Cut^2],repmat(xShrinkIndOffDiag_cut,1,2),[xIndOff_cut',xIndOffDiagCounter_cut']);
        At_br(ind) = -0.5; %it is important to make the constraint be symmetric       
    end
    
    prob1 = SedumiToMosek(Anew, bnew, cnew, Knew);

    time = [];
    for iter = 1:MaxIter
       
        
        tic;
        [~, res1] = mosekopt('minimize info', prob1);
        T = toc;
        time = [time,T];
        OBJ(iter) = res1.sol.itr.pobjval;
        %vecx = res1.sol.itr.xx(1:Knew.f);
        if Cone == true
            vecx = res1.sol.itr.xx(NumOfFreeBefore+1:NumOfFreeAfter);
                    X = ZERO;  %reset
            X(xIndSym) = vecx;
            X = X + tril(X,-1)';
            [eigvec,eigval] = eig(X);
            if (min(diag(eigval))>=0)
                return;
            end

            mnew = mnew + len_cut;
            PSDVarsCount = PSDVarsCount + Cut^2;
            TotalVarCount = TotalVarCount + Cut^2;
            Knew.s = [Knew.s;Cut];

            [~,I] = sort(diag(eigval));
            eigvec = eigvec(:,I);
            V = eigvec(:,1:Cut)';
            temp = kron(V,V);
            temp(:,xIndOffDiag) = temp(:,xIndOffDiag) + temp(:,xIndOffDiagCounter);
            tempcons = temp(:,xIndSym);
            NewCons = tempcons(xIndSym_cut,:);

            At_button = zeros(len_cut,TotalVarCount);
            %At_button(:,1:Knew.f) = NewCons;
            At_button(:,NumOfFreeBefore+1:Knew.f) = NewCons;
            At_button(:,end-Cut^2+1:end) = At_br;

            badd = zeros(len_cut,1);
            cadd = zeros(TotalVarCount,1); 
            prob2 = SedumiToMosek(At_button, badd, cadd, Knew);

            prob1.bardim = [prob1.bardim;Cut];
            prob1.blc = [prob1.blc;zeros(len_cut,1)];
            prob1.buc = [prob1.buc;zeros(len_cut,1)];
            prob1.a = [prob1.a;prob2.a];
            prob1.bara.subk = [prob1.bara.subk,prob2.bara.subk];
            prob1.bara.subl = [prob1.bara.subl,prob2.bara.subl];
            prob1.bara.subi = [prob1.bara.subi,prob2.bara.subi+mnew-len_cut];%%important
            prob1.bara.subj = [prob1.bara.subj,prob2.bara.subj];
            prob1.bara.val  = [prob1.bara.val,prob2.bara.val];

        elseif Cone == false
            
            y = res1.sol.itr.y;%dual variable
            S = c_sdp - At_sdp'*y;
            %S = c_sdp(Knew.f+1:end) - At_sdp(Knew.f+1:end,Knew.f+1:end)'*y(Knew.f+1:end);
            [eigvec,eigval] = eig(mat(S));
            if (min(diag(eigval))>=0)
                return;
            end
            [~,I] = sort(diag(eigval));
            eigvec = eigvec(:,I);
            V = eigvec(:,1:Cut)';      
            c_add = kron(V,V)*c_sdp;
            %c_add = kron(V,V)*c_sdp(Knew.f+1:end);
            cnew = [cnew;c_add];%should be replaced by mosek after verifiying the idea works
            a_add = At_sdp*kron(V,V)';
            %a_add = At_sdp(Knew.f+1:end,Knew.f+1:end)*kron(V,V)';
            Anew = [Anew, a_add];
            %Anew = [Anew,[zeros(1,width(a_add));a_add]];
            Knew.s = [Knew.s;Cut];
            
            %if iter >1 
            prob1 = SedumiToMosek(Anew, bnew, cnew, Knew);
            %end
        end

                
%         %Anew = [Anew,zeros(mnew-Cut^2,Cut^2)];
%         Anew = [Anew,zeros(mnew-len_cut,Cut^2)];
%         Anew = [Anew;At_button];
%         %bnew = [bnew;zeros(Cut^2,1)];
%         bnew = [bnew;zeros(len_cut,1)];
%         cnew = [cnew;zeros(Cut^2,1)];

    end
end