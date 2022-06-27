function [OBJ_Barr,OBJ] = PGD(At_sdp,b_sdp,c_sdp,K_sdp,opts,X0,opt)
    
    dx = opts.dx;
    if ~isfield(K_sdp,'f')
        K_sdp.f = 0;
    end
    n = K_sdp.s;
    OBJ_Barr = zeros(opts.Maxiter,1);
    OBJ =  zeros(opts.Maxiter,1);
    
    %NumOfComb = (n/dx)*(n/dx-1)/2;
    [xIndSym,~,xIndOffDiag,~,~,xIndOffDiagCounter] = SymmetricIndices(n);
    [xIndSym_cut,xIndDiag_cut,xIndOff_cut,xShrinkIndDiag_cut,xShrinkIndOffDiag_cut,xIndOffDiagCounter_cut] = SymmetricIndices(2*dx);
    Paras = PreInatilze(At_sdp,b_sdp,c_sdp,n,dx,opts);
    X = X0;
    for iter = 1:opts.Maxiter
        %disp(iter);
        dX = gradient(X,Paras);
        k = LineSearch(X,Paras);
        X = X - k*projection(dX,Paras);
        OBJ_Barr(iter) = c_sdp'*X - evaluation(X,Paras);
        OBJ(iter) = c_sdp'*X ;
        if mod(iter,50) == 0 || i == 1
            fprintf('%d  %5.2f  %5.2f  %5.2f\n', iter, OBJ(iter),norm(projection(dX,Paras)),(OBJ(iter)-opt)/opt);
        end
        if (norm(projection(dX,Paras),2)<=10^(-3))
            break
        end
    end
end

function dX = gradient(X,Paras)
    IndicesAll = Paras.IndicesAll;
    dX = zeros(Paras.n^2,1);
    Xt = zeros(Paras.n_new);
    for j = 1:Paras.NumOfP
        idx = IndicesAll(j,:);
        Xt(1:Paras.n_new^2) =  X(idx);
        if (min(eig(det(Xt)))<=0)
            dX = zeros(Paras.n^2,1);
            return
        end
        InvXt = inv(Xt);
        dX(idx) = dX(idx)+InvXt(1:Paras.n_new^2)';
    end
    dX = Paras.c_sdp-Paras.t*dX;
end

function Paras = PreInatilze(At_sdp,b_sdp,c_sdp,n,dx,opts)
    %Indices =BIGPSDposition(n,dx); %The nonzero indices in a n x n matrix (only symmetric part)
    IndicesAll =BIGPSDpositionAll(n,dx);%The nonzero indices in a n x n matrix (including lower and upper)
    n_new = dx*2; %dimension for small block
    Paras.IndicesAll = IndicesAll;
    Paras.dx = dx;
    Paras.n = n;
    Paras.n_new = n_new;
    Paras.NumOfP = nchoosek(n/dx,2);%Number of Blocks
    Paras.t = opts.t;
    invAAt = inv(At_sdp*At_sdp');
    Paras.ProjMatrix = eye(n^2) - At_sdp'*invAAt*At_sdp;
   % Paras.v = At_sdp'*invAAt*b_sdp;
  %  Paras.alpha = 0.5;
    Paras.b_sdp = b_sdp;
    Paras.At_sdp = At_sdp;
    Paras.a = opts.a;%0.2;
    Paras.beta = opts.beta;%0.75;
    Paras.c_sdp = c_sdp;
end

function z = projection(X0,Paras)
    z = Paras.ProjMatrix*X0;
end
function val = evaluation(X,Paras)
    val = 0;
    IndicesAll = Paras.IndicesAll;
    sX = zeros(Paras.n,1);
    for j = 1:Paras.NumOfP
        idx = IndicesAll(j,:);
        sX(1:Paras.n_new^2) =  X(idx);
        test = min(eig(mat(sX)));
        if test<=0
            val = val - inf;
        else
            val = val +log(det(mat(sX)));
        end
    end
    val = val*Paras.t;
   % val = +log(det(mat(X)))*Paras.t;
end

function K = LineSearch(X,Paras)
    K = 1;
    while true
       f1 = L(X,Paras,K);
       f2 = G(X,Paras,K);
       if f2<=f1
           break
       else
           K = K*Paras.beta;
       end
       if K < 1.e-8
            warning('Gradient method gets stuck with very small step size!');
            break;
        end
    end
end

function C = L(X,Paras,k)
    A = Cost(X,Paras);
    B = Paras.a*k*norm(projection(gradient(X,Paras),Paras),"fro")^2;
    C = A-B;
end

function C = Cost(X,Paras)
    A = Paras.c_sdp'*X;
    B = evaluation(X,Paras);
    C = A-B;
end

function C = G(X,Paras,k)
    C = Cost(X-k*projection(gradient(X,Paras),Paras),Paras);
end