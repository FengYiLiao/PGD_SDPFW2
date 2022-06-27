function Homo = ChechHomo(K)
%Check if the positive semidefinite matrices are in the same dimension 
    if ~isfield(K, 's')
        return;
    end
    Homo = true;
    len= length(K.s);
    w = K.s(1);
    for i = 1:len
        if K.s(i) ~= w
            Homo = false;
            break;
        end
    end
end