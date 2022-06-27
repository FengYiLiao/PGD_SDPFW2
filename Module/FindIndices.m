function Indices = FindIndices(idx,n)
    %find the indices of the correspoing indices in big matrix
    %len means length
    %n is the size of big matrix
    len = length(idx);
    Cols = zeros(1,len*(len+1)/2);
    Rows = zeros(1,len*(len+1)/2);
    %idx = 1:len;
    
    start = 1;
    for num = 1:len
        Rows(start:start+len-num) = idx(num:end);
        Cols(start:start+len-num) = idx(num)*ones(1,len-num+1);
        start = start + (len-(num-1));
    end
    Indices = sub2ind([n,n],Rows,Cols);
end