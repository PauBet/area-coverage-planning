function nullvar = iszero(M)

nullvar = true;
for i=1:size(M,1)
    for j=1:size(M,2)
        if (M(i,j) ~= 0)
            nullvar = false;
            break;
        end
    end
end

end