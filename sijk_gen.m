function[sijk]=sijk_gen(si)
[n,m]=size(si);
sijk=zeros(m,m,m);

for i = 1:m
    for j=1:m
        for k = 1:m
            if(k>j && j>i)
                sijk(i,j,k) = si(i)*si(j)*si(k);
            end
        end
    end
end
    