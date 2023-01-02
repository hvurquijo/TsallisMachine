function[sij]=sij_gen(si)
[n,m]=size(si);


if m>1 && n<m && n~=1
    sij = cell(1,m);
    for i=1:m
        sij{i} = triu(ones(n),1)*diag(si(:,i)).*triu((ones(n)*diag(si(:,i)))');
    end
elseif  m==1 && n>m
    sij = triu(ones(n),1)*diag(si).*triu((ones(n)*diag(si))');
    
elseif n>1 && n>m && m~=1
    sij = cell(1,n);
    for i=1:n
        sij{i} = triu(ones(m),1)*diag(si(i,:)).*triu((ones(m)*diag(si(i,:)))');
    end
elseif n==1 && n<m
    sij = triu(ones(m),1)*diag(si).*triu((ones(m)*diag(si))');
end
    

