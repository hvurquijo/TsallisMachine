function[C]=qsum(A,B,q)
[a1,a2] = size(A);
[b1,b2] = size(B);

if (a1==b1 && a2==b2)
    C = zeros(a1,a2);
    for i=1:a1
        for ii= 1:a2
            C(i,ii) = A(i,ii)+B(i,ii)+(1-q)*A(i,ii)*B(i,ii);
        end
    end
end