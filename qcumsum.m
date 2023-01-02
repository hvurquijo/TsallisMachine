function[C]=qcumsum(A,q)
[a1,a2]=size(A);
C = zeros(a1,a2);
if (a1==1 && a2>=1)
    C(1,1) = qsum(0,A(1,1),q);
    if(a2>1)
        for i=2:a2
            C(1,i)= qsum(C(1,i-1),A(1,i),q);
        end
    end
else
end