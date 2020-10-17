function y = theta23(kind, z , q , tol)
%THETA 	computes the theta functions
nmax = max(ceil(log(tol)/log(q)),20);
s=0;
if ( kind == 2)
    for j=nmax:-1:0
        s =s+q^((j+1/2)^2);
    end
    s = 2*s;
end
if ( kind == 3)
   for j=nmax:-1:1
     s =s+q^(j^2);
   end
   s = 1+2*s;
end
y = s;
end