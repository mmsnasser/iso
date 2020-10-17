function y = mu (xx)
%%
x          =  sym(xx);
[K,Kprime] =  ellipk(x);
yy         = (pi/2)*Kprime/K;
y          =  double(yy);
end