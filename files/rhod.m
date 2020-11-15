function d = rhod(x,y)
%
d  =  2*asinh(abs(x-y)/(sqrt(1-abs(x)^2)*sqrt(1-abs(y)^2)));
%
end