function [beta1,beta2,beta3] = hyp_tri_ang(b1,b2,b3)
% This function compute the angles beta=[beta(1),beta(2),beta(3)] of the
% hyperbolic triangle with the verices b=[b(1),b(2),b(3)].
T        =  @(z,a)((z-a)./(1-conj(a).*z));
omg      =   abs(carg(T(b2,b1)./T(b3,b1)));
beta1    =   Arg(T(b2,b1),angle(T(b3,b1)))-angle(T(b3,b1))
omg      =   abs(carg(T(b3,b2)/T(b1,b2)));
beta2    =   Arg(T(b3,b2),angle(T(b1,b2)))-angle(T(b1,b2))
omg      =   abs(carg(T(b1,b3)/T(b2,b3)));
beta3   =   Arg(T(b1,b3),angle(T(b2,b3)))-angle(T(b2,b3))
% Compute the are of the triangle
end