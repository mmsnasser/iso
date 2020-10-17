% Mohamed Nasser, 25-03-2019
% To compute the capacity of hyp polygon
% We use the MATLAB function hyptricap.m and annq.m
%%
% This code is used to generate the results in the Table in the paper. 
%%
clear;clc; 
addpath files fmm
% Choose the vertices of the hyp polygon (clockwise oriented)
format short g
ver{1} = [0.6       0.1-0.8i -0.5+0.6i ];
ver{2} = [0.601     0.0-0.6i -0.599      0.0+0.6i];
ver{3} = [0.6       0.1-0.8i -0.5-0.5i -0.5+0.6i 0.5+0.5i];
ver{4} = [0.6       0.1-0.8i -0.5-0.5i -0.8     -0.5+0.6i   0.5+0.5i];
ver{5} = [0.6       0.1-0.8i -0.5-0.5i -0.8     -0.5+0.6i   0.9i      0.5+0.5i];
ver{6} = [0.6       0.5-0.5i  0.1-0.8i -0.5-0.5i -0.8      -0.5+0.6i  0.9i      0.5+0.5i];
ver{7} = [0.7+0.2i  0.7-0.2i  0.4-0.5i -0.8i     -0.4-0.7i -0.7-0.4i -0.8      -0.7+0.3i...
    -0.4+0.7i  0.9i      0.3+0.8i  0.5+0.5i];
%
v  =  ver{7}; 
m  =  length(v);   
n  =  m*2^12;
% Choose alpha inside the unit circle and outside the polygon
alpha1    =  ( 1+max(real(v)))/2; alpha2    =  (-1+min(real(v)))/2;
if abs(alpha1)>abs(alpha2)
    alpha=alpha2;
else
    alpha=alpha1;
end
% Choose z2 inside the polygon
z2       =  0;
% compute the capacity of the domain G
cap      = hyppolycap(v,alpha,z2,n);
% The hyp distance function
rho = @(x,y)(2*asinh(abs(x-y)/(sqrt(1-abs(x)^2)*sqrt(1-abs(y)^2))));
% Find the perimeter of the polygon P
v(m+1)=v(1);
L = 0;
for k=1:m
    L=L+rho(v(k),v(k+1));
end
% Find the verices of the polygon P_0 with equal angles
R  = (-sin(pi/m)+sqrt(sin(pi/m)^2+sinh(L/(2*m))^2))/sinh(L/(2*m));
vs=R*exp(-i*2*pi*(0:m-1)/m); % The vertices must be clockwise oriented
% Choose alpha inside the unit circle and outside the symmetric polygon
% alphas = R*exp(i*pi/m);
alpha1    =  ( 1+max(real(vs)))/2; alpha2    =  (-1+min(real(vs)))/2;
if abs(alpha1)>abs(alpha2)
    alphas=alpha2;
else
    alphas=alpha1;
end
% Choose z2 inside the polygon
z2s      =  0;
% compute the capacity of the domain G_0
caps     = hyppolycap(vs,alphas,z2s,n);
%%
format short g
'vertices'
v.';
%%
format long g
'cap_of_polygon_P  cap_of_polygon_P_0'
[cap caps]
%%