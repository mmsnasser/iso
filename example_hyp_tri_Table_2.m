% Mohamed Nasser, 25-03-2019
% To compute the capacity of hyp triangle
% We use the MATLAB function hyptricap.m and annq.m
%%
% In this code, we generate the results in the table {tab:hyp_tri} in the
% paper.
%%
clc; clear all
addpath fmm files
% Choose three vertices of the hyp triangle (in clockwise orientation)
format short g
vert=[
 0.6    0.2-0.5i     -0.3-0.5i
%  0.7    0.2-0.5i     -0.3-0.5i
%  0.8    0.2-0.5i     -0.3-0.5i
 0.9    0.2-0.5i     -0.3-0.5i
%  0.3    0.2-0.5i     -0.3-0.5i
 0.3i   0.3-0.5i     -0.3-0.5i
%  0.3i   0.3-0.4i     -0.3-0.4i
%  0.3i   0.4-0.3i     -0.4-0.3i
 0.5i   0.25-0.4i    -0.25-0.4i
%  0.5i   0.43-0.25i   -0.43-0.25i
%  0.5i   0.44-0.25i   -0.44-0.25i
%  0.6i   0.52-0.3i    -0.52-0.3i 
%  0.7i   0.61-0.35i   -0.61-0.35i
%  0.8i   0.69-0.4i    -0.69-0.4i 
 0.9i   0.78-0.45i   -0.78-0.45i
 0.95i  0.7-0.4i     -0.5-0.8i   
%  0.4i   0.35-0.2i    -0.35-0.2i 
%  0.3i   0.26-0.15i   -0.26-0.15i
 0.2i   0.17-0.1i    -0.17-0.1i 
 0.1i   0.087-0.05i  -0.087-0.05i
 -0.1i  0.5-0.5i     -0.5-0.5i   
%  -0.1i  0.6-0.5i     -0.6-0.5i   
 -0.1i  0.7-0.5i     -0.7-0.5i   
 ];
b1v=vert(:,1);b2v=vert(:,2);b3v=vert(:,3);
for kk=1:length(b1v)
    b1=b1v(kk);b2=b2v(kk);b3=b3v(kk);
% Choose alpha inside the unit circle and outside the triangle
alpha    =  -b2;
% Choose z2 inside the triangle
z2       =  0.75*(b1+b2+b3)/3;
% Choose n
n        =   3*2^12;
% compute the capacity of the domain G
cap(kk,1)      = hyppolycap(vert(kk,:),alpha,z2,n);
% We find the three angles of the triangle, alphA, alpB, and alpC
beta1 = hyp_ang(b3,b1,b2);
beta2 = hyp_ang(b1,b2,b3);
beta3 = hyp_ang(b2,b3,b1);
% Compute the are of the triangle
Area(kk,1)     =   pi-(beta1+beta2+beta3);
% Compute theta, the average of the three angles
beta     =  (beta1+beta2+beta3)/3;
% To have an equilateral triangle with angle theta, compute first r
r        =   sqrt((1-sqrt(3)*tan(beta/2))/(1+sqrt(3)*tan(beta/2)));
% Then, we compute the vertices of the triangle
bo1      =   r*exp(0.0i*pi);
bo2      =   r*exp((4/3)*i*pi);
bo3      =   r*exp((2/3)*i*pi);
% Choose alpha inside the unit circle and outside the equilateral triangle
alphao   =  -bo2;
% Choose z2 inside the equilateral triangle
z2o      =   0;
% compute the capacity of the domain G0
verto    =  [bo1,bo2,bo3];
capo(kk,1)     = hyppolycap(verto,alphao,0,n);
%%
end
%%
format short g
[b1v b2v b3v]
%%
format long g
[cap capo]
%%

