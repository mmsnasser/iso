% Mohamed Nasser, 25-03-2019
% To compute the capacity of hyp polygon
% We use the MATLAB function hyptricap.m and annq.m
%%
% This code is used to generate the rsults in the Table {tab:hyp_pol_rm}
% in the paper.
%%
clc; clear all
addpath fmm files
% Choose m: the number of the vertices of the hyp polygon
format long
rv=[0.1:0.1:0.9].';
mv=[3:7];
n =  3*5*7*2^9
for kk=1:length(rv)
    for jj=1:length(mv)
r=rv(kk); m=mv(jj);
[r m]
vs=r*exp(-i*2*pi*(0:m-1)/m); % The vertices must be clockwise oriented
% Choose alpha inside the unit circle and outside the symmetric polygon
% alpha1    =  ( 1+max(real(vs)))/2; alpha2    =  (-1+min(real(vs)))/2;
% if abs(alpha1)>abs(alpha2)
%     alphar=alpha2;
% else
%     alphar=alpha1;
% end
alphar = (r+0.25*(0.8-r)).*exp(i*pi/6);
% Choose z2 inside the polygon
z2r      =  0;
% compute the capacity of the domain G_0
capr(kk,jj)     = hyppolycap(vs,alphar,z2r,n);
    end
end
%%
format short g
%%
[rv capr]
%%