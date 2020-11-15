% Mohamed Nasser, 10-05-2020
% cojucture: equilateral regular hyperbolic polygon
% We use the MATLAB function hyptricap.m and annq.m
%%
clc; clear all
addpath fmm files
format short g
%
c   =  10;
M2  =  sqrt(1+4*pi^2/(c^2))+2*pi/c;
%
mv = [3:50];
%
for kk=1:length(mv)
m  =  mv(kk); % number of vertices
n  =  m*2^10;
% Find the verices of the polygon with equal angles
r  = (-sin(pi/m)+sqrt(sin(pi/m)^2+sinh(c/(2*m))^2))/sinh(c/(2*m));
v  =   r*exp(-i*2*pi*(0:m-1)/m); % The vertices must be clockwise oriented
% Choose alpha inside the unit circle and outside the symmetric polygon
% alphas = r*exp(i*pi/m);
alpha1    =  ( 1+max(real(v)))/2; alpha2    =  (-1+min(real(v)))/2;
if abs(alpha1)>abs(alpha2)
    alpha=alpha2;
else
    alpha=alpha1;
end
% Choose z2 inside the polygon
z2       =  0;
% compute the capacity of the domain G_0
cap(kk) = hyppolycap(v,alpha,z2,n);
UB(kk)  = 2*pi/log(M2);
end
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
plot(mv,cap,'-r','LineWidth',1.5)
hold on
plot(mv,UB,'--b','LineWidth',1.5)
%
%
title({'$c=10$'},'Interpreter','LaTeX','FontSize',20) % c=10
xlabel('Number of vertices: $m$','Interpreter','LaTeX')
legend({'${\rm cap}({\bf D},T_m)$','$2\pi/\log(M_2)$'},'Interpreter','LaTeX',...
        'location','SouthEast')
% 
grid(gca,'minor')
grid on
set(gca, 'XMinorTick','on')
set(gca, 'YMinorTick','on')
ax=gca;
ax.GridAlpha=0.5;
ax.MinorGridAlpha=0.5;
%
axis([0 50 8.5 11])  % c=10

set(gca,'FontSize',22)

set(gca,'LooseInset',get(gca,'TightInset'))
%
print -depsc hyp_per_UB_10 % c10
%%