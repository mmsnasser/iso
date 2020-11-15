% Mohamed Nasser, 10-05-2020
% cojucture: equilateral regular hyperbolic polygon
% We use the MATLAB function hyptricap.m and annq.m
%%
clc; clear all
addpath fmm files
format short g
%
c   =  0.5;
M1  =  sqrt(1+4*pi/c);
%
mv = [3:50];
%
for kk=1:length(mv)
m  =  mv(kk); % number of vertices
n  =  m*2^10;
% Find the verices of the polygon with equal angles
omega=((m-2)*pi-c)/m; % the angles of the polygon
if omega<0.01
    error;
end
r  = tanh(0.5*acosh(cot(pi/m)*cot(omega/2)));
% Choose the vertices of the hyp polygon (clockwise oriented)
v  = r*exp(-i*2*pi*(0:m-1)/m); 
% Choose alpha inside the unit circle and outside the symmetric polygon
alpha1 = (1+max(real(v)))/2; alpha2  = (-1+min(real(v)))/2;
if abs(alpha1)>abs(alpha2)
    alphas=alpha2;
else
    alphas=alpha1;
end
% Choose z2 inside the polygon
z2  = 0;
% compute the capacity of the domain G_0
cap(kk) = hyppolycap(v,alphas,z2,n);
LB(kk)  = 2*pi/log(M1);
end
%%
dif = cap.'-LB.'
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
plot(mv,cap,'-r','LineWidth',1.5)
hold on
plot(mv,LB,'--b','LineWidth',1.5)
%
%
title({'$c=0.5$'},'Interpreter','LaTeX','FontSize',20) % c=0.5
xlabel('Number of vertices: $m$','Interpreter','LaTeX')
legend({'${\rm cap}({\bf D},T_m)$','$2\pi/\log(M_1)$'},'Interpreter','LaTeX')
% 
grid(gca,'minor')
grid on
set(gca, 'XMinorTick','on')
set(gca, 'YMinorTick','on')
ax=gca;
ax.GridAlpha=0.5;
ax.MinorGridAlpha=0.5;
%
axis([0 50 3.8 4.3]);yticks([3.8:0.1:4.3]) % c=0.5

set(gca,'FontSize',22)

set(gca,'LooseInset',get(gca,'TightInset'))
%
print -depsc hyp_area_LB_p5
%%