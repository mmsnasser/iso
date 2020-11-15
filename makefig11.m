% Mohamed Nasser, 05-05-2020
% To compute the capacity of hyp triangle
% We use the MATLAB function hyptricap.m and annq.m
%
% This code compute cap(D,T) alongside the upper bound and lower bound 
% This is ginve in terms of s, 
% the vertices of the triangle are: s, s*e^{ia}, s*e^{2ia}
% 
% 
clear;clc; 
addpath files fmm
% Choose three vertices of the hyp triangle (in clockwise orientation)
rv=[0.002,0.005,0.01,0.02,0.035,0.05:0.025:0.9,...
    0.91:0.01:0.95,0.96:0.0025:0.98,0.982].';
n        =   3*2^12;
rho = @(x,y)(2*asinh(abs(x-y)./(sqrt(1-abs(x).^2).*sqrt(1-abs(y).^2))));
%%
for kk=1:length(rv)
% r for the equilateral triangle with angle theta
r        =   rv(kk);
% Then, we compute the vertices of the triangle
trv1      =   r*exp(0.0i*pi);
trv2      =   r*exp((4/3)*i*pi);
trv3      =   r*exp((2/3)*i*pi);
% Choose alpha inside the unit circle and outside the equilateral
% triangle T
v = [trv1,trv2,trv3];
alpha1    =  ( 1+max(real(v)))/2; alpha2    =  (-1+min(real(v)))/2;
if abs(alpha1)>abs(alpha2)
    alphar=alpha2;
else
    alphar=alpha1;
end
% Choose z2 inside the equilateral triangle
z2r      =   0;
% compute the capacity cap(D,T) where D is the unit disk and T is the
% hyperbolic triangle with the vertices trv1,trv2,trv3 
capr(kk,1)     = hyptricap(trv1,trv2,trv3,alphar,z2r,n);
%
end
%%
% Compute upper bound and lower bound for the capacity cap(D,T)
for j=1:length(rv)
    Ubd(j,1) = 3*pi/mu(sqrt(3).*rv(j)./sqrt(rv(j).^4+rv(j).^2+1));
    Lbd(j,1) = 6*pi/mu(rv(j).^3);
end
%%
format short g
[rv Lbd capr Ubd]
%
%
figure
hold on; box on
plot(rv,capr,'-k','LineWidth',1.5)
hold on
plot(rv,Ubd,'-.b','LineWidth',1.5)
plot(rv,Lbd,'--r','LineWidth',1.5)
%
set(gca,'FontSize',18)
xlabel('{$s$}','FontSize',22,'Interpreter','latex');
% 
legend({'${\rm cap}({\bf D},T)$',...
        '$3\pi/\mu\left(\sqrt{3}s/\sqrt{s^4+s^2+1}\right)$',...
        '$6\pi/\mu\left(s^3\right)$'},...
        'FontSize',18,'Interpreter','latex','Location','northwest');
%
% set(gca,'XTick',[0:0.5:3.5],'FontSize',15);
% set(gca,'YTick',[0:2:16]);
% axis([0 3.5 0 17])
grid(gca,'minor')
grid on
set(gca, 'XMinorTick','on')
set(gca, 'YMinorTick','on')
ax=gca;
ax.GridAlpha=0.75;
ax.MinorGridAlpha=0.75;
axis([0 1 0 22])
set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc hyp_tri_cap_bd
% print -dpdf  hyp_tri_cap_bd
%%