% Mohamed Nasser, 25-03-2019
% To compute the capacity of hyp polygon
% We use the MATLAB function hyptricap.m and annq.m
%%
% This code is used to plot Figure 5 in the paper.
%%
clc; clear all
addpath fmm files
% Choose m: the number of the vertices of the hyp polygon
format long
rv=[0.05:0.025:0.95].';
mv=[3:7];
n =  3*5*7*2^9
for kk=1:length(rv)
    for jj=1:length(mv)
        r=rv(kk); m=mv(jj);
        [r m]
        vs=r*exp(-i*2*pi*(0:m-1)/m); 
        % The vertices must be clockwise oriented Choose alpha inside the 
        % unit circle and outside the symmetric polygon 
        alphar = (r+0.25*(0.8-r)).*exp(i*pi/6);
        % Choose z2 inside the polygon
        z2r      =  0;
        % compute the capacity of the domain G_0
        capr(kk,jj)     = hyppolycap(vs,alphar,z2r,n);
    end
end
%%
format long g
[rv capr]
format short g
%%
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

hold on; box on
plot(rv,capr(:,1),'r','LineWidth',1.5)
plot(rv,capr(:,2),'-.b','LineWidth',1.5)
plot(rv,capr(:,3),':k','LineWidth',1.5)
plot(rv,capr(:,4),'--r','LineWidth',1.5)
plot(rv,capr(:,5),'b','LineWidth',1.5)

legend({'$m=3$','$m=4$','$m=5$','$m=6$','$m=7$'},...
        'FontSize',22,'Interpreter','latex','Location','northwest');

xlabel('{$r$}','FontSize',22,'Interpreter','latex');
ylabel('${\rm cap}(D,P_0)$','FontSize',22,'Interpreter','latex');

set(gca,'XTick',[0:0.2:1],'FontSize',22);
set(gca,'YTick',[0:5:30]);
axis([0 1 0 31])
grid(gca,'minor')
grid on
set(gca, 'XMinorTick','on')
set(gca, 'YMinorTick','on')
ax=gca;
ax.GridAlpha=0.75;
ax.MinorGridAlpha=0.75;
set(gca,'FontSize',22)

set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc fig_hyp_pg_capr
%%