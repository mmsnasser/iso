% Mohamed Nasser, 10-07-2020
% 
% 
clear;clc; 
addpath files fmm
% Choose three vertices of the hyp triangle (in clockwise orientation)
n        =   5*2^10;
t        =  (0:2*pi/n:2*pi-2*pi/n).';
%%
ver_1    =  [3   ;  1+1.5i ; -1+1.5i ; -2 ;   -1.5i  ];
ver_2    =  [1   ; -i      ; -1      ;  i ];
%%
[et1,etp1]=polygonp(ver_1,n/length(ver_1));
[et2,etp2]=polygonp(ver_2,n/length(ver_2));
%
et   = [et1  ;  et2];
etp  = [etp1 ;  etp2];
%
alpha = 1.25+0.5i;
z2    = 0.25i;
%%
figure

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

hold on; box on
crv=et(1:n); crv(n+1)=crv(1);
fill(real(crv),imag(crv),[0.85 0.85 0.85])
plot(real(crv),imag(crv),'k','LineWidth',1.5)
crv=et(n+1:2*n); %crv(n+1)=crv(1);
fill(real(crv),imag(crv),'w')
plot(real(crv),imag(crv),'b','LineWidth',1.5)

plot(real(alpha),imag(alpha),'pk','MarkerSize',7,'MarkerFaceColor','k')
plot(real(z2),imag(z2),'ok','MarkerSize',6,'MarkerFaceColor','k')



text(real(alpha)+0.08,imag(alpha)+0.08,'$\alpha$','fontsize', 24,'Interpreter','latex');
text(real(z2)+0.08,imag(z2)+0.08,'$z_2$','fontsize', 24,'Interpreter','latex');

axis([-2.01 3.01 -1.51 1.51])
axis square
axis off

set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc domain_G
% print -dpdf  hyp_tri_equ_fig1
%%
[q,cap,zet] = annq (et,etp,n,alpha,z2,'b');
%%
falpha= fcau(et,etp,zet,alpha);
%%
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

hold on; box on
crv=zet(1:n); crv(n+1)=crv(1);
fill(real(crv),imag(crv),[0.85 0.85 0.85])
plot(real(crv),imag(crv),'k','LineWidth',1.5)
crv=zet(n+1:2*n); %crv(n+1)=crv(1);
fill(real(crv),imag(crv),'w')
plot(real(crv),imag(crv),'b','LineWidth',1.5)

plot(real(falpha),imag(falpha),'pk','MarkerSize',7,'MarkerFaceColor','k')
plot(0,0,'ok','MarkerSize',6,'MarkerFaceColor','k')

plot([0,q],[0,0],':k','LineWidth',1)
plot([0,0],[0,1],':k','LineWidth',1)

text(0.3,0.08,'$q$','fontsize', 24,'Interpreter','latex');
text(0  ,0.35,'$1$','fontsize', 24,'Interpreter','latex');

text(real(falpha)-0.18, 0.13,'$f(\alpha)$','fontsize', 24,'Interpreter','latex');

axis([-1.01 1.01 -1.01 1.01])
axis square
axis off

set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc domain_R
%%
