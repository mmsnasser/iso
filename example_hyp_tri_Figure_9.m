% Mohamed Nasser, 10-07-2020
% 
% 
clear;clc; 
addpath files fmm
% Choose three vertices of the hyp triangle (in clockwise orientation)
n        =   3*2^12;
t          =  (0:2*pi/n:2*pi-2*pi/n).';
%%
% r for the equilateral triangle with angle theta
s        =   0.7;
% Then, we compute the vertices of the triangle
trv1      =   s*exp(0.0i*pi);
trv2      =   s*exp((4/3)*i*pi);
trv3      =   s*exp((2/3)*i*pi);
% Choose alpha inside the unit circle and outside the equilateral
% triangle T
v = [trv1,trv2,trv3]; ver = v; ver(4)=ver(1); 
%%
et(1:n,1)  =   exp(i.*t);
etp(1:n,1) =  i.* exp(i.*t);
% parametrization the triangle
[spnt,spntp]     =   deltw(t,3,3);
for k=1:3
    aa        = ver(k);  bb = ver(k+1);
    if bb==0
        'error: division by zero'
        cap=[];
        return;
    end
    [cent(k),rd(k)] = my3Pts(aa,bb,bb/(abs(bb)^2));
    ang{k}    = carg([aa-cent(k),bb-cent(k)]);
end
for k=1:3
    alp(k)=ang{k}(1);
    bet(k)=ang{k}(2);
end
for k=1:3
    if rd(k)>10^10
        et (n+1+(k-1)*n/3:n+k*n/3,1) = ver(k)+(3/(2*pi))*(ver(k+1)-ver(k)).*(spnt (1+(k-1)*n/3:k*n/3)-spnt(1+(k-1)*n/3));
        etp(n+1+(k-1)*n/3:n+k*n/3,1) =        (3/(2*pi))*(ver(k+1)-ver(k)).*(spntp(1+(k-1)*n/3:k*n/3));
    else
        thet  =  (3*(bet(k)-alp(k))/(2*pi)).*spnt(1+(k-1)*n/3:k*n/3)+k*alp(k)-(k-1)*bet(k);
        thetp =  (3*(bet(k)-alp(k))/(2*pi)).*spntp(1+(k-1)*n/3:k*n/3);
        et (n+1+(k-1)*n/3:n+k*n/3,1) =   cent(k)+rd(k).*exp(i.*thet);
        etp(n+1+(k-1)*n/3:n+k*n/3,1) =           rd(k).*exp(i.*thet).*(i.*thetp);
    end
end
% plot the domain
%%
alpha = 0.5+0.5*s;
z2    = 0;
[q,cap,zet] = annq (et,etp,n,alpha,z2,'b');
%%
figure
sect=[0,exp(i*linspace(0,2*pi/3)),0];

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

pnt = [ver];
plot(real(pnt),imag(pnt),'ok','MarkerSize',4,'MarkerFaceColor','k')
plot(real(alpha),imag(alpha),'pk','MarkerSize',6,'MarkerFaceColor','k')

seg1 = [0;1];
plot(real(seg1),imag(seg1),':k','LineWidth',1)
seg2 = [0;exp(i*2*pi/3)];
plot(real(seg2),imag(seg2),':k','LineWidth',1)
seg3 = [0;exp(i*4*pi/3)];
plot(real(seg3),imag(seg3),':k','LineWidth',1)

ang1  =   v(1)+0.2.*exp(i.*linspace(2.93,3.35,100));
plot(real(ang1),imag(ang1),'k','LineWidth',1)
ang2  =   v(3)+0.2.*exp(i.*linspace(5,5.47,100));
plot(real(ang2),imag(ang2),'k','LineWidth',1)
ang3  =   v(2)+0.2.*exp(i.*linspace(0.8,1.3,100));
plot(real(ang3),imag(ang3),'k','LineWidth',1)

text( 0.68, 0.07,'$s$','fontsize', 22,'Interpreter','latex');
text(-0.37, 0.68,'$se^{2\pi {\rm i}/3}$','fontsize', 22,'Interpreter','latex');
text(-0.39,-0.68,'$se^{4\pi {\rm i}/3}$','fontsize', 22,'Interpreter','latex');

text( 0.39,-0.01,'$\beta$','fontsize', 20,'Interpreter','latex');
text(-0.26, 0.36,'$\beta$','fontsize', 20,'Interpreter','latex');
text(-0.26,-0.37,'$\beta$','fontsize', 20,'Interpreter','latex');

text( 0.10, 0.27,'$b$','fontsize', 22,'Interpreter','latex');
text( 0.10,-0.31,'$b$','fontsize', 22,'Interpreter','latex');
text(-0.33, 0.00,'$b$','fontsize', 22,'Interpreter','latex');

text(+0.48+0.5*s, 0.07,'$\alpha$','fontsize', 20,'Interpreter','latex');

axis([-1.01 1.01 -1.01 1.01])
axis square
axis off

set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc hyp_tri_equ_fig1
% print -dpdf  hyp_tri_equ_fig1
%%
rad1  = linspace(s+0.001,1-0.001,100);
rad1z = fcau(et,etp,zet,rad1);
rad2  = exp(4i*pi/3).*linspace(s+0.001,1-0.001,100);
rad2z = fcau(et,etp,zet,rad2);
rad3  = exp(2i*pi/3).*linspace(s+0.001,1-0.001,100);
rad3z = fcau(et,etp,zet,rad3);
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

cpnt = [zet(n+1);zet(n+n/3+1);zet(n+2*n/3+1)];
plot(real(cpnt),imag(cpnt),'ok','MarkerSize',4,'MarkerFaceColor','k')
plot(real(falpha),imag(falpha),'pk','MarkerSize',6,'MarkerFaceColor','k')

plot([0,real(rad1z)],[0,imag(rad1z)],':k','LineWidth',1)
plot([0,real(rad2z)],[0,imag(rad2z)],':k','LineWidth',1)
plot([0,real(rad3z)],[0,imag(rad3z)],':k','LineWidth',1)

text( 0.38, 0.02,'$q$','fontsize', 21,'Interpreter','latex');
text(-0.32, 0.32,'$qe^{2\pi {\rm i}/3}$','fontsize', 21,'Interpreter','latex');
text(-0.30,-0.30,'$qe^{4\pi {\rm i}/3}$','fontsize', 21,'Interpreter','latex');

text(real(falpha)-0.15, 0.08,'$f(\alpha)$','fontsize', 20,'Interpreter','latex');

axis([-1.01 1.01 -1.01 1.01])
axis square
axis off

set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc hyp_tri_equ_fig1r
% % print -dpdf  hyp_tri_equ_fig1r
%%
% plot the domain
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

sect=[et(5*n/3+1:2*n).',exp(i*linspace(0,2*pi/3,100))];

hold on; box on
fill(real(sect),imag(sect),[0.85 0.85 0.85])
crv=et(1:n); crv(n+1)=crv(1);
plot(real(crv),imag(crv),'k','LineWidth',1.5)
crv=et(n+1:2*n); %crv(n+1)=crv(1);
plot(real(crv),imag(crv),'b','LineWidth',1.5)

pnt = [ver];
plot(real(pnt),imag(pnt),'ok','MarkerSize',4,'MarkerFaceColor','k')
plot(real(alpha),imag(alpha),'pk','MarkerSize',8,'MarkerFaceColor','k')

seg1 = [0;1];
plot(real(seg1),imag(seg1),':k','LineWidth',1)
seg2 = [0;exp(i*2*pi/3)];
plot(real(seg2),imag(seg2),':k','LineWidth',1)
seg3 = [0;exp(i*4*pi/3)];
plot(real(seg3),imag(seg3),':k','LineWidth',1)

seg1 = [0;1];
plot(real(seg1),imag(seg1),':k','LineWidth',1)
seg2 = [0;exp(i*2*pi/3)];
plot(real(seg2),imag(seg2),':k','LineWidth',1)
seg3 = [0;exp(i*4*pi/3)];
plot(real(seg3),imag(seg3),':k','LineWidth',1)

ang1  =   v(1)+0.2.*exp(i.*linspace(2.93,3.35,100));
plot(real(ang1),imag(ang1),'k','LineWidth',1)
ang2  =   v(3)+0.2.*exp(i.*linspace(5,5.47,100));
plot(real(ang2),imag(ang2),'k','LineWidth',1)
ang3  =   v(2)+0.2.*exp(i.*linspace(0.8,1.3,100));
plot(real(ang3),imag(ang3),'k','LineWidth',1)

text( 0.68, 0.07,'$s$','fontsize', 24,'Interpreter','latex');
text(-0.37, 0.68,'$se^{2\pi {\rm i}/3}$','fontsize', 24,'Interpreter','latex');
text(-0.39,-0.68,'$se^{4\pi {\rm i}/3}$','fontsize', 24,'Interpreter','latex');

text( 0.37,-0.01,'$\beta$','fontsize', 22,'Interpreter','latex');
text(-0.26, 0.35,'$\beta$','fontsize', 22,'Interpreter','latex');
text(-0.26,-0.37,'$\beta$','fontsize', 22,'Interpreter','latex');

text( 0.10, 0.28,'$b$','fontsize', 24,'Interpreter','latex');
text( 0.10,-0.31,'$b$','fontsize', 24,'Interpreter','latex');
text(-0.34, 0.00,'$b$','fontsize', 24,'Interpreter','latex');

text(+0.48+0.5*s, 0.08,'$\alpha$','fontsize', 22,'Interpreter','latex');

axis([-1.01 1.01 -1.01 1.01])
axis square
axis off
set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc hyp_tri_equ_fig2
% print -dpdf  hyp_tri_equ_fig2
%%
% plot the domain
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

sect=[zet(5*n/3+1:2*n).',exp(i*linspace(0,2*pi/3,100))];

hold on; box on
fill(real(sect),imag(sect),[0.85 0.85 0.85])
crv=zet(1:n); crv(n+1)=crv(1);
plot(real(crv),imag(crv),'k','LineWidth',1.5)
crv=zet(n+1:2*n); %crv(n+1)=crv(1);
plot(real(crv),imag(crv),'b','LineWidth',1.5)

cpnt = [zet(n+1);zet(n+n/3+1);zet(n+2*n/3+1)];
plot(real(cpnt),imag(cpnt),'ok','MarkerSize',4,'MarkerFaceColor','k')
plot(real(falpha),imag(falpha),'pk','MarkerSize',8,'MarkerFaceColor','k')
 
plot([0,real(rad1z)],[0,imag(rad1z)],':k','LineWidth',1)
plot([0,real(rad2z)],[0,imag(rad2z)],':k','LineWidth',1)
plot([0,real(rad3z)],[0,imag(rad3z)],':k','LineWidth',1)
 
text( 0.38, 0.02,'$q$','fontsize', 23,'Interpreter','latex');
text(-0.32, 0.32,'$qe^{2\pi {\rm i}/3}$','fontsize', 23,'Interpreter','latex');
text(-0.30,-0.30,'$qe^{4\pi {\rm i}/3}$','fontsize', 23,'Interpreter','latex');

text(real(falpha)-0.15, 0.09,'$f(\alpha)$','fontsize', 22,'Interpreter','latex');

axis([-1.01 1.01 -1.01 1.01])
axis square
axis off
set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc hyp_tri_equ_fig2r
% print -dpdf  hyp_tri_equ_fig2r
%%