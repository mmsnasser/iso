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
s        =   0.5;
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

%%
figure
sect=[0,exp(i*linspace(0,2*pi/3)),0];
hold on; box on
fill(real(sect),imag(sect),[0.85 0.85 0.85])
crv=et(1:n); crv(n+1)=crv(1);
plot(real(crv),imag(crv),'k','LineWidth',1.5)

pnt = [ver];
plot(real(pnt),imag(pnt),'ok','MarkerSize',4,'MarkerFaceColor','k')

seg1 = [0;1];
plot(real(seg1),imag(seg1),':k','LineWidth',1)
seg2 = [0;exp(i*2*pi/3)];
plot(real(seg2),imag(seg2),':k','LineWidth',1)
seg3 = [0;exp(i*4*pi/3)];
plot(real(seg3),imag(seg3),':k','LineWidth',1)

seg1 = [0;s];
plot(real(seg1),imag(seg1),'b','LineWidth',1.5)
seg2 = [0;s.*exp(i*2*pi/3)];
plot(real(seg2),imag(seg2),'b','LineWidth',1.5)
seg3 = [0;s.*exp(i*4*pi/3)];
plot(real(seg3),imag(seg3),'b','LineWidth',1.5)

text( 0.48, 0.07,'$s$','FontSize',24,'Interpreter','latex');
text(-0.21, 0.46,'$se^{2\pi {\rm i}/3}$','FontSize',24,'Interpreter','latex');
text(-0.23,-0.45,'$se^{4\pi {\rm i}/3}$','FontSize',24,'Interpreter','latex');

axis([-1.01 1.01 -1.01 1.01])
axis square
axis off
set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc hyp_tri_equ_fig3
print -dpdf  hyp_tri_equ_fig3
%%

%%
figure
sect=[exp(i*linspace(0,pi))];
tcurv = linspace(0,1,100);
ycurv=0.01.*exp(10.*tcurv.*(1-tcurv)).*sin(8*pi.*tcurv);

hold on; box on
fill(real(sect),imag(sect),[0.85 0.85 0.85])
crv=et(1:n); crv(n+1)=crv(1);
plot(real(crv),imag(crv),'k','LineWidth',1.5)

plot(ycurv,tcurv,'r','LineWidth',1.5)

pnt = [-s^(3/2),s^(3/2)];
plot(real(pnt),imag(pnt),'ok','MarkerSize',4,'MarkerFaceColor','k')

seg1 = [-1;1];
plot(real(seg1),imag(seg1),':k','LineWidth',1)

seg1 = [-s^(3/2),s^(3/2)];
plot(real(seg1),imag(seg1),'b','LineWidth',1.5)

text( 0.1,0.45,'$\tilde\Delta$','FontSize',20,'Interpreter','latex');

text( 0.32,-0.12,'$s^{3/2}$','FontSize',24,'Interpreter','latex');
text(-0.52,-0.12,'$-s^{3/2}$','FontSize',24,'Interpreter','latex');

axis([-1.01 1.01 -1.01 1.01])
axis square
axis off
set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc hyp_tri_equ_fig4
print -dpdf  hyp_tri_equ_fig4
%%