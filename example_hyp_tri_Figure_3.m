% Mohamed Nasser, 25-03-2019
clc; clear all
addpath fmm files
% Choose three vertices of the hyp triangle (in clockwise orientation)
b1       =   0.95i;
b2       =   0.7-0.4i;
b3       =  -0.5-0.8i;

% b1       =  -0.1i;
% b2       =   0.7-0.5i;
% b3       =  -0.7-0.5i;



alpha    =  -b2;
% Choose z2 inside the triangle
z2       =  0.75*(b1+b2+b3)/3;

% Choose alpha inside the unit circle and outside the triangle
% Choose n
n          =   3*2^13;
t          =  (0:2*pi/n:2*pi-2*pi/n).';
et(1:n,1)  =   exp(i.*t);
etp(1:n,1) =  i.* exp(i.*t);
% parametrization the triangle
[s,sp]     =   deltw(t,3,3);
ver      =  [b1 ;  b2 ; b3]; ver(4)=ver(1); 
for k=1:3
    aa        = ver(k);  bb = ver(k+1);
    cent(k)   = my3Pts(aa,bb,bb/(abs(bb)^2));
    rd(k)     = abs(bb-cent(k));
    ang{k}    = carg([aa-cent(k),bb-cent(k)]);
end
for k=1:3
    alp(k)=ang{k}(1);
    bet(k)=ang{k}(2);
end
for k=1:3
    thet  =  (3*(bet(k)-alp(k))/(2*pi)).*s(1+(k-1)*n/3:k*n/3)+k*alp(k)-(k-1)*bet(k);
    thetp =  (3*(bet(k)-alp(k))/(2*pi)).*sp(1+(k-1)*n/3:k*n/3);
    et (n+1+(k-1)*n/3:n+k*n/3,1) =   cent(k)+rd(k).*exp(i.*thet);
    etp(n+1+(k-1)*n/3:n+k*n/3,1) =           rd(k).*exp(i.*thet).*(i.*thetp);
end
% plot the domain
%%
figure
hold on; box on
crv=et(1:n); crv(n+1)=crv(1);
plot(real(crv),imag(crv),'k','LineWidth',1.5)
crv=et(n+1:2*n); %crv(n+1)=crv(1);
plot(real(crv),imag(crv),'b','LineWidth',1.5)
plot(real(alpha),imag(alpha),'pk','markerFaceColor','k')
plot(real(z2),imag(z2),'ok','markerFaceColor','k')
text(-0.65,0.45,'{$\alpha$}','FontSize',26,'Interpreter','latex')
text(+0.10,0.0,'{$z_2$}','FontSize',26,'Interpreter','latex')
axis equal
axis([-1.05 1.05 -1.05 1.05])
axis off
drawnow
set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc fig_hyp_tri1
%%
beta1 = hyp_ang(b3,b1,b2);
beta2 = hyp_ang(b1,b2,b3);
beta3 = hyp_ang(b2,b3,b1);
Area = pi-sum(beta1+beta2+beta3);
%%
beta     =  (beta1+beta2+beta3)/3;
r        =   sqrt((1-sqrt(3)*tan(beta/2))/(1+sqrt(3)*tan(beta/2)));
bo1      =   r*exp(0.0i*pi);
bo2      =   r*exp((4/3)*i*pi);
bo3      =   r*exp((2/3)*i*pi);
alphao   =  -bo2;
z2o      =   0;
zet(1:n,1)  =   exp(i.*t);
zetp(1:n,1) =  i.* exp(i.*t);
% parametrization the triangle
[s,sp]     =   deltw(t,3,3);
ver      =  [bo1 ;  bo2 ; bo3]; ver(4)=ver(1); 
for k=1:3
    aa        = ver(k);  bb = ver(k+1);
    cent(k)   = my3Pts(aa,bb,bb/(abs(bb)^2));
    rd(k)     = abs(bb-cent(k));
    ang{k}    = carg([aa-cent(k),bb-cent(k)]);
end
for k=1:3
    alp(k)=ang{k}(1);
    bet(k)=ang{k}(2);
end
for k=1:3
    thet  =  (3*(bet(k)-alp(k))/(2*pi)).*s(1+(k-1)*n/3:k*n/3)+k*alp(k)-(k-1)*bet(k);
    thetp =  (3*(bet(k)-alp(k))/(2*pi)).*sp(1+(k-1)*n/3:k*n/3);
    zet (n+1+(k-1)*n/3:n+k*n/3,1) =   cent(k)+rd(k).*exp(i.*thet);
    zetp(n+1+(k-1)*n/3:n+k*n/3,1) =           rd(k).*exp(i.*thet).*(i.*thetp);
end
figure
hold on; box on
crv=zet(1:n); crv(n+1)=crv(1);
plot(real(crv),imag(crv),'k','LineWidth',1.5)
crv=zet(n+1:2*n); %crv(n+1)=crv(1);
plot(real(crv),imag(crv),'b','LineWidth',1.5)
plot(real(alphao),imag(alphao),'pk','markerFaceColor','k')
plot(real(z2o),imag(z2o),'ok','markerFaceColor','k')
text(+0.45,0.70,'{$\alpha$}','FontSize',26,'Interpreter','latex')
text(+0.05,0.05,'{$z_2$}','FontSize',26,'Interpreter','latex')
axis equal
axis([-1.05 1.05 -1.05 1.05])
axis off
drawnow
set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc fig_hyp_tri2
%%