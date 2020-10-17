% Mohamed Nasser, 25-03-2019
% To compute the capacity of hyp polygon
% We use the MATLAB function hyptricap.m and annq.m
%%
% This code used to plot the two figures "hyp_pol"
% commented in the paper.
%%
clc; clear all
addpath fmm files
% Choose the vertices of the hyp polygon (clockwise oriented)
format short g
v = [0.6 0.1-0.8i  -0.5-0.5i -0.5+0.6i 0.5+0.5i];
m=length(v);   n =  m*2^12;
% Choose alpha inside the unit circle and outside the polygon
alpha1    =  ( 1+max(real(v)))/2; alpha2    =  (-1+min(real(v)))/2;
if abs(alpha1)>abs(alpha2)
    alpha=alpha2;
else
    alpha=alpha1;
end
% Choose z2 inside the polygon
z2       =  0;
% compute the capacity of the domain G
t          =  (0:2*pi/n:2*pi-2*pi/n).';
et(1:n,1)  =   exp(i.*t);
etp(1:n,1) =  i.* exp(i.*t);
% parametrization the polygon
v(m+1)=v(1); 
[s,sp]     =   deltw(t,m,3);
for k=1:m
    aa        = v(k);  bb = v(k+1);
    [cent(k),rd(k)] = my3Pts(aa,bb,bb/(abs(bb)^2));
    ang{k}    = carg([aa-cent(k),bb-cent(k)]);
end
for k=1:m
    alp(k)=ang{k}(1);
    bet(k)=ang{k}(2);
end
for k=1:m
    thet  =  (m*(bet(k)-alp(k))/(2*pi)).*s(1+(k-1)*n/m:k*n/m)+k*alp(k)-(k-1)*bet(k);
    thetp =  (m*(bet(k)-alp(k))/(2*pi)).*sp(1+(k-1)*n/m:k*n/m);
    et (n+1+(k-1)*n/m:n+k*n/m,1) =   cent(k)+rd(k).*exp(i.*thet);
    etp(n+1+(k-1)*n/m:n+k*n/m,1) =           rd(k).*exp(i.*thet).*(i.*thetp);
end
% Computing the capacity using the MATLAB function annq
[~,cap]  =  annq (et,etp,n,alpha,z2,'b')
%%
figure
hold on; box on
crv=et(1:n); crv(n+1)=crv(1);
plot(real(crv),imag(crv),'k','LineWidth',1.5)
crv=et(n+1:2*n); %crv(n+1)=crv(1);
plot(real(crv),imag(crv),'b','LineWidth',1.5)
plot(real(alpha),imag(alpha),'pk','markerFaceColor','k')
plot(real(z2),imag(z2),'ok','markerFaceColor','k')
text(-0.725,0.075,'{$\alpha$}','FontSize',26,'Interpreter','latex')
text(+0.05,0.05,'{$z_2$}','FontSize',26,'Interpreter','latex')
axis equal
axis([-1.05 1.05 -1.05 1.05])
axis off
drawnow
set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc fig_hyp_pol1
%%
% The hyp distance function
rho = @(x,y)(2*asinh(abs(x-y)/(sqrt(1-abs(x)^2)*sqrt(1-abs(y)^2))));
% Find the perimeter of the polygon P
v(m+1)=v(1);
L = 0;
for k=1:m
    L=L+rho(v(k),v(k+1));
end
% Find the verices of the polygon P_0 with equal angles
R  = (-sin(pi/m)+sqrt(sin(pi/m)^2+sinh(L/(2*m))^2))/sinh(L/(2*m));
vs=R*exp(-i*2*pi*(0:m-1)/m); % The vertices must be clockwise oriented
% Choose alpha inside the unit circle and outside the symmetric polygon
% alphas = R*exp(i*pi/m);
alpha1    =  ( 1+max(real(vs)))/2; alpha2    =  (-1+min(real(vs)))/2;
if abs(alpha1)>abs(alpha2)
    alphas=alpha2;
else
    alphas=alpha1;
end
% Choose z2 inside the polygon
z2s      =  0;
% compute the capacity of the domain G_0
zet(1:n,1)  =   exp(i.*t);
zetp(1:n,1) =  i.* exp(i.*t);
% parametrization the polygon
vs(m+1)=vs(1); 
[s,sp]     =   deltw(t,m,3);
for k=1:m
    aa        = vs(k);  bb = vs(k+1);
    [cent(k),rd(k)] = my3Pts(aa,bb,bb/(abs(bb)^2));
    ang{k}    = carg([aa-cent(k),bb-cent(k)]);
end
for k=1:m
    alp(k)=ang{k}(1);
    bet(k)=ang{k}(2);
end
for k=1:m
    thet  =  (m*(bet(k)-alp(k))/(2*pi)).*s(1+(k-1)*n/m:k*n/m)+k*alp(k)-(k-1)*bet(k);
    thetp =  (m*(bet(k)-alp(k))/(2*pi)).*sp(1+(k-1)*n/m:k*n/m);
    zet (n+1+(k-1)*n/m:n+k*n/m,1) =   cent(k)+rd(k).*exp(i.*thet);
    zetp(n+1+(k-1)*n/m:n+k*n/m,1) =           rd(k).*exp(i.*thet).*(i.*thetp);
end
% Computing the capacity using the MATLAB function annq
[~,cap]  =  annq (zet,zetp,n,alpha,z2,'b')
%%
figure
hold on; box on
crv=zet(1:n); crv(n+1)=crv(1);
plot(real(crv),imag(crv),'k','LineWidth',1.5)
crv=zet(n+1:2*n); %crv(n+1)=crv(1);
plot(real(crv),imag(crv),'b','LineWidth',1.5)
plot(real(alphas),imag(alphas),'pk','markerFaceColor','k')
plot(real(z2s),imag(z2s),'ok','markerFaceColor','k')
text(-0.725,0.05,'{$\alpha$}','FontSize',26,'Interpreter','latex')
text(+0.05,0.05,'{$z_2$}','FontSize',26,'Interpreter','latex')
axis equal
axis([-1.05 1.05 -1.05 1.05])
axis off
drawnow
set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc fig_hyp_pol2
%%
