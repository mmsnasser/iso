function [cap,et,alpha,z2] = hyptricap(ver,alpha,z2,n)
% Mohamed Nasser, 25-03-2019
% This MATLAB function compute the capacity of the hyp triangle with the
% vertices a,b and c. 
% Let G be the domain inside the unit circle and outside the hyp polygon,
% alpha is a point in G
% z2 is a point inside the polygon.
%
% parametrization of the unit circle
t          =  (0:2*pi/n:2*pi-2*pi/n).';
et(1:n,1)  =   exp(i.*t);
etp(1:n,1) =  i.* exp(i.*t);
% parametrization the polygon
m          = length(ver);ver(m+1)=ver(1); 
[s,sp]     =   deltw(t,m,3);
for k=1:m
    aa        = ver(k);  bb = ver(k+1);
    [cent(k),rd(k)] = my3Pts2(aa,bb,bb/(abs(bb)^2));
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
[~,cap]  =  annq (et,etp,n,alpha,z2,'b');
% plot the domain
figure(1)
hold on; box on
crv=et(1:n); crv(n+1)=crv(1);
plot(real(crv),imag(crv),'k','LineWidth',1.5)
crv=et(n+1:2*n); crv(n+1)=crv(1);
plot(real(crv),imag(crv),'b','LineWidth',1.5)
plot(real(ver),imag(ver),'or')
plot(real(alpha),imag(alpha),'or','markerFaceColor','r')
plot(real(z2),imag(z2),'ob','markerFaceColor','b')
axis equal
axis([-1.05 1.05 -1.05 1.05])
drawnow
end