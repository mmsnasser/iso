function cap = hyptricap(a,b,c,alpha,z2,n)
% Mohamed Nasser, 25-03-2019
% This MATLAB function compute the capacity of the hyp triangle with the
% vertices a,b and c. 
% Let G be the domain inside the unit circle and outside the triangle and
% alpha is a point in G
% z2 is a point inside the triangle.
%
% parametrization of the unit circle
t          =  (0:2*pi/n:2*pi-2*pi/n).';
et(1:n,1)  =   exp(i.*t);
etp(1:n,1) =  i.* exp(i.*t);
% parametrization the triangle
[s,sp]     =   deltw(t,3,3);
ver      =  [a ;  b ; c]; ver(4)=ver(1); 
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
        et (n+1+(k-1)*n/3:n+k*n/3,1) = ver(k)+(3/(2*pi))*(ver(k+1)-ver(k)).*(s (1+(k-1)*n/3:k*n/3)-s(1+(k-1)*n/3));
        etp(n+1+(k-1)*n/3:n+k*n/3,1) =        (3/(2*pi))*(ver(k+1)-ver(k)).*(sp(1+(k-1)*n/3:k*n/3));
    else
        thet  =  (3*(bet(k)-alp(k))/(2*pi)).*s(1+(k-1)*n/3:k*n/3)+k*alp(k)-(k-1)*bet(k);
        thetp =  (3*(bet(k)-alp(k))/(2*pi)).*sp(1+(k-1)*n/3:k*n/3);
        et (n+1+(k-1)*n/3:n+k*n/3,1) =   cent(k)+rd(k).*exp(i.*thet);
        etp(n+1+(k-1)*n/3:n+k*n/3,1) =           rd(k).*exp(i.*thet).*(i.*thetp);
    end
end
% Computing the capacity using the MATLAB function annq
[~,cap]  =  annq (et,etp,n,alpha,z2,'b');
% plot the domain
figure(1)
hold on; box on
crv=et(1:n); crv(n+1)=crv(1);
plot(real(crv),imag(crv),'k','LineWidth',1.5)
crv=et(n+1:2*n); %crv(n+1)=crv(1);
plot(real(crv),imag(crv),'b','LineWidth',1.5)
plot(real(alpha),imag(alpha),'or','markerFaceColor','r')
plot(real(z2),imag(z2),'ob','markerFaceColor','b')
axis equal
axis([-1.05 1.05 -1.05 1.05])
drawnow
end