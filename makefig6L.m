% Mohamed Nasser, 10-05-2020
% cojucture: equilateral regular hyperbolic polygon
%%
clc; clear all
addpath fmm files
format short g
%
c  =  3; % a given hyp area
M1  =  sqrt(1+4*pi/c);
mv = [3,5,7];
for kk=1:length(mv)
m  =  mv(kk); % number of vertices
n{kk}  =  m*2^10;
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
    alpha=alpha2;
else
    alpha=alpha1;
end
% Choose z2 inside the polygon
z2  = 0;
% compute the capacity of the domain G_0
[cap,et{kk},alpha,z2]  = hyppolycap(v,alpha,z2,n{kk});
end
%%
disk = et{1}(1:n{1})./M1;
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
plot(real(et{1}(n{1}+1:2*n{1})),imag(et{1}(n{1}+1:2*n{1})),'--r','LineWidth',1.5)
hold on
plot(real(et{2}(n{2}+1:2*n{2})),imag(et{2}(n{2}+1:2*n{2})),'-.b','LineWidth',1.5)
plot(real(et{3}(n{3}+1:2*n{3})),imag(et{3}(n{3}+1:2*n{3})),'-m','LineWidth',1.5)
plot(real(disk),imag(disk),':k','LineWidth',2)
plot(real(et{1}(1:n{1})),imag(et{1}(1:n{1})),'-k','LineWidth',1.5)
set(gca,'FontSize',20)
% 
legend({'$P_3$','$P_5$','$P_7$','$|z|=\frac{1}{M_1}$'},...
         'Interpreter','LaTeX','location','north','Orientation','horizontal','FontSize',20)
% 
axis equal
axis([-1.4 1.4 -1.03 1.41])
axis off
set(gca,'LooseInset',get(gca,'TightInset'))
% %
print -depsc hyp_poly_fig_ca
%%