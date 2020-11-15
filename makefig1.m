clear;clc; 
addpath files fmm
%%
t = linspace(0.01,6,1000).';
% f1=@(t)(2*pi./log(1./tanh((1/2).*asinh(t./(2*pi)))));
% f2=@(t)(2*pi./mu(tanh(asinh(sqrt(t./(4*pi))))));
f1=@(t)(2*pi./log(sqrt(1+4*pi^2./t.^2)+2*pi./t));
f2=@(t)(2*pi./mu(tanh(t./4)));
for k=1:length(t)
    y(k,1) = f1(t(k))-f2(t(k));
end
%%
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

hold on; 
plot(t,y,'b','LineWidth',1.5)


axis equal
% axis([-0.4 6.2 -0.8 4.2])
axis([-0.5 6.1 -0.6 2.2])

grid(gca,'minor')
grid on
set(gca, 'XMinorTick','on')
set(gca, 'YMinorTick','on')
ax=gca;
ax.GridAlpha=0.5;
ax.MinorGridAlpha=0.5;

set(gca,'FontSize',20)

ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

text(3,-0.6,{'$c$'},'Interpreter','LaTeX','FontSize',20)
% xlabel('$c$','Interpreter','LaTeX','FontSize',20)
% ylabel('Capacity','Interpreter','LaTeX','FontSize',20)

set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc fig_diff_2
%%