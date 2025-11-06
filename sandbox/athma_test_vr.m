clear all
close all

rng(1)

dt=0.1; % resolution
T=100;	% length

N=10; %number of variables

t=0:dt:T;

nt=length(t);

ns1=0.5;
ns2=1.5;

x=sin(0.1*t+pi*ns2*rand(N,1))+ns1*rand(N,nt);
%x=[1e-3; 1e-2]*t+ns1*rand(2,nt);

cm=cov(x');
corr(x')

vr=(sum(sum(cm))-sum(diag(cm)))/sum(diag(cm))

y=x(:,1:10:nt);

cm=cov(y');
vr=(sum(sum(cm))-sum(diag(cm)))/sum(diag(cm))

figure(1)
plot(t,x,'*-')

figure(2)
plot(t(1:10:nt),y,'*-')
