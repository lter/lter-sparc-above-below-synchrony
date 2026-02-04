% Testing 3 synchrony metrics (Pearson, Kendall rank, Spearman) with 2 timeseries 

clear all
close all

T=10;
sg=0.1;

% Series 1
dt1=1;
t1=0:dt1:T;

%x1=t1+sg*randn(1,length(t1));
x1=sin(t1)+sg*randn(1,length(t1));

% Series 2
dt2=1;
t2=0:dt2:T;

x2=0.01*t2+sg*randn(1,length(t2));
%x2=0.9*sin(t2)+sg*randn(1,length(t2));

% Plot
plot(t1,x1,'*');
hold 'on'
plot(t2,x2,'*');
hold 'off'

% Syncrony metrics
cm=cov([x1; x2]');
vr=(sum(sum(cm))-sum(diag(cm)))/sum(diag(cm));

cp=corr(x1',x2','Type','Pearson');
ck=corr(x1',x2','Type','Kendall');
cs=corr(x1',x2','Type','Spearman');

fprintf("Variance ratio is %f\n Pearson correlation is %f\n Kendall rank correlation is %f\n Spearman rank correlation is %f\n",vr,cp,ck,cs);
