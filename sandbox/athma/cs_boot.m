% Uses cs_sims to bootstrap nc random communities and bootn replicates
clear all
close all

rng(1)

np=15;
nm=100;
n=np+nm;

nc=20; % 10 random communities
bootn=10;

summ=zeros(2,3,nc,bootn);

tic		
for i=1:nc
	intm=(i*0.01)*randn(n,n);
	for j=1:n
		intm(j,j)=1;
	end
	intm=intm*intm';

	for k=1:bootn
		if i==1 && k==bootn
			figure(1)
			summ(:,:,i,k)=cs_sims(np,nm,intm,1);
		else
			summ(:,:,i,k)=cs_sims(np,nm,intm,0);
		end
	end
	i
end
toc

figure(2)
plot(reshape(squeeze(summ(1,1,:,:)),[nc*bootn 1]),reshape(squeeze(summ(2,1,:,:)),[nc*bootn 1]),'*');
hold 'on'
plot(reshape(squeeze(summ(1,1,:,:)),[nc*bootn 1]),reshape(squeeze(summ(2,2,:,:)),[nc*bootn 1]),'*');
plot(reshape(squeeze(summ(1,1,:,:)),[nc*bootn 1]),reshape(squeeze(summ(2,3,:,:)),[nc*bootn 1]),'*');
xlabel("Coupling");
ylabel("Variance ratio")
legend({"pm","p","m"});
