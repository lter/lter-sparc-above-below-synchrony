% Uses cs_sims to bootstrap nc random communities and bootn replicates
clear all
close all

%rng(1)

np_max=10;
npi=length(2:2:np_max);	% np sampling

nm=10;

nc=10; % 10 random communities
bootn=10;

summ=zeros(3,3,npi,nc,bootn);

tic	
wb=waitbar(0,'Simulating random community #');
k1=1;
for np=2:2:np_max
	n=np+nm;
	for i=1:nc
		intm=(i*0.01)*randn(n,n);
		for j=1:n
			intm(j,j)=1;
		end
		intm=intm*intm';

			for k2=1:bootn
				waitbar((k2+bootn*(i-1))/(nc*bootn),wb);
				summ(:,:,k1,i,k2)=cs_sims(np,nm,intm,0);
			end
		end
	k1=k1+1;
end
toc
close(wb)

save("cs_boot_np_out.mat","summ","npi","nc","bootn");

%figure(1)
%violinplot(reshape(squeeze(summ(3,1,:,:,:)),[npi nc*bootn])');
%xticks(2:2:np_max);
%xlabel('Number of plant species','FontSize',20);
%ylabel('Coupling/Synchrony','FontSize',20);
%
%figure(2)
%violinplot(reshape(squeeze(summ(3,:,5,:,:)),[3 nc*bootn])');

%figure(2)
%plot(reshape(squeeze(summ(1,1,:,:)),[nc*bootn 1]),reshape(squeeze(summ(2,1,:,:)),[nc*bootn 1]),'*');
%hold 'on'
%plot(reshape(squeeze(summ(1,1,:,:)),[nc*bootn 1]),reshape(squeeze(summ(2,2,:,:)),[nc*bootn 1]),'*');
%plot(reshape(squeeze(summ(1,1,:,:)),[nc*bootn 1]),reshape(squeeze(summ(2,3,:,:)),[nc*bootn 1]),'*');
%xlabel("Coupling");
%ylabel("Variance ratio")
%legend({"pm","p","m"});
%
%figure(3)
%plot(reshape(squeeze(summ(1,1,:,:)),[nc*bootn 1]),reshape(squeeze(summ(3,1,:,:)),[nc*bootn 1]),'*');
%xlabel("Coupling initial");
%ylabel("Coupling observed")
