function out=cs_sims(np,nm,intm,viz)

% Testing syncrony in simulated communities where time series are drawn from multivariate normal

%viz: 0 no plot, 1 plot
%int_m interaction matrix - affects correlation of multivariate normal random sampling 

%	close all

	%rng(1)
	%tic()
	n=np+nm;

	Y=10;

	tp=0:1:Y;
	tm=0:1/400:Y;

	xp=zeros(np,length(tp));
	xm=zeros(nm,length(tm));

	if isempty(intm)
		intm=0.01*randn(n,n);
		for j=1:n
			intm(j,j)=1;
		end
		intm=intm*intm';
	end

	%assuming we don't see changes in plants intrannually but they are still sampled intraannually to show effect on microbe

	%samp=mvnrnd(ones(n,1),intm,1)

	tpi=1;

	ym=zeros(nm,length(tp));

	for i=2:length(tm)
		samp=mvnrnd(zeros(n,1),intm,1);
		xm(:,i)=xm(:,i-1)+(1/20)*samp(np+1:end)';
		if(tm(i)==tp(tpi+1))
			tpi=tpi+1;
			xp(:,tpi)=xp(:,tpi-1)+samp(1:np)';
			ym(:,tpi)=xm(:,i);
		end
	end

	if viz==1
		figure(1)
		p1=plot(tp,xp(1,:),'.-','Color','#0072BD','MarkerSize',15);
		hold 'on';
		plot(tp,xp(2:end,:),'.-','Color','#0072BD','MarkerSize',15);
		p2=plot(tm,xm(1,:),'-','Color','#D95319');
		plot(tm,xm(2:end,:),'-','Color','#D95319');
		xlabel('Time (Years)','FontSize',20);
		ylabel('Change in density','FontSize',20);
		legend([p1 p2],{'Plants','Soil microbes'},'Location','northoutside','Orientation','horizontal')
	end

	cov_pm=cov([xp; ym]');
	vr_pm=(sum(sum(cov_pm))-sum(diag(cov_pm)))/sum(diag(cov_pm));
	cpl_pm=mean(mean(intm));
	r_cpl_pm=mean(mean(cov_pm));

	cov_p=cov(xp');
	vr_p=(sum(sum(cov_p))-sum(diag(cov_p)))/sum(diag(cov_p));
	cpl_p=mean(mean(abs(intm(1:np,1:np))));
	r_cpl_p=mean(mean(cov_p));

	cov_m=cov(xm');
	vr_m=(sum(sum(cov_m))-sum(diag(cov_m)))/sum(diag(cov_m));
	cpl_m=mean(mean(abs(intm(np+1:end,np+1:end))));
	r_cpl_m=mean(mean(cov_m));

	%toc()

	out=[cpl_pm cpl_p cpl_m; vr_pm vr_p vr_m; r_cpl_pm r_cpl_p r_cpl_m];

end
