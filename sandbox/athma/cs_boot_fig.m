clear all
close all

load("cs_boot_np_out.mat");

figure(1)
violinplot(reshape(squeeze(summ(3,1,:,:,:)),[npi nc*bootn])');
xticklabels({'2','4','6','8','10'});
xlabel('Number of plant species','FontSize',15);
ylabel('Coupling (Ochoa-Hueso et al. 2021)','FontSize',15);

figure(2)
violinplot(reshape(squeeze(summ(3,:,5,:,:)),[3 nc*bootn])');
