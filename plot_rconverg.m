function plot_rconverg(results_str,k,r,legendOn,legtxt)

figure; semilogy(0:length(results_str(1).r_conv{k,r})-1,results_str(1).r_conv{k,r},'LineWidth',8); hold on
semilogy(0:length(results_str(2).r_conv{k,r})-1,results_str(2).r_conv{k,r},'--','LineWidth',8); hold on
semilogy(0:length(results_str(3).r_conv{k,r})-1,results_str(3).r_conv{k,r},'k:','LineWidth',8); hold on
if legendOn
     legend(legtxt)
end
xlabel('IR iteration');
ylabel('|| r_i - r^* || / || r^* ||');
ylim([1e-9 1e4])
yticks([1e-8 1e-4 1e0 1e4])
yticklabels({'1e-8','1e-4','1e0','1e4'})
set(gca, 'FontSize',50)