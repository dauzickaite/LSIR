function plot_convergence_criteria_augmented(kappa_list,noise_list,xvalues,yvalues,udinv)

kappa_list_sq = kappa_list.^2;

nk = kappa_list_sq'.*noise_list;
convtbl = nk < udinv;

xg = linspace(0.5,8.5,9);
yg = xg;

figure
imagesc(convtbl); hold on
hm = mesh(xg,yg,zeros(9));
hm.FaceColor = 'none';
hm.EdgeColor = 'k';

%imagesc(convtbl)
%axis ij
colormap('gray')
xticklabels(xvalues)
yticklabels(yvalues)
xlabel('|| r^* ||');
ylabel('\kappa(A)');
set(gca, 'FontSize',50)

