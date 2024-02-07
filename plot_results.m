function plot_results(x_error,r_error,ir_iter,xvalues,yvalues)

figure
hx = heatmap(xvalues,yvalues,x_error);

hx.Title = '|| x - x^* || / || x^* ||';
hx.XLabel = '|| r^* ||';
hx.YLabel = '\kappa(A)';
hx.CellLabelFormat = '%.0e';
set(gca,'ColorScaling','log')
set(gca, 'FontSize',50)


%
figure
hr = heatmap(xvalues,yvalues,r_error);
%hr.CellLabelFormat = '%.0e';

hr.Title = '|| r - r^* || / || r^* ||';
hr.XLabel = '|| r^* ||';
hr.YLabel = '\kappa(A)';
hr.CellLabelFormat = '%.0e';
set(gca,'ColorScaling','log')
set(gca, 'FontSize',50)

%
figure
hit = heatmap(xvalues,yvalues,ir_iter);

hit.Title = 'IR iterations';
hit.XLabel = '|| r^* ||';
hit.YLabel = '\kappa(A)';
set(gca, 'FontSize',50)
