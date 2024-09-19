function plot_results(results_str,xvalues,yvalues,plot_inner)

figure
hx = heatmap(xvalues,yvalues,results_str.x_error);

hx.Title = '|| x - x^* || / || x^* ||';
hx.XLabel = '|| r^* ||';
hx.YLabel = '\kappa(A)';
hx.CellLabelFormat = '%.0e';
set(gca,'ColorScaling','log')
set(gca, 'FontSize',50)


%
figure
hr = heatmap(xvalues,yvalues,results_str.r_error);
%hr.CellLabelFormat = '%.0e';

hr.Title = '|| r - r^* || / || r^* ||';
hr.XLabel = '|| r^* ||';
hr.YLabel = '\kappa(A)';
hr.CellLabelFormat = '%.0e';
set(gca,'ColorScaling','log')
set(gca, 'FontSize',50)

%
figure
hit = heatmap(xvalues,yvalues,results_str.ir_iter);

hit.Title = 'IR iterations';
hit.XLabel = '|| r^* ||';
hit.YLabel = '\kappa(A)';
set(gca, 'FontSize',50)

if plot_inner
    figure
    hit = heatmap(xvalues,yvalues,results_str.inner_iter);

    hit.Title = 'Inner iterations';
    hit.XLabel = '|| r^* ||';
    hit.YLabel = '\kappa(A)';
    set(gca, 'FontSize',50)
end


