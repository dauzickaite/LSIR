% script to obtain all the plots
% requires advanpix toolbox to simulate quadruple and 64 digit precision
clear
addpath '/AdvanpixMCT-4.8.5.14607'
mp.Digits(34);

m = 1e3;
n = 1e1;
ir_it_max = 30;

plot_recog_xtrue_LS_condition(m,n)

%% double + quad precision
kappa_list = [1e0,1e2,1e4,1e6,1e8,1e10,1e12,1e14];
noise_list = [1e0,1e-2,1e-4,1e-6,1e-8,1e-10,1e-12,1e-14];

xvalues = {'1e0','1e-2','1e-4','1e-6','1e-8','1e-10','1e-12','1e-14'};
yvalues = {'1e0','1e2','1e4','1e6','1e8','1e10','1e12','1e14'};


% augmented system
% direct solve via QR
solver = 'QR';
lsir_augmented_quad_double(m,n,ir_it_max,kappa_list,noise_list,solver,[],xvalues,yvalues)

% gmres, no preconditioning
solver = 'iter';
precond = false;
lsir_augmented_quad_double(m,n,ir_it_max,kappa_list,noise_list,solver,precond,xvalues,yvalues)

% gmres, preconditioned
solver = 'iter';
precond = true;
lsir_augmented_quad_double(m,n,ir_it_max,kappa_list,noise_list,solver,precond,xvalues,yvalues)

% LS approach
% direct solve via QR
solver = 'QR';
lsir_LS_quad_double(m,n,ir_it_max,kappa_list,noise_list,solver,[],xvalues,yvalues);

% lsqr, no preconditioning
solver = 'iter';
precond = false;
lsir_LS_quad_double(m,n,ir_it_max,kappa_list,noise_list,solver,precond,xvalues,yvalues)

% lsqr, preconditioned
solver = 'iter';
precond = true;
lsir_LS_quad_double(m,n,ir_it_max,kappa_list,noise_list,solver,precond,xvalues,yvalues)

% seminormal
lsir_seminormal_quad_double(m,n,ir_it_max,kappa_list,noise_list,xvalues,yvalues)

% combined approach
precond = false;
lsir_combined_quad_double(m,n,ir_it_max,kappa_list,noise_list,precond,xvalues,yvalues)

precond = true;
lsir_combined_quad_double(m,n,ir_it_max,kappa_list,noise_list,precond,xvalues,yvalues)


%% single + double precision
kappa_list = [1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7];
noise_list = [1e0,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7];

xvalues = {'1e0','1e-1','1e-2','1e-3','1e-4','1e-5','1e-6','1e-7'};
yvalues = {'1e0','1e1','1e2','1e3','1e4','1e5','1e6','1e7'};


% augmented system
% direct solve via QR
solver = 'QR';
lsir_augmented_double_single(m,n,ir_it_max,kappa_list,noise_list,solver,[],xvalues,yvalues)

% gmres, no preconditioning
solver = 'iter';
precond = false;
lsir_augmented_double_single(m,n,ir_it_max,kappa_list,noise_list,solver,precond,xvalues,yvalues)

% gmres, preconditioned
solver = 'iter';
precond = true;
lsir_augmented_double_single(m,n,ir_it_max,kappa_list,noise_list,solver,precond,xvalues,yvalues)

% LS approach
% direct solve via QR
solver = 'QR';
lsir_LS_double_single(m,n,ir_it_max,kappa_list,noise_list,solver,[],xvalues,yvalues);

% lsqr, no preconditioning
solver = 'iter';
precond = false;
lsir_LS_double_single(m,n,ir_it_max,kappa_list,noise_list,solver,precond,xvalues,yvalues)

% lsqr, preconditioned
solver = 'iter';
precond = true;
lsir_LS_double_single(m,n,ir_it_max,kappa_list,noise_list,solver,precond,xvalues,yvalues)

% seminormal
lsir_seminormal_double_single(m,n,ir_it_max,kappa_list,noise_list,xvalues,yvalues)

% combined approach
precond = false;
lsir_combined_double_single(m,n,ir_it_max,kappa_list,noise_list,precond,xvalues,yvalues)

precond = true;
lsir_combined_double_single(m,n,ir_it_max,kappa_list,noise_list,precond,xvalues,yvalues)
