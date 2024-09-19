% solve problems via iterative methods
clear
addpath '../AdvanpixMCT-4.8.5.14607'
mp.Digits(34);

m = 1e3;
n = 1e1;

%% double + quad precision
kappa_list = [1e0,1e2,1e4,1e6,1e8,1e10,1e12,1e14];
noise_list = [1e0,1e-2,1e-4,1e-6,1e-8,1e-10,1e-12,1e-14];

ir_parameters = struct('ir_it_max', 30, 'u_val', 4*eps('double'),...
    'u', 'double', 'ur', 'quad');

lsqr_parameters = struct('tol',1e-14,'it',n);

fgmres_parameters = struct('tol',1e-12,'it',50,'u','double','uA','double');

kappa_no = length(kappa_list);
noise_no = length(noise_list);

% create result data structures
solver_no = 3;

results(1:solver_no) = struct('x_error',zeros(kappa_no,noise_no),...
    'r_error',zeros(kappa_no,noise_no),'ir_iter',zeros(kappa_no,noise_no),...
    'inner_iter',zeros(kappa_no,noise_no));


Adata = struct();
Adata.m = m;
Adata.n = n;


% create the problem data 
rng(1234)
xt = rand(n,1);
xt = xt/norm(xt);

s = 4;


seed = 2023;
for kappa_ind = 1:kappa_no
    
    kappa = kappa_list(kappa_ind);
    
    [Adata.A,Adata.Q,Adata.R,nv,Adata.Rs] = genAdata(m,n,kappa,'iter',...
        s,seed,ir_parameters);
    
     % begin loop for rhs
     for noise_ind = 1:noise_no

            converged = false;

            noise = noise_list(noise_ind);

            rt =  noise.*nv;
            Adata.b = double(mp(Adata.A,64)*mp(xt,64) + mp(rt,64));
            
            Adata.xtrue = mp(Adata.A,64)\mp(Adata.b,64);
            Adata.xtruen = norm(Adata.xtrue);
            Adata.rtrue = mp(Adata.b,64) - mp(Adata.A,64)*mp(Adata.xtrue);
            Adata.rtruen = norm(mp(Adata.rtrue,64));  
           
            % x_0
            x = lsqr(Adata.A,Adata.b,lsqr_parameters.tol,...
                lsqr_parameters.it,Adata.Rs);
           

           % r = Adata.b - Adata.A*x;
            r = mp(Adata.b) - mp(Adata.A)*mp(x);

            x_relerror = norm(mp(x,64) - mp(Adata.xtrue,64))/Adata.xtruen;
            r_relerror = norm(mp(r,64) - mp(Adata.rtrue,64))/Adata.rtruen;

            for sind = 1:solver_no
                   results(sind).x_conv{kappa_ind,noise_ind}(1) = x_relerror;
                   results(sind).r_conv{kappa_ind,noise_ind}(1) = r_relerror;
            end

            if x_relerror <= ir_parameters.u_val %&& r_relerror <= u_val
                converged = true;
            end

            % if not converged then solve with all the methods
            if converged
                % save the x and r errors
                for sind = 1:solver_no
                    results(sind).x_error(kappa_ind,noise_ind) = x_relerror;
                    results(sind).r_error(kappa_ind,noise_ind) = r_relerror;
                end
            else
                results(1) = lsir_augmented(Adata,ir_parameters,'iter',...
                   true,results(1),kappa_ind,noise_ind,x,fgmres_parameters);
                
                
                 results(2) = lsir_LS(Adata,ir_parameters,'iter',...
                     true,results(2),kappa_ind,noise_ind,x,r,lsqr_parameters);
                
                 results(3) = lsir_combined(Adata,ir_parameters,'iter',...
                     true,results(3),kappa_ind,noise_ind,x,r,lsqr_parameters);
            end

    end

end
% plot the results
xvalues = {'1e0','1e-2','1e-4','1e-6','1e-8','1e-10','1e-12','1e-14'};
yvalues = {'1e0','1e2','1e4','1e6','1e8','1e10','1e12','1e14'};

plot_results(results(1),xvalues,yvalues,true)
plot_results(results(2),xvalues,yvalues,true)
plot_results(results(3),xvalues,yvalues,true)


%% single + double precision
kappa_list = [1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7];
noise_list = [1e0,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7];

ir_parameters = struct('ir_it_max', 30, 'u_val', 4*eps('single'),...
    'u', 'single', 'ur', 'double');

lsqr_parameters = struct('tol',1e-7,'it',n);
fgmres_parameters = struct('tol',1e-6,'it',50,'u','single','uA','single');

kappa_no = length(kappa_list);
noise_no = length(noise_list);

% create result data structures
solver_no = 3;

results_us(1:solver_no) = struct('x_error',zeros(kappa_no,noise_no),...
    'r_error',zeros(kappa_no,noise_no),'ir_iter',zeros(kappa_no,noise_no),...
    'inner_iter',zeros(kappa_no,noise_no));


Adatas = struct();
Adatas.m = m;
Adatas.n = n;


% create the problem data 
rng(1234)
xt = rand(n,1);
xt = xt/norm(xt);

seed = 2023;
for kappa_ind = 1:kappa_no
    
    kappa = kappa_list(kappa_ind);

    [Adatas.A,Adatas.Q,Adatas.R,nvs,Adatas.Rs] = genAdata(m,n,kappa,'iter',...
        s,seed,ir_parameters);
    
     % begin loop for rhs
     for noise_ind = 1:noise_no

            converged = false;

            noise = noise_list(noise_ind);

            rt =  noise.*nvs;
            Adatas.b = single(mp(Adatas.A,64)*mp(xt,64) + mp(rt,64));
            
            Adatas.xtrue = mp(Adatas.A,64)\mp(Adatas.b,64);
            Adatas.xtruen = norm(Adatas.xtrue);
            Adatas.rtrue = mp(Adatas.b,64) - mp(Adatas.A,64)*mp(Adatas.xtrue);
            Adatas.rtruen = norm(mp(Adatas.rtrue,64));        


            % x_0
            x = lsqr(Adatas.A,Adatas.b,lsqr_parameters.tol,...
                lsqr_parameters.it,Adatas.Rs);

           % r = Adata.b - Adata.A*x;
            r = double(Adatas.b) - double(Adatas.A)*double(x);

            x_relerror = norm(mp(x,64) - mp(Adatas.xtrue,64))/Adatas.xtruen;
            r_relerror = norm(mp(r,64) - mp(Adatas.rtrue,64))/Adatas.rtruen;

            for sind = 1:solver_no
                   results_us(sind).x_conv{kappa_ind,noise_ind}(1) = x_relerror;
                   results_us(sind).r_conv{kappa_ind,noise_ind}(1) = r_relerror;
            end

            if x_relerror <= ir_parameters.u_val %&& r_relerror <= u_val
                converged = true;
            end

            % if not converged then solve with all the methods
            if converged
                % save the x and r errors
                for sind = 1:solver_no
                    results_us(sind).x_error(kappa_ind,noise_ind) = x_relerror;
                    results_us(sind).r_error(kappa_ind,noise_ind) = r_relerror;
                end
            else
                results_us(1) = lsir_augmented(Adatas,ir_parameters,'iter',...
                   true,results_us(1),kappa_ind,noise_ind,x,fgmres_parameters);
                
                
                 results_us(2) = lsir_LS(Adatas,ir_parameters,'iter',...
                     true,results_us(2),kappa_ind,noise_ind,x,r,lsqr_parameters);
                
                results_us(3) = lsir_combined(Adatas,ir_parameters,'iter',...
                     true,results_us(3),kappa_ind,noise_ind,x,r,lsqr_parameters);
            end

    end

end
% plot the results
xvalues = {'1e0','1e-1','1e-2','1e-3','1e-4','1e-5','1e-6','1e-7'};
yvalues = {'1e0','1e1','1e2','1e3','1e4','1e5','1e6','1e7'};

plot_results(results_us(1),xvalues,yvalues,true)
plot_results(results_us(2),xvalues,yvalues,true)
plot_results(results_us(3),xvalues,yvalues,true)

% plot x convergence
plot_xconverg(results_us,2,1,true,{'augmented','LS','combined'})
plot_xconverg(results_us,7,1,false,[])
plot_xconverg(results_us,2,7,false,[])
plot_xconverg(results_us,7,7,false,[])

% plot r convergence
plot_rconverg(results_us,2,1,true,{'augmented','LS','combined'})
plot_rconverg(results_us,7,1,false,[])
plot_rconverg(results_us,2,7,false,[])
plot_rconverg(results_us,7,7,false,[])

