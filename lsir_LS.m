function results_str = lsir_LS(Adata,ir_parameters,solver,precond,...
    results_str,kappa_ind,noise_ind,x,r,lsqr_parameters)

A = Adata.A;
b = Adata.b;

n = Adata.n;

Q = Adata.Q;
switch solver
    case 'QR'
        Q = Q(:,1:n);
end
R = Adata.R;

Rs = Adata.Rs;

xtrue = Adata.xtrue;
rtrue = Adata.rtrue;
xtruen = Adata.xtruen;
rtruen = Adata.rtruen;

ir_it_max = ir_parameters.ir_it_max;
u_val = ir_parameters.u_val;
ur = ir_parameters.ur;


switch solver
    case 'iter' 
        lsqrit = lsqr_parameters.it;
        lsqrtol = lsqr_parameters.tol;               
end

ind = 0;
converged = false;

        
while ~converged && ind < ir_it_max
    ind = ind+1;

    switch ur
        case 'quad'
            rl = double(r);
        case 'double'
            rl = single(r);
    end

    
    % solve
    switch solver
        case 'QR'
             QTr = Q'*rl;
             xn = R\QTr;
        case 'iter'
            if precond
                [xn,~,~,lsit] = lsqr(A,rl,lsqrtol,lsqrit,Rs);
            else
                [xn,~,~,lsit] = lsqr(A,rl,lsqrtol,lsqrit);
            end

            results_str.inner_iter(kappa_ind,noise_ind) = ...
                results_str.inner_iter(kappa_ind,noise_ind) + lsit;

    end


    % update and compute the residual
    switch ur
        case 'quad'
            x = double(x) + double(xn);
            r = mp(b) - mp(A)*mp(x);
        case 'double'
            x = single(x) + single(xn);
            r = double(b) - double(A)*double(x);
    end


   r_relerror = norm(mp(r,64) - mp(rtrue,64))/rtruen;
   x_relerror = norm(mp(x,64) - mp(xtrue,64))/xtruen;
   
   results_str.x_conv{kappa_ind,noise_ind}(ind+1) = x_relerror;
   results_str.r_conv{kappa_ind,noise_ind}(ind+1) = r_relerror;


    if x_relerror <= u_val %&& r_relerror <= u_val
        converged = true;
    end

end

results_str.x_error(kappa_ind,noise_ind) = x_relerror;
results_str.r_error(kappa_ind,noise_ind) = r_relerror;

if converged
    results_str.ir_iter(kappa_ind,noise_ind) = ind;
else
    results_str.ir_iter(kappa_ind,noise_ind) = nan;
    if strcmp(solver,'iter')
        results_str.inner_iter(kappa_ind,noise_ind) = nan;
    end
end

end


