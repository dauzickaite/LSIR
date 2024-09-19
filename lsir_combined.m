function results_str = lsir_combined(Adata,ir_parameters,solver,precond,...
    results_str,kappa_ind,noise_ind,x,r,lsqr_parameters)

A = Adata.A;
b = Adata.b;

Rs = Adata.Rs;
ATprec = @(x,tchoice) RTinvAT(x,tchoice,A,Rs);

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

    % compute the residuals for augmented system
     switch ur
        case 'quad'
            f = mp(b) - mp(A)*mp(x) - mp(r);
            g = -mp(A)'*mp(r);
            
            f = double(f);
            g = double(g);
            r = double(r);
         case 'double'
            f = double(b) - double(A)*double(x) - double(r);
            g = -double(A)'*double(r);
            
            f = single(f);
            g = single(g);
            r = single(r);
     end


    % solve
    if precond
        [xn1,~,~,lsit1] = lsqr(A,f,lsqrtol,lsqrit,Rs);
        [xn2,~,~,lsit2] = lsqr(A,r,lsqrtol,lsqrit,Rs);

        RTinvg = (Rs')\g;
        [rn1,~,~,lsit3] = lsqr(ATprec,RTinvg,lsqrtol,lsqrit);
    else
        [xn1,~,~,lsit1] = lsqr(A,f,lsqrtol,lsqrit);
        [xn2,~,~,lsit2] = lsqr(A,r,lsqrtol,lsqrit);

        [rn1,~,~,lsit3] = lsqr(A',g,lsqrtol,lsqrit);
    end
    %rn2 = double(f) - Q*Q'*double(f);
    rn2 = f - A*xn1;

    results_str.inner_iter(kappa_ind,noise_ind) = ...
        results_str.inner_iter(kappa_ind,noise_ind) +...
        lsit1 + lsit2 + lsit3;

    % update, ensure the right precisions
    switch ur
        case 'quad'
            x = double(x) + double(xn1) + double(xn2);
            r = double(r) + double(rn1) + double(rn2);
        case 'double'
            x = single(x) + single(xn1) + single(xn2);
            r = single(r) + single(rn1) + single(rn2);
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
    results_str.inner_iter(kappa_ind,noise_ind) = nan;
end
    
end

function y = RTinvAT(x,tchoice,A,R)
    if strcmp(tchoice,'notransp')
            y = (R')\(A'*x);
    else
            y = A*(R\x);
    end
end
