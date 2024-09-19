function results_str = lsir_seminormal(Adata,ir_parameters,...
    results_str,kappa_ind,noise_ind,x,r)

A = Adata.A;
b = Adata.b;

R = Adata.R;

xtrue = Adata.xtrue;
rtrue = Adata.rtrue;
xtruen = Adata.xtruen;
rtruen = Adata.rtruen;

ir_it_max = ir_parameters.ir_it_max;
u_val = ir_parameters.u_val;
ur = ir_parameters.ur;

ind = 0;
converged = false;

while ~converged && ind < ir_it_max
    ind = ind+1;
        % solve
        switch ur
            case 'quad'
                ATr = mp(A')*r;
                xn1 = (R')\double(ATr);
            case 'double'
                ATr = double(A')*r;
                xn1 = (R')\single(ATr);
        end

        xn = R\xn1;

        % update
        x = x + xn;

        % compute the residual
        switch ur
            case 'quad'
                r = mp(b) - mp(A)*mp(x);
            case 'double'
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
    end
    
end


