function lsir_seminormal_quad_double(m,n,ir_it_max,kappa_list,noise_list,xvalues,yvalues)

kappa_no = length(kappa_list);
noise_no= length(noise_list);

x_error = zeros(kappa_no,noise_no);
r_error = zeros(kappa_no,noise_no);
ir_iter = zeros(kappa_no,noise_no);


rng(1234)
xt = rand(n,1);
xt = xt/norm(xt);


for kappa_ind = 1:kappa_no
    for noise_ind = 1:noise_no
        
        converged = false;

        kappa = kappa_list(kappa_ind);
        noise = noise_list(noise_ind);


        rng(2023)

        A = gallery('randsvd',[m,n],kappa,3);
        [Qn,~] = qr([mp(A,64),mp(randn(m,1),64)],0);
        rt =  noise.*Qn(:,end);
        b = double(mp(A,64)*mp(xt,64) + mp(rt,64));

        xtrue = mp(A,64)\mp(b,64);
        xtruen = norm(xtrue);
        rtrue = mp(b,64) - mp(A,64)*mp(xtrue);
        rtruen = norm(mp(rtrue,64));


        [Q,R] = qr(A);
        R = R(1:n,1:n);
        Q = Q(:,1:n);

        x = R\(Q'*b);

        u_val = 4*eps('double'); % eps('double')=2u

         % compute the residual
         r = mp(b) - mp(A)*mp(x);


        x_relerror = norm(mp(x,64) - mp(xtrue,64))/xtruen;
        r_relerror = norm(mp(r,64) - mp(rtrue,64))/rtruen;

        if x_relerror <= u_val %&& r_relerror <= u_val
            converged = true;
        end
         
        ind = 0;
        
        while ~converged && ind < ir_it_max
        ind = ind+1;
            % solve
            ATr = mp(A')*r;
            xn1 = (R')\double(ATr);
            xn = R\xn1;

            % update
            x = x + xn;
            
            % compute the residual
            r = mp(b) - mp(A)*mp(x);


           r_relerror = norm(mp(r,64) - mp(rtrue,64))/rtruen;
           x_relerror = norm(mp(x,64) - mp(xtrue,64))/xtruen;

            if x_relerror <= u_val %&& r_relerror <= u_val
                converged = true;
            end

        end

        x_error(kappa_ind,noise_ind) = x_relerror;
        r_error(kappa_ind,noise_ind) = r_relerror;
        
        

        if converged
            ir_iter(kappa_ind,noise_ind) = ind;
        else
            ir_iter(kappa_ind,noise_ind) = nan;
        end
    end
end

%% plots

plot_results(x_error,r_error,ir_iter,xvalues,yvalues,[])

end


