function lsir_combined_quad_double(m,n,ir_it_max,kappa_list,noise_list,precond,xvalues,yvalues)

kappa_no = length(kappa_list);
noise_no= length(noise_list);

x_error = zeros(kappa_no,noise_no);
r_error = zeros(kappa_no,noise_no);
ir_iter = zeros(kappa_no,noise_no);


lsqrit = n;
lsqrtol = 1e-14;


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


        x = lsqr(A,b,lsqrtol,lsqrit);
        r = b - A*x;

        u_val = 4*eps('double'); % eps('double')=2u
        
        ATprec = @(x,tchoice) RTinvAT(x,tchoice,A,R);

        for ind = 1:ir_it_max

            % compute the residuals of the augmented system               
            f = mp(b) - mp(A)*mp(x) - mp(r);
            g = -mp(A)'*mp(r);

            % solve
            if precond
                xn1 = lsqr(A,double(f),lsqrtol,lsqrit,R);
                xn2 = lsqr(A,double(r),lsqrtol,lsqrit,R);
                
                RTinvg = (R')\double(g);
                rn1 = lsqr(ATprec,RTinvg,lsqrtol,lsqrit);
            else
                xn1 = lsqr(A,double(f),lsqrtol,lsqrit);
                xn2 = lsqr(A,double(r),lsqrtol,lsqrit);

                rn1 = lsqr(A',double(g),lsqrtol,lsqrit);
            end
            rn2 = double(f) - Q*Q'*double(f);
            

            % update
            x = x + xn1 + xn2;
            r = r + rn1 + rn2;


           r_relerror = norm(mp(r,64) - mp(rtrue,64))/rtruen;
           x_relerror = norm(mp(x,64) - mp(xtrue,64))/xtruen;

            if x_relerror <= u_val %&& r_relerror <= u_val
                converged = true;
                break
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

plot_results(x_error,r_error,ir_iter,xvalues,yvalues)

end
function y = RTinvAT(x,tchoice,A,R)
    if strcmp(tchoice,'notransp')
            y = (R')\(A'*x);
    else
            y = A*(R\x);
    end
end
