function lsir_combined_double_single(m,n,ir_it_max,kappa_list,noise_list,precond,xvalues,yvalues)

kappa_no = length(kappa_list);
noise_no= length(noise_list);

x_error = zeros(kappa_no,noise_no);
r_error = zeros(kappa_no,noise_no);
ir_iter = zeros(kappa_no,noise_no);

inner_it = zeros(kappa_no,noise_no);

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
        A = single(A);
        [Qn,~] = qr([mp(A,64),mp(randn(m,1),64)],0);
        rt =  noise.*Qn(:,end);
        b = single(mp(A,64)*mp(xt,64) + mp(rt,64));

        xtrue = mp(A,64)\mp(b,64);
        xtruen = norm(xtrue);
        rtrue = mp(b,64) - mp(A,64)*mp(xtrue);
        rtruen = norm(mp(rtrue,64));


        [Q,R] = qr(A);
        R = R(1:n,1:n);
        Q = Q(:,1:n);


        x = lsqr(A,b,lsqrtol,lsqrit);
        r = b - A*x;

        u_val = 4*eps('single'); % eps('double')=2u
        
        ATprec = @(x,tchoice) RTinvAT(x,tchoice,A,R);

        x_relerror = norm(mp(x,64) - mp(xtrue,64))/xtruen;
        r_relerror = norm(mp(r,64) - mp(rtrue,64))/rtruen;
        
        if x_relerror <= u_val %&& r_relerror <= u_val
            converged = true;
        end
         
        ind = 0;
        
        while ~converged && ind < ir_it_max
        ind = ind+1;

            % compute the residuals of the augmented system               
            f = double(b) - double(A)*double(x) - double(r);
            g = -double(A)'*double(r);

            % solve
            if precond
                [xn1,~,~,lsit1] = lsqr(A,single(f),lsqrtol,lsqrit,R);
                [xn2,~,~,lsit2] = lsqr(A,single(r),lsqrtol,lsqrit,R);
                
                RTinvg = (R')\single(g);
                [rn1,~,~,lsit3] = lsqr(ATprec,RTinvg,lsqrtol,lsqrit);
            else
                [xn1,~,~,lsit1] = lsqr(A,single(f),lsqrtol,lsqrit);
                [xn2,~,~,lsit2] = lsqr(A,single(r),lsqrtol,lsqrit);

                [rn1,~,~,lsit3] = lsqr(A',single(g),lsqrtol,lsqrit);
            end
            rn2 = single(f) - Q*Q'*single(f);
            

            inner_it(kappa_ind,noise_ind) = inner_it(kappa_ind,noise_ind) +...
                lsit1 + lsit2 + lsit3;
            
            % update
            x = x + xn1 + xn2;
            r = r + rn1 + rn2;


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
            inner_it(kappa_ind,noise_ind) = nan;
        end
    end
end

%% plots

plot_results(x_error,r_error,ir_iter,xvalues,yvalues,inner_it)

end
function y = RTinvAT(x,tchoice,A,R)
    if strcmp(tchoice,'notransp')
            y = (R')\(A'*x);
    else
            y = A*(R\x);
    end
end
