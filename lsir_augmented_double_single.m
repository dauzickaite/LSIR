function lsir_augmented_double_single(m,n,ir_it_max,kappa_list,noise_list,solver,precond,xvalues,yvalues)

kappa_no = length(kappa_list);
noise_no= length(noise_list);

x_error = zeros(kappa_no,noise_no);
r_error = zeros(kappa_no,noise_no);
ir_iter = zeros(kappa_no,noise_no);

switch solver
    case 'iter' 
        gmres_iter = zeros(kappa_no,noise_no);

        lsqrit = n;
        lsqrtol = 1e-7;

        gmresit = 50; 
        gmrestol = 1e-6;
        
        Pright = @(x,prchoice) x;
                
        u = 'single';
        uA = 'single';

case 'QR'
        gmres_iter= [];        
end

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

        % minsvd=min(svd(R));
        alpha = 1; % 2^(-1/2)*minsvd; 

        switch solver
            case 'QR'
                x = R\(Q(:,1:n)'*b);
            case 'iter'
                x = lsqr(A,b,lsqrtol,lsqrit);
               
                Aa = [alpha*eye(m) A; A' zeros(n)];
                P = @(x,prchoice) precAug(x,prchoice,Q,R,alpha);
                
        end

        r = b - A*x;

        u_val = 4*eps('single'); % eps('single')=2u


        x_relerror = norm(mp(x,64) - mp(xtrue,64))/xtruen;
        r_relerror = norm(mp(r,64) - mp(rtrue,64))/rtruen;

        if x_relerror <= u_val %&& r_relerror <= u_val
            converged = true;
        end
         
        ind = 0;
        
        while ~converged && ind < ir_it_max
        ind = ind+1;


            % compute the residuals for augmented system
            f = double(b) - double(A)*double(x) - double(r);
            g = -double(A)'*double(r);

            % solve
            switch solver
                case 'QR'
                    hx = (R')\single(g);
                    k = Q'*single(f);

                    rn = Q*[hx; k(n+1:end)];
                    xn = R\(k(1:n)-hx);

                case 'iter'
                    if precond
                        [z,~,~,git] = mpgmres(Aa,single([alpha*f;g]),...
                            zeros(m+n,1),gmrestol,1,gmresit,Pright,P,u,uA);
                    else
                        Pleft = @(x,prchoice) x;
                        [z,~,~,git] = mpgmres(Aa,single([alpha*f;g]),...
                            zeros(m+n,1),gmrestol,1,gmresit,Pright,Pleft,u,uA);
                        
                    end

                    gmres_iter(kappa_ind,noise_ind) = gmres_iter(kappa_ind,noise_ind) + git;
                    

                    rn = z(1:m);
                    xn = z(m+1:end)/alpha;

            end

            % update
            x = single(x) + single(xn);
            r = single(r) + single(rn);


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
            if strcmp(solver,'iter')
                gmres_iter(kappa_ind,noise_ind) = nan;
            end
        end
    end
end

%% plots
switch solver
    case 'QR'
        plot_convergence_criteria_augmented(kappa_list,noise_list,xvalues,yvalues,2^(24))
end
plot_results(x_error,r_error,ir_iter,xvalues,yvalues,gmres_iter)

end
%% preconditioner
function Px = precAug(x,prchoice,Q,R,alpha)

[n,~] = size(R);
[m,~] = size(Q); 
Q = Q(:,1:n);

x1 = x(1:m);
x2 = x(m+1:end);

Px1 = alpha^(-1)*(x1 - Q*Q'*x1) + Q*(R'\x2);
Px2 = R\(Q'*x1) - alpha*R\(R'\x2);

Px = [Px1; Px2];
end


