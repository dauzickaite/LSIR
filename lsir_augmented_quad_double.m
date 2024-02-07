function lsir_augmented_quad_double(m,n,ir_it_max,kappa_list,noise_list,solver,precond,xvalues,yvalues)

kappa_no = length(kappa_list);
noise_no= length(noise_list);

x_error = zeros(kappa_no,noise_no);
r_error = zeros(kappa_no,noise_no);
ir_iter = zeros(kappa_no,noise_no);

switch solver
    case 'iter' 
        gmres_iter = zeros(kappa_no,noise_no);
        gmres_flag = zeros(kappa_no,noise_no);

        lsqrit = n;
        lsqrtol = 1e-14;

        gmresit = n+m;
        gmrestol = 1e-14;
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
        [Qn,~] = qr([mp(A,64),mp(randn(m,1),64)],0);
        rt =  noise.*Qn(:,end);
        b = double(mp(A,64)*mp(xt,64) + mp(rt,64));

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
                P = @(x) precAug(x,Q,R,alpha);
        end

        r = b - A*x;

        u_val = 4*eps('double'); % eps('double')=2u


        for ind = 1:ir_it_max

            % compute the residuals for augmented system
            f = mp(b) - mp(A)*mp(x) - mp(r);
            g = -mp(A)'*mp(r);

            % solve
            switch solver
                case 'QR'
                    hx = (R')\double(g);
                    k = Q'*double(f);

                    rn = Q*[hx; k(n+1:end)];
                    xn = R\(k(1:n)-hx);

                case 'iter'
                    if precond
                        [z,flag,~,git] = gmres(Aa,double([alpha*f;g]),gmresit,gmrestol,[],P);
                    else
                        [z,flag,~,git] = gmres(Aa,double([alpha*f;g]),gmresit,gmrestol);
                    end

                    gmres_iter(kappa_ind,noise_ind) = git(2);
                    gmres_flag(kappa_ind,noise_ind) = flag;

                    rn = z(1:m);
                    xn = z(m+1:end)/alpha;

            end

            % update
            x = x + xn;
            r = r + rn;


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
switch solver
    case 'QR'
        plot_convergence_criteria_augmented(kappa_list,noise_list,xvalues,yvalues,2^(53))
end
plot_results(x_error,r_error,ir_iter,xvalues,yvalues)

end
%% preconditioner
function Px = precAug(x,Q,R,alpha)

[n,~] = size(R);
[m,~] = size(Q); 
Q = Q(:,1:n);

x1 = x(1:m);
x2 = x(m+1:end);

Px1 = alpha^(-1)*(x1 - Q*Q'*x1) + Q*(R'\x2);
Px2 = R\(Q'*x1) - alpha*R\(R'\x2);

Px = [Px1; Px2];
end


