function results_str = lsir_augmented(Adata,...
    ir_parameters,solver,precond,results_str,kappa_ind,noise_ind,x,fgmres_parameters)


A = Adata.A;
b = Adata.b;

m = Adata.m;
n = Adata.n;

Q = Adata.Q;
R = Adata.R;

Rs = Adata.Rs;

xtrue = Adata.xtrue;
rtrue = Adata.rtrue;
xtruen = Adata.xtruen;
rtruen = Adata.rtruen;

ir_it_max = ir_parameters.ir_it_max;
u_val = ir_parameters.u_val;
ur = ir_parameters.ur;

alpha = 1;

switch solver
    case 'iter' 
        gmresit = fgmres_parameters.it; 
        gmrestol = fgmres_parameters.tol;
                
        u = fgmres_parameters.u;
        uA = fgmres_parameters.uA;
        
        Aa = [alpha*eye(m) A; A' zeros(n)];
        Pr = @(x,prchoice) precAugsketched(x,prchoice,Rs,m);
        Pl = @(x,prchoice) precAugsketched(x,prchoice,Rs',m);

end
         
ind = 0;
converged = false;

r = b - A*x;

while ~converged && ind < ir_it_max
    ind = ind+1;

    % compute the residuals for augmented system
     switch ur
        case 'quad'
            f = mp(b) - mp(A)*mp(x) - mp(r);
            g = -mp(A)'*mp(r);
            rhs = double([alpha*f;g]);
         case 'double'
            f = double(b) - double(A)*double(x) - double(r);
            g = -double(A)'*double(r);
            rhs = single([alpha*f;g]);
     end

    % solve
    switch solver
        case 'QR'
            switch ur
                case 'quad'
                    hx = (R')\double(g);
                    k = Q'*double(f);
                case 'double'
                    hx = (R')\single(g);
                    k = Q'*single(f);
            end

            rn = Q*[hx; k(n+1:end)];
            xn = R\(k(1:n)-hx);

        case 'iter'

            if precond
                  [z,~,~,git] = mpgmres(Aa,rhs,...
                        zeros(m+n,1),gmrestol,1,gmresit,Pr,Pl,u,uA);
            else
                Pright = @(x,prchoice) x;
                Pleft = @(x,prchoice) x;
                [z,~,~,git] = mpgmres(Aa,rhs,...
                    zeros(m+n,1),gmrestol,1,gmresit,Pright,Pleft,u,uA);

            end

            results_str.inner_iter(kappa_ind,noise_ind) = ...
                results_str.inner_iter(kappa_ind,noise_ind) + git;

            rn = z(1:m);
            xn = z(m+1:end)/alpha;

    end

    % update, ensure the right precisions
    switch ur
        case 'quad'
            x = double(x) + double(xn);
            r = double(r) + double(rn);
        case 'double'
            x = single(x) + single(xn);
            r = single(r) + single(rn);
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

%% sketched preconditioner
function Px = precAugsketched(x,prchoice,R,m)

Px1 = x(1:m);
Px2 = R\x(m+1:end);

Px = [Px1; Px2];
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


