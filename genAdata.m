function [A,Q,R,nv,Rs] = genAdata(m,n,kappa,solver,s,seed,ir_parameters)

rng(seed)
A = gallery('randsvd',[m,n],kappa,3);
if strcmp(ir_parameters.u,'single')
    A = single(A);
    vect_tol = 1e-6;
else
    vect_tol = 1e-15;
end

[Qn,~] = qr([mp(A,64),mp(randn(m,1),64)],0);
nv = Qn(:,end);


switch solver
    case 'QR'
        Rs = [];
        
        [Q,R] = qr(A);
        R = R(1:n,1:n);
                
    case 'iter'
        Q = [];
        R = [];
        
        d = s*n;
        G = (1/sqrt(d))*randn(m,d);
        [Qs,Rs] = qr(G'*A);
        Rs = Rs(1:n,1:n);
        
%         nv = G*Qs(:,end); % vector in the nullspace of A   
%         normATnv = norm(A'*nv);
% 
%         if normATnv <= vect_tol
%                 done = true;
%             else
%                 tries = 0;
%                 done = false;
%         end   
% 
%         while ~done
% 
%             tries = tries +1;
%             G = (1/sqrt(d))*randn(m,d);
%             [Qs,~] = qr(G'*[A,randn(m,1)],0);
% 
%             nv = G*Qs(:,end); % vector in the nullspace of A   
%             normATnv = norm(A'*nv);
% 
%             if normATnv <= vect_tol
%                 done = true;
%             elseif tries <= 3
%                 d = d + n;
%             else
%                 error('Cannot obtain a vector in a null space')
%             end   
%         end


end








