function plot_recog_xtrue_LS_condition(m,n)

kappa_list = [1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12,1e13,1e14];
noise_list = [1e0,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14];

kappa_no = length(kappa_list);
noise_no= length(noise_list);


rng(1234)
xt = rand(n,1);
xt = xt/norm(xt);

cond_ls_LHS= zeros(kappa_no,noise_no);
cond_ls_RHS= zeros(kappa_no,noise_no);


for kappa_ind = 1:kappa_no
    for noise_ind = 1:noise_no


    kappa = kappa_list(kappa_ind);
    noise = noise_list(noise_ind);


    rng(2023)

    
    A = gallery('randsvd',[m,n],kappa,3); % norm(A) = 1;
    [Qn,~] = qr(mp([mp(A),mp(randn(m,1))]),0);
    rt =  noise.*Qn(:,end);
    b = double(mp(A,64)*mp(xt,64) + mp(rt,64));

    xtrue = mp(A,64)\mp(b,64);
    xtruen = norm(xtrue);
    rtruen = noise;
    
    cond_ls_LHS(kappa_ind,noise_ind) = rtruen/xtruen;
    cond_ls_RHS(kappa_ind,noise_ind) = kappa^(-2);
    

    end
end

cond_sat = cond_ls_RHS > cond_ls_LHS;


xvalues = {'1e0','1e-2','1e-4','1e-6','1e-8','1e-10','1e-12','1e-14'};
yvalues = {'1e0','1e2','1e4','1e6','1e8','1e10','1e12','1e14'};


xg = linspace(0.5,15.5,16);
yg = xg;
figure
imagesc(cond_sat); hold on
hm = mesh(xg,yg,zeros(16));
hm.FaceColor = 'none';
hm.EdgeColor = 'k';

%imagesc(convtbl)
%axis ij
colormap('gray')
xticks(1:2:15)
yticks(1:2:15)
xticklabels(xvalues)
yticklabels(yvalues)
xlabel('|| r^* ||');
ylabel('\kappa(A)');
set(gca, 'FontSize',50)

