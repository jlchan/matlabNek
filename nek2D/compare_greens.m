one = @(x) ones(size(x));
zero = @(x) zeros(size(x));

N = 4; % order
K = 16;
CFL = 100;
dt = CFL/(K*N^2);
eps = .01;
c = 1;
% beta = {one,one,zero,zero};
beta = {one,one,one,one};

% clean this part up - consolidate
kx = K; ky = K; Nq = N+1; 
Nqkx = kx*(Nq-1) + 1; Nqky = ky*(Nq-1) + 1; % num global dofs along line

[Mg, Kg, Cg, bcInds] = assemble(N,K,beta); % get global SEM matrices
[X Y] = get_physical_points(N,kx,ky);

A = (1/dt)*Mg + eps*Kg + c*Cg;

middle = round(Nqkx*Nqky/2);
e = zeros(Nqkx*Nqky,1);
e(middle) = .01;
G = A\e;

% agmg inverse
levels = agmg_setup(A);
g_agmg = precond_agmg(levels, e);

figure(1)
surf(X,Y,reshape(G,Nqkx,Nqky));
shading interp
view(2)
cax1 = caxis;
colorbar

figure(2)
surf(X,Y,reshape(g_agmg,Nqkx,Nqky));
shading interp
view(2)
cax2 = caxis;
colorbar

cax = [min([cax1 cax2]) max([cax1 cax2])];
hv = get(0,'Children');
for i = 1:length(hv)
    caxis(cax);
    set(gca,'fontsize',14);    
    xlabel('x','fontsize',14)
    ylabel('y','fontsize',14)
end

% print(figure(1),'-depsc','G_CFL10.eps')
% print(figure(2),'-depsc','agmgG_CFL10.eps')
