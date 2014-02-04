addpath('./agmg')
one = @(x) ones(size(x));
zero = @(x) zeros(size(x));

N = 16; % order
K = 2;
CFL = 20;
dt = CFL/(K*N^2);
dt = 1.0;
eps = .001;
c = 1;
% beta = {one,one,zero,zero};
%beta = {one,one,zero,one};
beta = {zero,zero,one,one}; % for automatic "downwind ordering"

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

% wathen inverse
PGS = tril((A));
start = 0;
for k = 1:Nqky
    inds = start + [1:Nqkx];    
    PGS(inds,inds) = A(inds,inds); %grab upper blocks
    start = start + Nqkx;
end
g_wathen = PGS\e; %block Gauss seidel matrix

cax = [];
figure(1)
surf(X,Y,reshape(G,Nqkx,Nqky));
shading interp
cax = [cax caxis];
colorbar

figure(2)
surf(X,Y,reshape(g_agmg,Nqkx,Nqky));
shading interp
cax = [cax caxis];
colorbar

figure(3)
surf(X,Y,reshape(g_wathen,Nqkx,Nqky));
shading interp
cax = [cax caxis];
colorbar

% cax = [min(cax) max(cax)];
% hv = get(0,'Children');
% for i = 1:length(hv)
%     caxis(cax);
%     set(gca,'fontsize',14);    
%     xlabel('x','fontsize',14)
%     ylabel('y','fontsize',14)
% end

% print(figure(1),'-depsc','G_CFL10.eps')
% print(figure(2),'-depsc','agmgG_CFL10.eps')
