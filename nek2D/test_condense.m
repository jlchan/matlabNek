addpath('./agmg')

% constants
one = @(x) ones(size(x));
zero = @(x) zeros(size(x));

N = 8; % order
K = 8;
% CFL = 1; dt = CFL/(K*N^2);
dt = 1.0;
eps = .01;
c = 1;

beta = {one,one,zero,zero};
%beta = {one,one,one, one};
%f = {one, one};

% manufactured solns
ue = {@(x) sin(pi*x), @(y) sin(pi*y)};
f = @(x,y) zeros(size(x));
%ue = {@(x) x.*(1-x), @(y) y.*(1-y)};
%f = @(x,y) -2*x.*(x-1) - 2*y.*(y-1);

% clean this part up - consolidate
Kx = K;Ky = K;Nq = N+1;Nqkx = Kx*(Nq-1) + 1; Nqky = Ky*(Nq-1) + 1; % num global dofs along line
[Mg, Kg, Cg, bcInds, fvec, u0, galnums] = assemble(N,K,beta,f); % get global SEM matrices
[X Y] = get_physical_points(N,Kx,Ky);

% initial condtion
x0 =0.5; y0=0.5; delta = 0.10; R=(X-x0).^2+(Y-y0).^2;
U0 = exp(-((R./(delta^2)).^1)).*X.*(1-X).*Y.*(1-Y); % pulse * bubble
u0 = reshape(U0,Nqkx*Nqky,1); 
fg_t =  fvec + (1/dt)*Mg*u0;
A = (1/dt)*Mg + eps*Kg + c*Cg;

lids = [];
gids = [];
for kx = 1:Kx
    for ky = 1:Ky
        for i = 1:Nq
            for j = 1:Nq
                elem_offset = ((ky-1)*Kx + (kx-1))*Nq*Nq;
                local_offset = ((i-1) + (j-1)*Nq) + 1;
                if ((i>1 && i<Nq) && (j>1 && j<Nq)) % if edge dof
                    lids = [lids galnums(elem_offset + local_offset)];
                else
                    gids = [gids galnums(elem_offset + local_offset)];                
                end
            end
        end
    end
end
gids = unique(gids);
As = A(lids,lids);Bs = A(lids,gids);Cs = A(gids,lids);Ds = A(gids,gids);
S = Ds - Cs*(As\Bs); % form schur complement
Rhs_S = fg_t(gids) - Cs*(As\fg_t(lids)); % form new load
Rhs = fg_t;

% agmg for original system
levels = agmg_setup(A);
tol = 1e-7; maxiter = 100;
% x = A\Rhs;
[x flag relres iter resvec] = agmg_solve(levels, Rhs, maxiter, tol);

% agmg for new system
levels_S = agmg_setup(S);
tol = 1e-7;maxiter = 100;
[x_edge flag relres iter_S resvec_S] = agmg_solve(levels_S, Rhs_S, maxiter, tol);
% x_edge = S\Rhs_S;

levels_As = agmg_setup(As);
[x_bub f rr it rv_b] = agmg_solve(levels_As, fg_t(lids)-Bs*x_edge, maxiter, tol);
%x_bub = As\(fg_t(lids)-Bs*x_edge);
xsc = zeros(size(x));
xsc([lids gids]) = [x_bub; x_edge];
disp('difference in solns')
diff = norm(xsc-x);

figure
semilogy(resvec);hold on
semilogy(resvec_S,'r')
title(['Difference in solution = ', num2str(diff)])

break

% plot change in spectra over eps
A1 = (1/dt)*Mg + .1*Kg + c*Cg;
A2 = (1/dt)*Mg + .01*Kg + c*Cg;
A3 = (1/dt)*Mg + .001*Kg + c*Cg;
lam1 = eig(full(A1));
lam2 = eig(full(A2));
lam3 = eig(full(A3));
figure;plot(lam1,'o','markersize',8)
hold on;
plot(lam2,'rs','markersize',8)
plot(lam3,'kd','markersize',8)
ax = axis;
As = A1(lids,lids);Bs = A1(lids,gids);Cs = A1(gids,lids);Ds = A1(gids,gids);
S1 = Ds - Cs*(As\Bs); % form schur complement
As = A2(lids,lids);Bs = A2(lids,gids);Cs = A2(gids,lids);Ds = A2(gids,gids);
S2 = Ds - Cs*(As\Bs); % form schur complement
As = A3(lids,lids);Bs = A3(lids,gids);Cs = A3(gids,lids);Ds = A3(gids,gids);
S3 = Ds - Cs*(As\Bs); % form schur complement
lam1 = eig(full(S1));
lam2 = eig(full(S2));
lam3 = eig(full(S3));
figure;plot(lam1,'o','markersize',8)
hold on;
plot(lam2,'rs','markersize',8)
plot(lam3,'kd','markersize',8)
axis(ax)

% % examining green's functions
% middle = gids(round(length(gids)/2));
% middle = lids(round(length(gids)/2));
% %middle = round(Nqkx*Nqky/2);
% e = zeros(Nqkx*Nqky,1);
% e(middle) = .01;
% G = A\e;
% G_agmg = precond_agmg(levels, e);
% 
% e_s = e(gids)-Cs*(As\e(lids));
% G_edge = precond_agmg(levels_S,e_s);
% G_bub = As\(fg_t(lids)-Bs*G_edge);
% Gsc = zeros(size(x));
% Gsc([lids gids]) = [G_bub; G_edge];
% 
% figure
% surf(X,Y,reshape(G,Nqkx,Nqky));
% shading interp;view(2);cax1 = caxis;colorbar
% 
% figure
% surf(X,Y,reshape(G_agmg,Nqkx,Nqky));
% shading interp;view(2);caxis(cax1);colorbar
% 
% figure
% surf(X,Y,reshape(Gsc,Nqkx,Nqky));
% shading interp;view(2);caxis(cax1);colorbar
