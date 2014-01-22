% test_agmg

function [results ref_results] = test_agmg

Kvec = 2.^(2:4);
Nvec = 8:16;
%Nvec = 4:8;
%Kvec = 2.^(2:4);

%Dt = 1*ones(length(Kvec),length(Nvec));   
%ref = run_test(Dt,1.0,0,Nvec,Kvec)

results = containers.Map;
ref_results = containers.Map;
[NN KK] = meshgrid(Nvec,Kvec);
for CFL = [1 5 10]
    for Re = [1e2 5e2 1e3]
        Dt = CFL./(NN.^2.*KK);
        key = ['CFL=',num2str(CFL),', Re=',num2str(Re)];
        results(key) = run_test(Dt,Re,1,Nvec,Kvec);
        ref_results(key) = run_test(Dt,Re,0,Nvec,Kvec);
    end
end

function Iters = run_test(Dt,Re,c,Nvec,Kvec)
%CFL = 10;
%nu = .0001;
nu = 1/Re;
if nargin<5
    Nvec = 1:2;
    Kvec = 2.^(2:3);
end
Iters = zeros(length(Kvec),length(Nvec));

for i = 1:length(Kvec)    
    Its = zeros(length(Nvec),1);
    for j = 1:length(Nvec)
        K = Kvec(i);N = Nvec(j);                      
        [Mg, Kg, Cg, bcInds, f, u0] = assemble(N,K); % get global SEM matrices

        %dt = CFL/(K*N^2);
        dt = Dt(i,j);
        
        % full matrix
        A = (1/dt)*Mg + nu*Kg + c*Cg;
        
        % time-rhs
        Rhs = (1/dt)*Mg*u0 + f;
        
        % bcs
        A(bcInds,:) = zeros(size(A(bcInds,:)));
        A(:,bcInds) = zeros(size(A(:,bcInds)));
        A(bcInds,bcInds) = eye(length(bcInds));
        Rhs(bcInds) = 0;
        
        levels = agmg_setup(A);
        
        tol = 1e-7;
        maxiter = 100;
        
        % solve
        [x flag relres iter resvec] = agmg_solve(levels, Rhs, maxiter, tol);
        %Iters(i,j) = iter;        
        Its(j) = iter;
    end    
    Iters(i,:) = Its;
end

