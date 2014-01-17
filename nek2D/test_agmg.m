% test_agmg

function results = test_agmg

results = containers.Map;
ref_results = containers.Map;
for CFL = 5:1:6
    for Re = 1000:2000:3000
        key = ['CFL=',num2str(CFL),', Re=',num2str(Re)];
        results(key) = run_test(CFL,Re,1);
        ref_results(key) = run_test(CFL,Re,0);
    end
end

function Iters = run_test(CFL,Re,c)
%CFL = 10;
%nu = .0001;
nu = 1/Re;
if nargin<5
    Nvec = 1:2;
    Kvec = 2.^(2:3);
end
Iters = zeros(length(Kvec),length(Nvec));

for i = 1:length(Kvec)
    for j = 1:length(Nvec)
        K = Kvec(i);N = Nvec(j);                      
        [Mg, Kg, Cg, bcInds, f, u0] = assemble(N,K); % get global SEM matrices

        dt = CFL/(K*N^2);
        
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
        Iters(i,j) = iter;
    end    
end

