function quadNodalExample

addpath('BookCodes1D')

% set up 1D operations
N = 12; Np = N+1; Np2 = Np^2;

[r] = JacobiGL(0,0,N);
[r s] = meshgrid(r);r = r(:);s = s(:);

% get derivative/mass matrices
V = Vander2D(r,s,N); invV = inv(V); 
[Dx Dy] = Grad2D(r,s,invV,N);
M = invV'*invV;

% define stiffness matrices
% K = Dx'*M*Dx + Dy'*M*Dy; % Poisson
S = (Dx+.5*Dy)'*M*(Dx+.5*Dy);
S = S + (.1*Dx+Dy)'*M*(.1*Dx+Dy);

f = zeros(Np2,1);
f = ones(Np2,1);

A = S;
b = f;

% bcs
tol = 1e-14;
onB = @(x,y) (abs(1-x.^2) < tol) | (abs(1-y.^2) < tol);
u0 = (1-r.*s).*onB(r,s);
bMap = find(onB(r,s));
b = b-A*u0;
for i = 1:length(bMap)
    A(bMap(i),:) = 0; A(:,bMap(i)) = 0;
    A(bMap(i),bMap(i)) = 1; b(bMap(i)) = u0(bMap(i));
end

% solve
u = A\b;

% % get smallest generalized eigval for coercivity
% gam = min(eig(A,M));

% uniform plotting grid
[ru] = linspace(-1,1,120); [ru su] = meshgrid(ru);ru = ru(:);su = su(:);
Vu = Vander2D(ru,su,N); Vp = Vu*invV; 
color_line3(ru,su,Vp*u,Vp*u,'.');
hold on
color_line3(r(bMap),s(bMap),u(bMap),u(bMap),'o');
% title(['Coercivity const = ', num2str(gam)])
% figure;
% color_line3(ru,su,Vp*f,Vp*f,'.')

function V = Vander2D(x,y,N)

Np = N+1;
V = zeros(length(x),Np^2);
k = 1;

for i = 0:N % over x
    Px = JacobiP(x,0,0,i);
    for j = 0:N % over y        
        V(:,k) = Px.*JacobiP(y,0,0,j);
        k = k+1;
    end
end

function [Dx Dy] = Grad2D(x,y,invV,N)

Np = N+1;
Dx = zeros(length(x),Np^2);
Dy = zeros(length(x),Np^2);
k = 1;
for i = 0:N % over x
    Px = JacobiP(x,0,0,i);
    DPx = GradJacobiP(x,0,0,i);
    for j = 0:N % over y        
        Dx(:,k) = DPx.*JacobiP(y,0,0,j);
        Dy(:,k) = Px.*GradJacobiP(y,0,0,j);
        k = k+1;
    end
end
Dx = Dx*invV;
Dy = Dy*invV;
