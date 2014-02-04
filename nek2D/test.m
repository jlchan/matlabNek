clear
% constants
one = @(x) ones(size(x));
zero = @(x) zeros(size(x));

N = 8; % order
K = 16;
CFL = 10; 
dt = 1.0;%CFL/(K*N^2);
%dt = 1e-2;
%dt = 100;
eps = .001;
c = 1;

%beta = {one,one,zero,zero};
beta = {one,one,zero,zero};
%f = {one, one};

% manufactured solns
ue = {@(x) sin(pi*x), @(y) sin(pi*y)};
f = @(x,y) ones(size(x));
%ue = {@(x) x.*(1-x), @(y) y.*(1-y)};
%f = @(x,y) -2*x.*(x-1) - 2*y.*(y-1);

% clean this part up - consolidate
Kx = K;Ky = K;Nq = N+1;Nqkx = Kx*(Nq-1) + 1; Nqky = Ky*(Nq-1) + 1; % num global dofs along line
[Mg, Kg, Cg, bcInds, fvec, u0, galnums] = assemble(N,K,beta,f,true); % get global SEM matrices
[X Y] = get_physical_points(N,Kx,Ky);

% initial condtion
x0 =0.6; y0=0.3; delta = 0.10; R=(X-x0).^2+(Y-y0).^2;
U0 = exp(-((R./(delta^2)).^1)).*X.*(1-X).*Y.*(1-Y); % pulse * bubble
u0 = reshape(U0,Nqkx*Nqky,1); 

fg_t =  fvec + (1/dt)*Mg*u0;
A = (1/dt)*Mg + eps*Kg + c*Cg;

lids = []; % bubbles
gids = []; % edge dofs
vids = []; % purely vertex dofs
for kx = 1:Kx
    for ky = 1:Ky
        for i = 1:Nq
            for j = 1:Nq
                elem_offset = ((ky-1)*Kx + (kx-1))*Nq*Nq;
                local_offset = ((i-1) + (j-1)*Nq) + 1;
                if ((i>1 && i<Nq) && (j>1 && j<Nq)) % if not edge dof
                    lids = [lids galnums(elem_offset + local_offset)];
                else
                    % vertex ids
                    if ((i==1 && j==1) || (i==1 && j==Nq)...
                            || (i==Nq && j==1) || (i==Nq && j==Nq))
                        vids = [vids galnums(elem_offset+local_offset)];
                    end                    
                    gids = [gids galnums(elem_offset + local_offset)];                
                end
            end
        end
    end
end
gids = unique(gids);
vids = unique(vids); 
eids = setdiff(gids,vids);
if N>1
    As = A(lids,lids);Bs = A(lids,gids);Cs = A(gids,lids);Ds = A(gids,gids);
    S = Ds - Cs*(As\Bs); % form schur complement
    fg_t = fvec;
    Rhs_S = fg_t(gids) - Cs*(As\fg_t(lids)); % form new load    
    ugids = S\Rhs_S;    
else
    uglob = A\fvec;
    ug = uglob;
end
ug(lids) = nan;
ug(gids) = ugids;
ug(setdiff(1:length(ug),vids)) = nan;

scatter3(X(:),Y(:),ug,ones(size(ug))*24,ug,'filled')

if (N>1)
    si = zeros(size(S,1),1);
    i = 1;
    for g = gids
        si(g) = i;
        i=i+1;
    end
    p = [si(eids); si(vids)];
    figure
    spy(S(p,p));
end

break
% time-dependent
figure
pcolor(X,Y,U0)
ax = axis;
cax = caxis;
pause
Nsteps = 1;
[L U] = lu(A);
for i = 1:Nsteps
    % implicit
    ug = U\(L\fg_t);
    fg_t = fvec + (1/dt)*Mg*ug; % next timestep
    fg_t(bcInds) = 0;
    
    % explicit
    %ug = dt*(1./diag(Mg)).*fg_t;
    %fg_t = fg + (1/dt)*Mg*ug - Kg*ug;
    
    if (mod(i,Nsteps/10)==0)
        surf(X,Y,reshape(ug,Nqkx,Nqky))
%         caxis(cax)
        colorbar
        title(['Time = ',num2str(i*dt)])
%         axis(ax)        
        drawnow
    end
end
