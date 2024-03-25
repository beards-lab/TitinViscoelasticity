function g  = dXdT(t,x,Nx,Ng,ds,kA,kD,kd,Fp,RU,RF,mu,L_0,nd,kDf, V)
s  = (0:1:Nx-1)'.*ds;

pu = reshape( x(1:(Ng+1)*Nx), [Nx,Ng+1]);
if size(x, 1) > Nx*(Ng+1) + 1 % we have pa
    pa = reshape( x((Ng+1)*Nx+1:2*(Ng+1)*Nx), [Nx,Ng+1]);
    g = zeros(2*(Ng+1)*Nx + 1,1);
else
    pa = [];
    g = zeros((Ng+1)*Nx + 1,1);
end
L = max(0,x(end));

% Calculate the un-attached chain velocities for every pu(s,n) entry
deltaF = kd*max(0,(L-s)/L_0).^nd - Fp;
Vp = deltaF/mu;

ij = (1:Nx)' + Nx*(0:Ng); % matrix of Ng X Nx indices over all elements

% Attach/dettach rectifier connector - only for Ca
if ~isempty(pa)
    detach = kD*x(ij + (Ng+1)*Nx).*(1 + kDf*max(0, deltaF+Fp));
    g(ij)             = g(ij)             - kA*x(ij) + detach;
    g(ij + (Ng+1)*Nx) = g(ij + (Ng+1)*Nx) + kA*x(ij) - detach;

% g(ij)             = g(ij)             - kA*x(ij) + kD*x(ij + (Ng+1)*Nx).*max(0,deltaF);
% g(ij + (Ng+1)*Nx) = g(ij + (Ng+1)*Nx) + kA*x(ij) - kD*x(ij + (Ng+1)*Nx).*max(0,deltaF);

end

% UPWIND differencing for sliding (+ direction) for pu
g(ij(1,:)) = g(ij(1,:)) - (1/ds)*(pu(1,:).*max(0,Vp(1,:)))';

% positive velocities - extending
indxs = ij(2:Nx,:);
g(indxs) = g(indxs) ...
    - (1/ds)*pu(indxs).*max(0,Vp(2:Nx, :)) + (1/ds)*pu(indxs-1).*max(0,Vp(1:Nx-1,:));

% like:
% k(1) = k(1) - pu(1)*V(1)/ds;
% k(2:Nx) = k(2:Nx) - (pu(2:Nx)*V(2:Nx) - pu(2:Nx)*Vp(1:Nx-1))/ds

%%
% negative velocities - shrinking
indxs = ij(1:Nx-1,:);
g(indxs) = g(indxs) ...
    - (1/ds)*pu(indxs+1).*min(0,Vp(2:Nx,:)) + (1/ds)*pu(indxs).*min(0,Vp(1:Nx-1,:));

g(ij(end,:)) = g(ij(end,:)) + (1/ds)*(pu(end,:).*min(0,Vp(end,:)))';

% unfolding rate for pu states
UR = RU.*pu(ij(:,1:Ng)); % rate of probabilty transitions from n to n+1 states

% refolding rate - speed up 
if RF > 0
    FR = RF.*pu(ij(:,2:Ng+1)); % rate of probabilty transitions from n+1 to n states
else
    FR = 0;
end

if t > 9 
    a = 1;
end
g(ij(:,2:(Ng+1))) = g(ij(:,2:(Ng+1))) + UR - FR;
g(ij(:,1:Ng))     = g(ij(:,1:Ng))     - UR + FR;

% unfolding for pa states
if false || ~isempty(pa)
    NxNg = Nx*(Ng+1);
    UR = RU.*pa(ij(:,1:Ng)); % rate of probabilty transitions from n to n+1 states
    g(ij(:,2:(Ng+1))+NxNg) = g(ij(:,2:(Ng+1))+NxNg) + UR - FR;
    g(ij(:,1:Ng)+NxNg)     = g(ij(:,1:Ng)+NxNg)     - UR + FR;
end

if isa(V, 'function_handle')
    g(end) = V(t);
else
    g(end) = V;
end
