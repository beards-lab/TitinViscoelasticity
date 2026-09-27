function g  = dXdTvar(t,x,Nx,Ng,ds,kA,kD,kd,Fp,RU,RF,mu,L_0,nd,kDf, V, cComp, muRec, kDslack)
% dXdTvar  Copy of dXdT.m for structural experiments (FitFirstStretch.m):
% mu may be an Nx or Nx x (Ng+1) array (strain-dependent viscosity); cComp in
% [0,1] scales the distal element's compressive branch (1 = dXdT, 0 =
% tension-only, slack under compression); muRec (optional) is the drag for
% recoiling chains (Vp < 0); kDslack (optional) is an extra detachment rate for
% attached chains whose distal segment is slack (L < s). Otherwise identical.
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
if nargin >= 17 && cComp ~= 1
    deltaF = kd*(max(0, (L-s)/L_0).^nd - cComp*max(0, (s-L)/L_0).^nd) - Fp;
else
    deltaF = kd*sign(L-s).*abs((L-s)/L_0).^nd - Fp;
end
Vp = deltaF./mu;
if nargin >= 18 && ~isnan(muRec)
    rec = deltaF < 0;
    Vp(rec) = deltaF(rec)/muRec;
end

ij = (1:Nx)' + Nx*(0:Ng); % matrix of Ng X Nx indices over all elements
% ij_att = ij + (Ng+1)*Nx; % matrix of attached elements

% Attach/dettach rectifier connector - only for Ca
if ~isempty(pa)
    % detach = kD*x(ij + (Ng+1)*Nx).*(1 + kDf*(deltaF + Fp));
    if kDf > 0 % force-dependent (slip-bond) detachment, distal force deltaF + Fp
        detach = kD*x(ij + (Ng+1)*Nx).*(1 + kDf*max(0, deltaF + Fp));
    else
        detach = kD*x(ij + (Ng+1)*Nx);
    end
    if nargin >= 19 && kDslack > 0 % bonds fail when the distal segment is unloaded
        detach = detach + kDslack*x(ij + (Ng+1)*Nx).*(L < s);
    end
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
if ~isscalar(RF)
    % smooth (force-suppressed) refolding: RF is an Nx x Ng rate matrix
    FR = RF.*pu(ij(:,2:Ng+1));
elseif RF > 0
    FR = RF.*pu(ij(:,2:Ng+1)); % rate of probabilty transitions from n+1 to n states
    FR(Fp(:, 2:11) > 0) = 0;
else
    FR = 0;
end

if t > 150
    a = 1;
end
g(ij(:,2:(Ng+1))) = g(ij(:,2:(Ng+1))) + UR - FR;
g(ij(:,1:Ng))     = g(ij(:,1:Ng))     - UR + FR;

% unfolding for pa states

% TODO attached states do not ahve to unfold - ok simplification?
if ~isempty(pa)
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
