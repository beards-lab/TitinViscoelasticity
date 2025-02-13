%% parameter conversion

%% original
Fss  = mod(10);		% Fss   Steady state level
n_ss = mod(12);		% n_ss  Steady state level exponent
kp   = mod(9) ;		% kp    proximal chain force constantkS   =
np   = mod(3) ;		% np    proximal chain force exponent
kd   = mod(22);		% kd    proximal chain force constant high Cabist
nd   = mod(5) ;		% nd    distal chain force exponent
a_U  = mod(6) ;		% a_U   chain unfolding rate constant
nU   = mod(4) ;		% nU    unfolding rate exponent
mu   = mod(14);		% mu    small enough not to affect the result
delU = mod(19)/Ng;  % delU 
kA   = mod(7) ;     % kA   
kD   = mod(8) ; 	% kD    PEVK detachment rate
%% convert a mod line
param = convertModsToParams(mod);
fprintf('%10.4g ', param);  % Prints each row with 4 significant figures

%% go on all mods
ms = modSet;

for i = 1:size(ms, 1)
    paramSet(i, :) = convertModsToParams(ms(i, :));
end

for i = 1:size(paramSet,1)
    fprintf('%10.4g ', paramSet(i,:));  % Prints each row with 4 significant figures
    fprintf('\n');  % New line after each row
end
    
%% convert modsels

% modSel = [7 8 9]
paramSel = convertModSels(modSel)

%% convert modNames
mn = 1:23;
mnn = convertModsToParams(mn);
modNames(mnn)'

%% print out

for i = 1:size(modSet,1)
    fprintf('%10.4g ', modSet(i,:));  % Prints each row with 4 significant figures
    fprintf('\n');  % New line after each row
end


function paramsSel = convertModSels(modSel)
    modSelOld = zeros(23, 1);
    modSelOld(modSel) = 1;
    paramsSel =convertModsToParams(modSelOld);
    paramsSel = find(paramsSel == 1);
end

function param = convertModsToParams(mod)
    %% conversion table
    param(1 ) = mod(10);	% Fss   Steady state level
    param(2 ) = mod(12);	% n_ss  Steady state level exponent
    param(3 ) = mod(9) ;	% kp    proximal chain force constantkS   =
    param(4 ) = mod(3) ;	% np    proximal chain force exponent
    param(5 ) = mod(22);	% kd    proximal chain force constant high Cabist
    param(6 ) = mod(5) ;	% nd    distal chain force exponent
    param(7 ) = mod(6) ;	% a_U   chain unfolding rate constant
    param(8 ) = mod(4) ;	% nU    unfolding rate exponent
    param(9 ) = mod(14);	% mu    small enough not to affect the result
    param(10) = mod(19);    % delU 
    param(11) = mod(7) ;    % kA   
    param(12) = mod(8) ; 	% kD    PEVK detachment rate    
end