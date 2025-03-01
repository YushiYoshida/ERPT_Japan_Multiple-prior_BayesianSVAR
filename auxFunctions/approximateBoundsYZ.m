function [rmin,rmax,rminERPT,rmaxERPT] = approximateBoundsYZ(restr,phi,opt)

%function [rmin,rmax,rminERPT,rmaxERPT,Q0] = approximateBoundsYZ(restr,phi,opt)
%@YZ2025JIMF, The 2nd line is for retaining all Q0 by Yoshida

% Approximate bounds of the conditional identified set via simulation.
% Inputs:
% - restr: structure containing information about restrictions
% - phi: structure containing reduced-form VAR parameters
% - opt: structure containing model information and options

H = opt.H;
ivar = opt.ivar;
jshock = opt.jshock;
KTilde = opt.KTilde;
vma = phi.vma;
%cumIR = opt.cumIR;

% Draw from space of orthonormal matrices satisfying the restrictions 
% KTilde times.
%Q0 = drawQ0(restr,phi,opt);
Q0 = drawQ0eq(restr,phi,opt); %@YZ2025JIMF, with zero resrictions
% Keep columns of Q needed to compute impulse responses.
Q0 = Q0(:,jshock,:); 

% Compute impulse responses for each draw of Q.
etaDraw = zeros(H+1,length(ivar),length(jshock),KTilde);

for hh = 1:H+1 % For each horizon

    % Extract required rows of horizon-h VMA coefficient matrix.
    Cphi = vma(ivar,:,hh);
    % Multiply by relevant columns of Q0 (for each draw of Q0) to obtain
    % impulse response.
    etaDraw(hh,:,:,:) = pagemtimes(Cphi,Q0);

end
 
%if ~isempty(cumIR) %suppressed by Yoshida
    
    % Compute cumulative impulse responses for relevant variables.
%    etaDraw(:,cumIR,:,:) = cumsum(etaDraw(:,cumIR,:,:),1);

%end

% Compute minimum and maximum impulse response over draws of Q.
rmin = min(etaDraw,[],4);
rmax = max(etaDraw,[],4);

%@YZ2025JIMF, compute min, max of ERPT
A_yy1=ones(opt.H+1);
A_yy2=triu(A_yy1);
etaCum=pagemtimes(permute(etaDraw,[2 1 3 4]),A_yy2);
etaERPT(1,:,:,:)=etaCum(1,:,:,:)./etaCum(3,:,:,:);
rminERPT = min(etaERPT,[],4);
rmaxERPT = max(etaERPT,[],4);


end

