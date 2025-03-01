function [rmin,rmax] = numericalBoundsShockRank(restr,phi,q0,opt,optimOptions)
% Function uses a constrained optimisation routine to find the upper and
% lower bounds of the conditional identified set.
% Inputs are:
% - restr: structure representing restrictions
% - phi: structure containing reduced-form VAR parameters
% - opt: structure containing model information and options
% - optimOptions: structure containing optimisation option

SA = restr.SA;
vma = phi.vma;
H = opt.H;
ivar = opt.ivar;

% Set equality constraint for optimisation (Aeq*Q(:) = beq).
Aeq = [];
beq = [];

% Set inequality constraints for optimisation (Aineq*q_1 <= b).
Aineq = -SA;
bineq = zeros(size(SA,1),1);

rmin = zeros(H+1,length(opt.ivar));
rmax = rmin;
LB = [];
UB = [];

parfor hh = 1:H+1 % For each horizon (change parfor to for if no Parallel Computing Toolbox)
    
    [rmin(hh,:),rmax(hh,:)] = genBoundsShockRank(q0,ivar,vma(:,:,hh),...
        Aineq,bineq,Aeq,beq,LB,UB,optimOptions)

end

rmax = -rmax; % Max obtained by minimising negative of objective

end