function [rmin,rmax] = genBoundsShockRank(q0,ivar,vma,Aineq,bineq,Aeq,...
                                            beq,LB,UB,optimOptions)

rmin = zeros(1,length(ivar));
rmax = zeros(1,length(ivar));

for ii = 1:length(ivar) % For each variable of interest

    % Minimise impulse response subject to restrictions.
    [~,rmin(ii)] = fmincon(@(q) genIRFShockRank(q,vma(ivar(ii),:),1),...
        q0,Aineq,bineq,Aeq,beq,LB,UB,@(q) qcon(q),optimOptions);
    % Maximise impulse response subject to restrictions.
    [~,rmax(ii)] = fmincon(@(q) genIRFShockRank(q,vma(ivar(ii),:),-1),...
        q0,Aineq,bineq,Aeq,beq,LB,UB,@(q) qcon(q),optimOptions);

end