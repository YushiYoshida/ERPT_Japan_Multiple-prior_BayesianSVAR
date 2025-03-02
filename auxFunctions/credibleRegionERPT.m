function [credlb,credub] = credibleRegionERPT(rmin,rmax,opt)
% Compute robustified credible region with credibility level alpha.
% Inputs:
% - rmin: posterior draws of lower bound
% - rmax: posterior draws of upper bound
% - opt: structure containing options
%@YZ2025JIMF modified by Yoshida

aalpha = opt.aalpha; % Credibility level
gridLength = opt.gridLength; % Number of point on discrete grid

%@YZ2025JIMF, The codes below in variable placed are changed because ERPT does not have 'ni'
% dimension.
[K,H,nj] = size(rmin);

%@YZ2025JIMF, Storage. below modified by Yoshida
credlb = zeros(H,nj);
credub = credlb;

%@YZ2025JIMF, suppressed by Yoshida
%for ii = 1:ni % For each variable
    
    for jj = 1:nj % For each shock
    
        Cent = zeros(H,1);
        Rad = zeros(H,1);

        for hh = 1:H % For each horizon

            % Construct discrete grid.
            r = linspace(min(rmin(:,hh,jj)),max(rmax(:,hh,jj)),...
                gridLength);
            gridr = kron(r,ones(K,1));

            % d(r,phi) = max{|r-l(phi)|,|r-u(phi)|}.
            d = max(abs(gridr-kron(rmin(:,hh,jj),ones(1,gridLength))),...
                abs(gridr-kron(rmax(:,hh,jj),ones(1,gridLength))));
            zhat = quantile(d,aalpha,1);
            [Rad(hh),ind] = min(zhat,[],2);
            Cent(hh) = gridr(1,ind);

        end

        % Approximate robustified credible region is an interval centered at 
        % Cent with radius Rad.
        credlb(:,jj) = Cent - Rad;
        credub(:,jj) = Cent + Rad;
    
    end

%end % suppressed by Yoshida for ii
    
end
