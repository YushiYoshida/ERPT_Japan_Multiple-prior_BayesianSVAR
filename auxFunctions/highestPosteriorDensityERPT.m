function [hpdlb,hpdub] = highestPosteriorDensityERPT(rPost,opt)
% Compute highest posterior density interval with credibility level alpha.
% Inputs:
% - rPost: posterior draws
% - opt: structure containing options
%@YZ2025JIMF, modified by Yoshida to omit (ni) for ERPT.

aalpha = opt.aalpha; % Credibility level
gridLength = opt.gridLength; % Number of point on discrete grid

[K,H,nj] = size(rPost);


hpdlb = zeros(H,nj);
hpdub = hpdlb;

%for ii = 1:ni % For each variable %YZ2025JIMF

    for jj = 1:nj % For each shock
    
        % Storage matrices
        Cent = zeros(H,1);
        Rad = zeros(H,1);

        for hh = 1:H % For each horizon

            % Construct grid over which to search for bounds.
            r = linspace(min(rPost(:,hh,jj)),max(rPost(:,hh,jj)),...
                gridLength);
            gridr = kron(r,ones(K,1));

            % d(r,phi) = |r-r(phi)|
            d = abs(gridr-kron(rPost(:,hh,jj),ones(1,gridLength)));
            zhat = quantile(d,aalpha,1);
            [Rad(hh),ind] = min(zhat,[],2);
            Cent(hh) = gridr(1,ind);

        end

        % Approximate highest posterior density region is an interval centered 
        % at Cent with radius Rad.
        hpdlb(:,jj) = Cent - Rad;
        hpdub(:,jj) = Cent + Rad;
    
    %end

end
