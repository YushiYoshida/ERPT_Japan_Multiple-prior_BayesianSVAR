tic

%% Construct data for estimating VAR for y_t.
YY = data(opt.p+1:end,:); % y_t 
XX = lagmatrix(data,1:opt.p); % Matrix of regressors in VAR for y_t
XX = XX(opt.p+1:end,:); % Drop initial missing observations
if opt.const == 1 % Add constant to matrix of regressors
    XX = [XX ones(size(XX,1),1)];
end
% Add exogenous variables to matrix of regressors
XX = [XX exog(opt.p+1:end,:)]; 

n = size(YY,2); % Number of variables
ni = length(opt.ivar); % Number of variables of interest
m = size(XX,2); % Number of parameters in each equation for y_t
nExog = opt.const + size(exog,2); % Number of exogenous variables
T = length(YY); % Number of observations used in estimating VAR

% Adjust indices representing narrative restrictions to account for
% losing the first p observations.
if ~isempty(restr.shockRankRestr)
    restr.shockRankRestr(:,1) = restr.shockRankRestr(:,1) - opt.p;
end

%% Conduct posterior inference on impulse responses.
% j-th column of B gives MLE of coefficients in j-th reduced form VAR 
% equation for y_t. Last rows correspond to exogenous terms (if included).
phiHat.B = (XX'*XX)\XX'*YY;
phiHat.S = (YY - XX*phiHat.B)'*(YY-XX*phiHat.B);
phiHat.Sigma = (1/T)*phiHat.S; % MLE of innovation covariance matrix
phiHat.P = (XX'*XX)\eye(m);
phiHat.cholP = chol(phiHat.P,'lower');

% Storage arrays.
rMinPost = zeros(opt.phiDraws,opt.H+1,ni);
rMaxPost = zeros(opt.phiDraws,opt.H+1,ni);
B = zeros(n*m,opt.phiDraws);
Sigma = zeros(n,n,opt.phiDraws);

phiDraw = 0; % Counter for no. of draws from posterior of phi
nonStableCount = 0; % Counter for no. of draws with nonstable VAR
nonEmptyDraw = 0; % Counter for no. of draws with non-empty IS

while nonEmptyDraw < opt.phiDraws

    % Posterior sampler with independent improper (Jeffreys) prior.
    phi.Sigma = iwishrnd(phiHat.S,T-m);
    phi.Sigmatr = chol(phi.Sigma,'lower');
    phi.Sigmatrinv = phi.Sigmatr\eye(n);
    phi.B = phiHat.B(:) + kron(phi.Sigmatr,phiHat.cholP)*randn(m*n,1);

    % Generate coefficients in orthogonal reduced-form VMA representation 
    % and check stability of VAR.
    [phi.vma,nonStable] = genVMA(phi,opt,nExog);
    nonStableCount = nonStableCount + nonStable;
    
    % If VAR is stable, proceed. Otherwise, redraw phi.
    
    if nonStable == 0
        
    phiDraw = phiDraw + 1;
    
    % Compute reduced-form VAR innovations.
    restr.U = (YY - XX*reshape(phi.B,m,n))';
        
    % Use approach in Amir-Ahmadi and Drautzburg (2021) to determine
    % whether Q(phi|N,S) is nonempty and save a value of q_1 satisfying
    % restrictions.
    [q0,restr.SA,Qempty(phiDraw)] = ...
        drawQShockRank(restr,phi,lpoptimOptions);
               
    if Qempty(phiDraw) == 0 % If Q(phi|N,S) is non-empty
            
            % Store reduced-form parameters.
            nonEmptyDraw = nonEmptyDraw + 1;
            Sigma(:,:,nonEmptyDraw) = phi.Sigma;
            B(:,nonEmptyDraw) = phi.B;
            
           % Find upper and lower bounds of conditional identified set for
           % impulse responses.
           [rMinPost(nonEmptyDraw,:,:,:),rMaxPost(nonEmptyDraw,:,:,:)] = ...
               numericalBoundsShockRank(restr,phi,q0,opt,optimOptions);              
                         
           if mod(opt.phiDraws-nonEmptyDraw,opt.dispIter) == 0
              
               fprintf('\n%d draws with non-empty identified set remaining...',...
               opt.phiDraws-nonEmptyDraw);
               
           end

    end
    
    end
    
end

% Compute lower and upper bound of set of posterior means.
meanlb = permute(mean(rMinPost,1),[2 3 1]);
meanub = permute(mean(rMaxPost,1),[2 3 1]);

% Compute robustified credible regions.
[credlb,credub] = credibleRegion(rMinPost,rMaxPost,opt);

% Compute posterior plausibility of the identifying restrictions.
postPlaus = (1 - sum(Qempty)/length(Qempty));
fprintf('\nPosterior plausibility of identifying restrictions: %0.4g.\n',...
    postPlaus)

% Proportion of original draws where VAR is not stable.
nonStablePc = nonStableCount/(nonStableCount+phiDraw)*100;
fprintf('\nVAR representation unstable in %2.2f per cent of draws. \n',...
    nonStablePc);

runTime = toc