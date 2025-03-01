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
nj = length(opt.jshock); % Number of shocks of interest
m = size(XX,2); % Number of parameters in each equation for y_t
nExog = opt.const + size(exog,2); % Number of exogenous variables
T = length(YY); % Number of observations used in estimating VAR

% Adjust indices representing narrative restrictions to account for
% losing the first p observations. If the restrictions are imposed in the 
% first p periods of the original dataset, drop the restrictions.
if ~isempty(restr.shockSignRestr)
    restr.shockSignRestr(:,3) = restr.shockSignRestr(:,3) - opt.p;
    restr.shockSignRestr(restr.shockSignRestr(:,3) <= 0) = [];
end
if ~isempty(restr.hdSignRestr)
    restr.hdSignRestr(:,3) = restr.hdSignRestr(:,3) - opt.p;
    restr.hdSignRestr(restr.hdSignRestr(:,3) <= 0) = [];
end

%% Conduct posterior inference on impulse responses.
% j-th column of B gives MLE of coefficients in j-th reduced form VAR 
% equation for y_t. Last rows correspond to exogenous terms (if included).
phiHat.B = (XX'*XX)\XX'*YY;
phiHat.S = (YY - XX*phiHat.B)'*(YY-XX*phiHat.B);
phiHat.Sigma = (1/T)*phiHat.S; % MLE of innovation covariance matrix
phiHat.P = (XX'*XX)\eye(m);
phiHat.cholP = chol(phiHat.P,'lower');

%@YZ2025JIMF, copied by Yoshida from GK (2018)
opt.whichAlgo = 0;
[opt.calcLRCIR,oneColQ,restr.f,restr.s,convInfo] = initialComp_YZ(n,opt,restr);


% Storage arrays.
rMinPost = zeros(opt.phiDraws,opt.H+1,ni,nj);
rMaxPost = zeros(opt.phiDraws,opt.H+1,ni,nj);
rSinglePriorPost = zeros(opt.phiDraws,opt.H+1,ni,nj);
B = zeros(n*m,opt.phiDraws);
Sigma = zeros(n,n,opt.phiDraws);
rMinERPT=zeros(opt.phiDraws,opt.H+1,nj);
rMaxERPT=zeros(opt.phiDraws,opt.H+1,nj);
Q00=zeros(opt.phiDraws, ni,nj,opt.KTilde);

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
    [irs,phi.vma,nonStable] = genVMA_YZ(phi,opt,nExog);
    nonStableCount = nonStableCount + nonStable;
    
    % If VAR is stable, proceed. Otherwise, redraw phi.
    
    if nonStable == 0
        
    phiDraw = phiDraw + 1;
    
    % Compute reduced-form VAR innovations.
    restr.U = (YY - XX*reshape(phi.B,m,n))';
        
    % Use Algorithm 1 (Step 2) to determine whether Q(phi|N,S) is nonempty
    % and draw Q satisfying restrictions.
    % By Yoshida, drawQ is changed to drawQeq
    [Qeq,restr.F,restr.SS,restr.S,Qempty(phiDraw)] = drawQeq(restr,irs,phi,opt);
               
    if Qempty(phiDraw) == 0 % If Q(phi|N,S) is non-empty
            
            % Store reduced-form parameters.
            nonEmptyDraw = nonEmptyDraw + 1;
            Sigma(:,:,nonEmptyDraw) = phi.Sigma;
            B(:,nonEmptyDraw) = phi.B;

            % Compute draw from posterior of impulse responses given 
            % conditionally uniform prior for Q|phi.
            
            for hh = 1:opt.H+1 % For each horizon
    
                % Compute impulse response.
                rSinglePriorPost(nonEmptyDraw,hh,:,:) = ...
                    phi.vma(opt.ivar,:,hh)*Qeq(:,opt.jshock);

            end
            
           % Find upper and lower bounds of conditional identified set for
           % impulse responses.
%           [rMinPost(nonEmptyDraw,:,:,:),rMaxPost(nonEmptyDraw,:,:,:),rMinERPT(nonEmptyDraw,:,:),rMaxERPT(nonEmptyDraw,:,:),Q00(nonEmptyDraw,:,:,:)] = ...
%               approximateBoundsYZSZ(restr,phi,opt); % for retaining all Q0
           [rMinPost(nonEmptyDraw,:,:,:),rMaxPost(nonEmptyDraw,:,:,:),rMinERPT(nonEmptyDraw,:,:),rMaxERPT(nonEmptyDraw,:,:)] = ...
               approximateBoundsYZ(restr,phi,opt);
           
           if mod(opt.phiDraws-nonEmptyDraw,opt.dispIter) == 0
              
               fprintf('\n%d draws with non-empty identified set remaining...',...
               opt.phiDraws-nonEmptyDraw);
               
           end

    end
    
    end
    
end

% Compute cumulative impulse response for draws under single prior where
% necessary.
if ~isempty(opt.cumIR)
    rSinglePriorPost(:,:,opt.cumIR,:) = ...
        cumsum(rSinglePriorPost(:,:,opt.cumIR,:),2);
end

% Compute posterior mean under single prior.
rSinglePriorMean = permute(mean(rSinglePriorPost,1),[2 3 4 1]);

% Compute lower and upper bound of set of posterior means.
meanlb = permute(mean(rMinPost,1),[2 3 4 1]);
meanub = permute(mean(rMaxPost,1),[2 3 4 1]);

%@YZ2025JIMF, compute lower and upper bounds of set of ERPT posterior means
meanERPTlb =permute(mean(rMinERPT,1),[2 3 1]);
meanERPTub =permute(mean(rMaxERPT,1),[2 3 1]);

% Compute robustified credible regions.
[credlb,credub] = credibleRegion(rMinPost,rMaxPost,opt);

% Compute highest posterior density (HPD) interval under single prior.
[hpdlb,hpdub] = highestPosteriorDensity(rSinglePriorPost,opt);

% Compute posterior plausibility of the identifying restrictions.
postPlaus = (1 - sum(Qempty)/length(Qempty));
fprintf('\nPosterior plausibility of identifying restrictions: %0.4g.\n',...
    postPlaus)

postMeanBoundWidth = meanub - meanlb; % Width of posterior mean bounds.
hpdWidth = hpdub - hpdlb; % Width of highest posterior density regions
credWidth = credub - credlb; % Width of robustified credible region

priorInformativeness = 1 - hpdWidth./credWidth; % Informativeness of prior

% Proportion of original draws where VAR is not stable.
nonStablePc = nonStableCount/(nonStableCount+phiDraw)*100;
fprintf('\nVAR representation unstable in %2.2f per cent of draws. \n',...
    nonStablePc);

runTime = toc