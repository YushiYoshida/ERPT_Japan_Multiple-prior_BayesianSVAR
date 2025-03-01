function [q1,SA,empty] = drawQShockRank(restr,phi,lpoptimOptions)
% Function checks whether there exists a value of q_1 satisfying the
% shock-rank restriction and traditional sign restrictions using the
% approach in Amir-Ahmadi and Drautzburg (2021).
% Inputs:
% - restr: structure containing information about restrictions
% - phi: structure containing reduced-form VAR parameters
% - opt: structure containing model information and options

shockRankRestr = restr.shockRankRestr;
U = restr.U;
signRestr = restr.signRestr;
Sigmatrinv = phi.Sigmatrinv;
vma = phi.vma;

n = size(Sigmatrinv,1); % Number of variables in VAR
T = size(restr.U,2); % Length of time series
mS = size(signRestr,1); % No. of traditional sign restrictions

%% Construct matrix representing shock-rank restriction.
% Each row of SR represents a restriction on the first structural shock
% such that SR(phi)*q_1 >= 0.

if ~isempty(shockRankRestr)
    
    inds = shockRankRestr(1); % Index of restricted period    
    otherInds = setdiff(1:T,inds); % Indices of other periods   
    SR = (Sigmatrinv*(U(:,inds)-U(:,otherInds)))'; 
               
    if shockRankRestr(1,2) == 1 % If absolute shock-rank restriction
        
        SR = [SR; (Sigmatrinv*(U(:,inds(1))+U(:,otherInds)))'];

    end
    
    SR = [SR; (Sigmatrinv*U(:,inds(1)))']; % Add shock-sign restriction
    SR = shockRankRestr(3)*SR; % Adjust for sign of shock    
           
else
    
    SR = [];
    
end

%% Construct matrix representing traditional sign restrictions.
% Restrictions are represented in S(phi), where S(phi)*q_1 >= 0.

S = zeros(mS,n);

for ii = 1:mS
    
    if signRestr(ii,5) == 1 % Sign restriction on impulse response
    
        S(ii,:) = vma(signRestr(ii,1),:,signRestr(ii,3)+1)*signRestr(ii,4);
    
    elseif signRestr(ii,5) == 2 % Sign restriction on A0
        
        S(ii,:) = Sigmatrinv(:,signRestr(ii,2))'*signRestr(ii,4);
        
    end
    
end

S = [S; Sigmatrinv(:,1)']; % Add sign normalisation

%% Check if set of q_1 satisfying restrictions is nonempty.
% This uses the approach in Amir-Ahmadi and Drautzburg (2021). We solve for
% the Chebyshev center of the polyhedron defined by the restrictions and
% the unit n-cube. If this ball has positive radius, the set is nonempty.

SA = [SR; S]; % Collect restrictions in single matrix
mSA = size(SA,1); % Total number of restrictions

A = zeros(mSA+1,n+1);
A(1:end-1,2:end) = -SA;

for ii = 1:(mSA+1)

    A(ii,1) = norm(A(ii,2:end));

end

A(end,1) = -1; % Add constraint that radius is positive

% Add additional restrictions (including that ball lies in unit n-cube).
A = [A; [ones(n,1), eye(n)]; [ones(n,1), -eye(n)]];

% Problem is to maximise f*x s.t. A*x <= b, with x = (R,q_1')'.
f = [-1, zeros(1,n)]; % 
b = [zeros(mSA+1,1); ones(2*n,1)];
x = linprog(f,A,b,[],[],[],[],lpoptimOptions);

if x(1) > 0 % Check if radius is positive (i.e. if set is nonempty)
    
    empty = 0;
    % Rescale Chebychev center to have unit norm. This is a value of 
    % q_1 that satisfies the restrictions.
    q1 = x(2:end)./norm(x(2:end));
    
else
    
    empty = 1;
    q1 = [];
    
end

end