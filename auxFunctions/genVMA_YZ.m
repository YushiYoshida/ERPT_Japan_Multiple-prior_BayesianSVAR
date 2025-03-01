function [irs,vma,nonStable] = genVMA_YZ(phi,opt,nExog)
%@YZ2025JIMF, Modified by Yushi Yoshida

% Given reduced-form VAR parameters, generate coefficients in orthogonal
% reduced-form vector moving average (VMA) representation and check 
% stability of the VAR representation.
% Inputs:
% - phi: structure containing reduced-form VAR parameters
% - opt: structure containing model information and options
% - nExog: number of exogenous variables.

Sigmatr = phi.Sigmatr;
B = phi.B;
p = opt.p;
H = opt.H;
%cumIR=opt.cumIR
    
n = size(Sigmatr,1); % Number of variables in VAR
B_mat = reshape(B,size(B,1)/n,n); % Reshape coefficients into matrix
B_mat = B_mat(1:end-nExog,:); % Drop coefficients on exogenous variables

% Put VAR(p) in companion form (i.e. VAR(1) for (y_t,y_t-1,...,y_t-p+1)').
B_c = zeros(n*p);
B_c(1:n,:) = B_mat';   
B_c(n+1:end,1:end-n) = eye(n*p-n);

% Compute coefficients in vector moving average representation.
vma = zeros(n,n,H+1); % Storage array 
CC = eye(n*p);
vma(:,:,1) = Sigmatr;

for jj = 2:H+1 % For each horizon

    CC = CC*B_c;
    % Extract upper-left nxn block of coefficient matrix.
    vma(:,:,jj) = CC(1:n,1:n)*Sigmatr;
    
end

% The following two if-end sets are copied with sligth modifications from GK(2018) by Yoshida
if opt.calcLRCIR == 1 % Compute long-run IR

    lrcir = (eye(n)-sum(reshape(B_c(1:n,:),[n,n,p]),3))\eye(n);
    lrcir = lrcir*Sigmatr;

else
    
    lrcir = [];

end

%if ~isempty(cumIR) %@YZ2025JIMF, suppressed by Yoshida
    
    % Compute cumulative IRF for relevant variables
%    vma(cumIR,:,:) = cumsum(vma(cumIR,:,:),3);
    
%end

%@YZ2025JIMF, Copied by Yoshida
% Add computed IRs to structure
irs.vma = vma;
irs.lrcir = lrcir;


% Compute eigenvalues of companion matrix.
E = eig(B_c);
% Check if any eigenvalues are greater than one in modulus.
nonStable = any(abs(E) >= 1);

end