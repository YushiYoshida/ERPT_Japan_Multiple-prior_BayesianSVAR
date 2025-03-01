function [Q,Ftilde,SS,S,empty] = drawQeq(restr,irs,phi,opt)
% Function attempts to draw Q from the space of orthonormal matrices 
% satisfying traditional and narrative restrictions.
% Inputs:
% - restr: structure containing information about restrictions
% - phi: structure containing reduced-form VAR parameters
% - opt: structure containing model information and options

%@YZ2025JIMF, Yoshida added the folloiwng as input
% - irs: structure containing IRs and possibly long-run (cumulative) IRs
% the original function is 'drawQ'

%@YZ2025JIMF, Yoshida added the following 4 lines.
f=restr.f;
Sigmatr = phi.Sigmatr;
eqRestr = restr.eqRestr;
lrcir = irs.lrcir;

shockSignRestr = restr.shockSignRestr;
hdSignRestr = restr.hdSignRestr;
signRestr = restr.signRestr;
Sigmatrinv = phi.Sigmatrinv;
vma = phi.vma;
U = Sigmatrinv*restr.U; % Multiply by Sigmatrinv to avoid doing this repeatedly
Sigmatrinvp = Sigmatrinv'; % Transpose to avoid doing this repeatedly

n = size(Sigmatrinv,1); % Number of variables in VAR
mSS = size(shockSignRestr,1); % No. of shock-sign restrictions
mHD = size(hdSignRestr,1); % No. of restrictions on historical decomp.
mS = size(signRestr,1); % No. of traditional sign restrictions

%@YZ2025JIMF, The following is added by Yoshida
mZ=size(eqRestr,1);

%% Construct matrix representing equality restrictions.
% F(phi) = [F_1(phi); ...; F_N(phi)], where F_i(phi) represent equality 
% restrictions on the ith column of Q, so F_i(phi)*q_i = 0.
%@YZ2025JIMF, Yoshida copied & pasted this section from GK's (2018) drawQ.m code.
%@YZ2025JIMF, 'N' is changed to 'n'

if mZ == 0 % If no equality restrictions

    F = [];

else % If equality restrictions

    F = zeros(mZ,n);

    for ii = 1:mZ

        if eqRestr(ii,4) == 1 % Restriction on A0

            F(ii,:) = Sigmatrinv(:,eqRestr(ii,2))';

        elseif eqRestr(ii,4) == 2 % Restriction on A0^(-1)

            F(ii,:) = Sigmatr(eqRestr(ii,1),:);

        elseif eqRestr(ii,4) == 3 % Restriction on A_l
            
            B = reshape(B,[n*p+const,n]);
            Bl = B(const+(eqRestr(ii,3)-1)*n+1:const+(eqRestr(ii,3))*n,...
                eqRestr(ii,2));
            F(ii,:) = (Sigmatrinv*Bl)';
            
        elseif eqRestr(ii,4) == 4 % Restriction on LRCIR

            F(ii,:) = lrcir(eqRestr(ii,1),:);

        end

    end

end

% Order rows of F in terms of the column of Q restricted.

Ftilde = [];

if isempty(F)
    
    Ftilde = F;
    
else

    for jj = 1:n % For each column of Q

       % Find rows of F restricting jth column of Q.
       Fj = F(eqRestr(:,1)==jj & (eqRestr(:,4)==1 | eqRestr(:,4)==3)...
           | eqRestr(:,2)==jj & (eqRestr(:,4)==2 | eqRestr(:,4)==4),:);
       Ftilde = [Ftilde; Fj];

    end
    
end


%% Construct matrix representing shock-sign restrictions.
% Restrictions are represented as SS(phi,U)*vec(Q) >= 0.
SS = zeros(mSS,n^2);

for ii = 1:mSS % For each restriction
    
    SS(ii,(shockSignRestr(ii,1)-1)*n+1:shockSignRestr(ii,1)*n) = ...
        U(:,shockSignRestr(ii,3))'*shockSignRestr(ii,2);
    
end

%% Construct matrix representing traditional sign restrictions.
% Restrictions are represented as S(phi)*vec(Q) >= 0.

S = zeros(mS,n^2);

for ii = 1:mS % For each restriction
    
    if signRestr(ii,5) == 1 % Sign restriction on impulse response
    
        S(ii,(signRestr(ii,2)-1)*n+1:signRestr(ii,2)*n) = ...
            vma(signRestr(ii,1),:,signRestr(ii,3)+1)*signRestr(ii,4);
    
    elseif signRestr(ii,5) == 2 % Sign restriction on A0
        
        S(ii,(signRestr(ii,2)-1)*n+1:signRestr(ii,2)*n) = ...
            Sigmatrinv(:,signRestr(ii,2))'*signRestr(ii,4);
        
    end
    
end

%% Attempt to draw Q satisfying restrictions.

iter = 0;
empty = 1;
A_yy1=ones(opt.H+1); %@YZ2025JIMF, added by Yoshida
A_yy2=triu(A_yy1);  %@YZ2025JIMF

while iter <= opt.maxDraw && empty == 1
    
    % YZ2025JIMF
    %The following is copied and pasted by Yoshida for zero sign
    %restrictions. 'N' is changed to 'n'.
    %
    % Draw Q from space of orthonormal matrices satisfying equality 
    % restrictions.

    % Generate vector of standard normal random variables.
    z = randn(n,1); 

    Qtilde = zeros(n);
    
    if isempty(F)

        Qtilde(:,1) = z;
        
    else
    
        % Find rows of F restricting first column of Q (F1).
        F1 = Ftilde(1:f(1),:);

        % Compute residual from linear projection of z on F1'.
        [q,r] = qr(F1',0);
        Qtilde(:,1) = z-F1'*(r \ (q'*z));

    end

    for jj = 2:n % For each column of Q

        if isempty(F)
            
            regMat = Qtilde(:,1:jj-1);
            
        else
        
           % Find rows of F restricting jth column of Q (Fj).
           Fj = Ftilde((sum(f(1:jj-1))+1):sum(f(1:jj)),:);
           regMat = [Fj' Qtilde(:,1:jj-1)];
        
        end

       % Generate vector of independent standard normal random variables.
       z = randn(n,1); 

       % Compute residual from linear projection of z on regMat.
       [q,r] = qr(regMat,0);
       Qtilde(:,jj) = z - regMat*(r \ (q'*z));

    end

%@YZ2025JIMF, suppressed by Yoshida
    % Draw Q from space of orthonormal matrices.
%    z = randn(n);
%    [Qtilde,~] = qr(z);
    
    % Normalise diagonal elements of A0 to be positive. Note that Matlab 
    % is implicitly expanding arrays to be compatible with elementwise 
    % array operations.
    Q = ((sign(diag(Sigmatrinvp*Qtilde))').*Qtilde)./vecnorm(Qtilde);
    
    % Check whether proposed draw satisfies shock-sign restrictions.
    shockSignRestrSat = all(SS*Q(:) >= 0);

    if shockSignRestrSat == 0  
        % If shock-sign restrictions not satisfied, increment iteration
        % counter and return to beginning of loop.
        iter = iter + 1;
        continue;          
    end
     
    % Check whether proposed draw satisfies traditional sign restrictions.
    signRestrSat = all(S*Q(:) >= 0);
        
    if signRestrSat == 0  
        % If sign restrictions not satisfied, increment iteration
        % counter and return to beginning of loop.
        iter = iter + 1;
        continue;          
    end
    
    % Check whether proposed draw satisfies restrictions on historical
    % decomposition.
    hdCheck = zeros(mHD,1); 

    for ii = 1:mHD % For each restriction  
    
        % Pre-compute q_j*q_j' for j=1,...,n.
        QQ = pagemtimes(reshape(Q,[n,1,n]),reshape(Q,[1,n,n]));

        % Extract row of hdSignRestr corresponding to ith restriction.
        hdRestr = hdSignRestr(ii,:);

        % Compute contribution of each shock.
        HH = zeros(1,1,n);    

        for hh = 1:hdSignRestr(ii,4)+1 % Sum contribution over horizons

            HH = HH + pagemtimes(pagemtimes(vma(hdRestr(1),:,hh),QQ),...
                U(:,hdRestr(3)+hh-1));

        end

        if hdRestr(5) == 1 && hdRestr(6) == 1

            % Type A - most important contributor
            hdCheck(ii) = abs(HH(hdRestr(2))) == max(abs(HH));

        elseif hdRestr(5) == 1 && hdRestr(6) == -1

            % Type A - least important contributor
            hdCheck(ii) = abs(HH(hdRestr(2))) == min(abs(HH));

        elseif hdRestr(5) == 2 && hdRestr(6) == 1

            % Type B - overwhelming contributor
            hdCheck(ii) = abs(HH(hdRestr(2))) - ...
                sum(abs(HH((1:n) ~= hdRestr(2)))) >= 0;                

        elseif hdRestr(5) == 2 && hdRestr(6) == -1

            % Type B - negligible contributor
            hdCheck(ii) = abs(HH(hdRestr(2))) - ...
                sum(abs(HH((1:n) ~= hdRestr(2)))) <= 0;

        end
        
        if hdCheck(ii) == 0              
            % If ith restriction on historical decomposition not 
            % satisfied, exit loop (saves computing contributions 
            % relating to subsequent restrictions).
            continue;
        end

    end

    hdSignRestrSat = all(hdCheck);
    
    % If all restrictions satisfied, terminate while loop. Otherwise,
    % increment counter and attempt to draw again.
    empty = 1 - shockSignRestrSat*hdSignRestrSat*signRestrSat;
    iter = iter + 1;
    
    % @YZ2025JIMF, added by Yoshida for ERPT range restrictions
    if empty==0
        for hh=1:opt.H+1
            temp_irf(hh,:,:)=phi.vma(opt.ivar,:,hh)*Q(:,opt.jshock);
        end
        temp_cirf=pagemtimes(permute(temp_irf,[2 1 3]),A_yy2);
        temp_erpt=temp_cirf(1,:,:)./temp_cirf(3,:,:);
        if temp_erpt>=-1&temp_erpt<=1
            empty=0;
        else
            empty=1;
        end
    end

end

end