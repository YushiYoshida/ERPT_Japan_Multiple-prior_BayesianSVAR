function Q0 = drawQ0eq(restr,phi,opt)
% Function draws Q KTilde times from the space of orthonormal matrices 
% satisfying the traditional and narrative restrictions. 
% Inputs:
% - restr: structure containing information about restrictions
% - phi: structure containing reduced-form VAR parameters
% - opt: structure containing model information and options

%@YZ2025JIMF, modifed by Yoshida to add zero-restrictions
% original function is drawQ0

hdSignRestr = restr.hdSignRestr;
SS = restr.SS;
S = restr.S;
Ftilde = restr.F; %@YZ2025JIMF, added by Yoshida for zero restrictions
f = restr.f; %@YZ2025JIMF, by Yoshida
Sigmatrinv = phi.Sigmatrinv;
vma = phi.vma;
U = Sigmatrinv*restr.U; % Multiply by Sigmatrinv to avoid doing this repeatedly
Sigmatrinvp = Sigmatrinv'; % Transpose to avoid doing this repeatedly

n = size(Sigmatrinv,1); % Number of variables in VAR
mHD = size(hdSignRestr,1); % No. of restrictions on historical decomp

%% Draw KTilde values of Q satisfying sign restrictions.
A_yy1=ones(opt.H+1);%@YZ2025JIMF
A_yy2=triu(A_yy1); %@YZ2025JIMF, added by Yoshida to construct cumulated IRs
ni=length(opt.ivar);
nj=length(opt.jshock);

Q0 = zeros(n,n,opt.KTilde); % Storage

parfor kk = 1:opt.KTilde % Change from parfor to for if Parallel Computing Toolbox not installed

    flag = 0;

    while flag == 0

        % Draw Q from space of orthonormal matrices.
        %z = randn(n);       %@YZ2025JIMF, suppresed by Yoshida
        %[Qtilde,~] = qr(z); %@YZ2025JIMF, suppresed by Yoshida
    % the following block of codes are copied from drawQeq.m which is also copied from GK files.
    z = randn(n,1); 

    Qtilde = zeros(n);
    
    if isempty(Ftilde) %@YZ2025JIMF, changed from 'F' to 'Ftilde' by Yoshida

        Qtilde(:,1) = z;
        
    else
    
        % Find rows of F restricting first column of Q (F1).
        F1 = Ftilde(1:f(1),:);

        % Compute residual from linear projection of z on F1'.
        [q,r] = qr(F1',0);
        Qtilde(:,1) = z-F1'*(r \ (q'*z));

    end

    for jj = 2:n % For each column of Q

        if isempty(Ftilde) %@YZ2025JIMF, changed from 'F' to 'Ftilde' by Yoshida
            
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





        % Normalise diagonal elements of A0 to be positive. Note that Matlab 
        % is implicitly expanding arrays to be compatible with elementwise 
        % array operations.
        Q = ((sign(diag(Sigmatrinvp*Qtilde))').*Qtilde)./vecnorm(Qtilde);

        % Check whether proposed draw satisfies shock-sign restrictions.
        shockSignRestrSat = all(SS*Q(:) >= 0);

        if shockSignRestrSat == 0  
            % If shock-sign restrictions not satisfied, return to 
            % beginning of loop.
            continue;          
        end

        % Check whether proposed draw satisfies traditional sign restrictions.
        signRestrSat = all(S*Q(:) >= 0);

        if signRestrSat == 0  
            % If sign restrictions not satisfied, return to beginning of loop.
            continue;          
        end
        
        % Check whether proposed draw satisfies restrictions on historical
        % decomposition.
        hdCheck = zeros(mHD,1); 

        % Pre-compute q_j*q_j' for j=1,...,n.
        QQ = pagemtimes(reshape(Q,[n,1,n]),reshape(Q,[1,n,n]));

        for ii = 1:mHD % For each restriction

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
        
        % If restrictions satisfied, save value of Q and terminate 
        % while loop.
        flag = shockSignRestrSat*hdSignRestrSat*signRestrSat;
        
        if flag == 1
            temp_irf=zeros(opt.H+1,ni,nj);
            for hh=1:opt.H+1
                temp_irf(hh,:,:)=phi.vma(opt.ivar,:,hh)*Q(:,opt.jshock);
            end
            temp_cirf=pagemtimes(permute(temp_irf,[2 1 3]),A_yy2);
            temp_erpt=temp_cirf(1,:,:)./temp_cirf(3,:,:);
            if temp_erpt>=-1&temp_erpt<=1
                flag=1;
                Q0(:,:,kk) = Q;
            else
                flag=0;
            end                
         end
        
    end
    
end

end