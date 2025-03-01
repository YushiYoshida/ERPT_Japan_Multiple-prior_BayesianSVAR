function [calcLRCIR,oneColQ,f,s,convInfo] = initialComp(n,opt,restr)
%@YZ2025JIMF, This function of GK (2018) is copied by Yoshida

% This function: computes flag to indicate whether long-run cumulative IR 
% (LRCIR) should be computed; checks that the order of variables satisfies 
% Definition 3; checks that the SVAR is, in fact, set identified using 
% condition 4.12; checks whether the restrictions apply to only one column 
% of Q (so that the approach of GMM18 is applicable); and does some 
% initial computations related to the convexity of the identified set.
% Inputs:
% - N: number of variables in VAR
% - opt: structure containing model information and options
% - restr: structure containing matrices representing restrictions            

eqRestr = restr.eqRestr;
signRestr = restr.signRestr;
mZ = size(eqRestr,1); % No. of equality restrictions
mS = size(signRestr,1); % No. of sign restrictions (excluding normalisation)
jshock = opt.jshock;
whichAlgo = opt.whichAlgo;
% Generate flag to indicate whether LRCIR should be computed..
calcLRCIR = 0;

if mZ > 0
    calcLRCIR = any(eqRestr(:,4) == 4); % Flag to compute LRCIR
end

%% Check that order of variables satisfies Definition 3.
% Also check that the SVAR is actually set identified.
f = zeros(n,1);

if mZ > 0 % If there are equality restrictions

    for ii = 1:n % For each column of Q

        % Find number of equality restrictions on ith column of Q.
        f(ii) = sum((eqRestr(:,1)==ii & ...
        (eqRestr(:,4)==1 | eqRestr(:,4)==3)) | ...
        (eqRestr(:,2)==ii & (eqRestr(:,4)==2 | eqRestr(:,4)==4))); 

    end

end
    
% Check that order of variables satisfies Definition 3.
if any(diff(flipud(f)) < 0)

    fprintf('\nOrder of variables does not satisfy f_1 >= f_2 >= ... >= f_N >= 0.\n');
    error('Error: Order of variables violates Definition 3!');

elseif find(f == f(jshock),1) ~= jshock

    fprintf('\nj*-th variable is not ordered first among ties.\n');
    error('Error: Order of variables violates Definition 3!');

else

    fprintf('\nOrder of variables satisfies Definition 3.\n');

end

% Check that SVAR is, in fact, set identified using condition 4.12:
% f_i <= N-i for i = 1,...,N with strict inequality for at least one i.
if ~(all(f <= n - (1:n)') && any(f < n - (1:n)'))

    fprintf('\nEquality restrictions do not satisfy f_i <= N-i for i = 1,...,N, with strict inequality for at least one i.\n');
    error('Error:SVAR is not set identified according to condition 4.12!.');

end

%% Check whether restrictions relate to only one column of Q.
s = zeros(n,1);

if mS > 0 % If there are sign restrictions
       
    for ii = 1:n % For each column of Q.
        
        % Count number of inequality restrictions on ith column of Q
        s(ii) = sum((signRestr(:,2)==ii & signRestr(:,5)==1) | ...
            (signRestr(:,1)==ii & signRestr(:,5)==2));

    end
    
end

totRestr = f + s; % Total number of restrictions on each column of Q
    
% Construct indicator for restrictions on one column of Q only.
oneColQ = (sum(totRestr > 0) <= 1); 

if oneColQ == 1 && whichAlgo == 3

    fprintf('\nRestrictions apply to one column of Q.'); 
    fprintf('\nComputing bounds of identified set using analytical approach of Gafarov, Meier and Montiel-Olea (2018)...\n');
    
elseif whichAlgo == 3 && oneColQ~=1
    
    fprintf('\nRestrictions apply to multiple columns of Q.\n');
    error('Approach of Gafarov, Meier and Montiel-Olea (2018) not applicable!');
    
elseif whichAlgo == 1
    
    fprintf('\nComputing bounds of identified set using Algorithm 1...\n');
    
elseif whichAlgo == 2
    
    fprintf('\nComputing bounds of identified set using Algorithm 2...\n');
    
end

%% Assess convexity of identified set using Proposition 3.
%convInfo.propCase = 0;

%@YZ2025JIMF Yoshida suppress this section
%if jshock == 1 && f(1) < n-1 
    
    % If equality restrictions satisfy Proposition 3 (I)(i).
%    convInfo.propCase = 1;

%elseif jshock >= 2 && ...
%        all(f(1:jshock-1) < n-(1:(jshock-1))') 
        
    % If equality restrictions satisfy Proposition 3 (I)(ii).
%    convInfo.propCase = 2;
    
%elseif jshock >= 2
    
    convInfo.propCase = 3;
    
%end

% If equality restrictions only and restrictions satisfy Proposition
% 3(I)(i) or (ii), IS is convex.
if mZ >= 0 && mS == 0
    
    if convInfo.propCase == 1 || convInfo.propCase == 2
        
        convInfo.isConvex = 1;
        propStr = {'i','ii'};
        fprintf('\nIdentified set is convex based on Proposition 3(I)%s.\n',...
            propStr{convInfo.propCase});
        
    elseif convInfo.propCase == 3
        
        fprintf('\nChecking convexity of identified set using Proposition 3(I)(iii) at each draw of phi.\n');
        
    end
    
end

% If sign restrictions only on the jshock-th column of Q.
if mS > 0 && all(s(setdiff(1:n,jshock))==0)
    
    if convInfo.propCase == 1 || convInfo.propCase == 2
        
        convInfo.propCase = 4;
        fprintf('\nChecking convexity of identified set using Proposition 3(II)(iv) at each draw of phi.\n');
        
    elseif convInfo.propCase == 3
        
        convInfo.propCase = 5;
        fprintf('\nChecking convexity of identified set using Proposition 3(II)(v) at each draw of phi.\n');

    end
    
end

end