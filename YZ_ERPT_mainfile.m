% The original code is "MP_Ext_mainfile.m" from "Identification and Inference
% under Narrative Restrictions" by Raffaella Giacomini, Toru Kitagawa, and
% Matthew Read, 2021, arXiv:2102.06456
% 
% The code is modified by Yushi Yoshida for the
% research paper, "Can exchange rate pass-throughs be perverse? 
% A robust multiple-prior Bayesian SVAR approach," Yushi Yoshida and
% Weiyang Zhai, available online Feb 2025, Journal of International Money and
% Finance, 103312
%
% The modified codes are indicated by the signature, @YZ2025JIMF
% The following 8 codes in auxFunctions are also modified
% (i)approximationBoundsYZ.m (ii)mainfileYZ.m (iii)intialComp_YZ.m
% (iv)credibleRegionERPT.m (v)HighestPoseteriorDensityERPT.m
% (vi)genVMA_YZ.m (vii)drawQ0eq.m (viii)drawQeq.m
% Other 6 codes in auxFunctions are not modified.

%%
% Estimate SVAR via Bayesian methods. Consider both a
% uniform prior over Q|phi and the prior-robust approach.


clear variables
close all
clc

oldFolder = pwd;
cd ..
addpath([oldFolder,'/auxFunctions']);
cd(oldFolder);

%% @YZ2025JIMF Choose narrative restrictions
% nr_yz =1 for Monetary shock in March 2002
% =2 for Monetary shock in Nov 2008, 
% =3 for Transitory Global shock in Nov 2008
% =4 for Supply shock in Mar 2011
% =5 for ER shock in Jan 2013
% =6 for Transitory Global shock in Mar 2022
% =7 for two shocks Mar 2011 & Jan 2013
nr_yz=4;

%% Import data.
% @YZ2025JIMF L30-L42
% The variables are:
% 1) GDP 2) CPI 3) shadow rate 4) exchange rate 5) import price 
% 6) world export price

[num,txt,raw] =xlsread('Database_final.xlsx','JPNdata_log_only_m');%seasonally unadjusted
dates = datenum(txt(2:end,1),'yyyy/mm/dd');
first = find(dates == datenum('01/01/95','dd/mm/yy'));
last = size(num,1)-find(dates == datenum('01/07/23','dd/mm/yy'));
DateNumber = datenum(txt(2+first:end-last,1),'dd/mm/yy');
data = num(first:end-last,:);
Detrended =[1,2,3,4,5,6];
%[Trend,HPData] = hpfilter(data(:,1:6),[14400,14400,Inf,14400,14400,14400]); % MATLB2022 version detrend all variables
[Trend,HPData] = hpfilter(data(:,1:6),Smoothing=[14400,14400,Inf,14400,14400,14400]); % detrend all variables
data(:,Detrended) = HPData;

% @YZ2025JIMF, re-ordered to be consistent with Definition 3 of GK (2018)
data=data(:,[2 3 4 1 6 5]);

% Re-order variables so that the federal funds rate is ordered first (the 
% monetary policy shock is the shock of interest).
%data = 100*data(:,[end 1:end-1]);

% Declare any exogenous variables (other than the constant; leave empty 
% if none).
exog = [];

% Set variable names for plotting figures.
%@YZ2025JIMF
varnames = {'CPI','Shadow Rate','Exchange Rate','JPN GDP','World Export Price','JPN Import Price'};

%@YZ2025JIMF for multiple types of structural shocks, re-ordered to match
%the ordering of variables
shocknames={'Demand shock', 'Monetary policy shocsk','Exogenous ER shock','Supply shock', 'Transitory global shock','Persistent global shock'};

% Set units for impulse responses. 0 = %, 1 = ppt, 2 = bps.
%irfUnits = [2 0 0 0 0 0];
irfUnits = [0 2 0 0 0 0];

%% Options.
opt.p = 6; %@YZ2025JIMF, origninally 12 No. of lags in VAR
opt.const = 1; % const = 1 if constant in VAR, = 0 otherwise
opt.ivar = [1 2 3 4 5 6];  % Indices of variables of interest
opt.cumIR = []; % Indices of variables for which cumulative impulse response should be used
opt.jshock = [1 2 3 4 5 6]; %@YZ2025JIMF, originally =1, Indices of structural shocks of interest
opt.H = 60; % Terminal horizon for impulse responses
opt.phiDraws = 1000; % No. of draws from posterior of phi with non-empty identified set
opt.aalpha = 0.68; % Credibility level for credible intervals
opt.maxDraw = 10000;%@YZ2025JIMF for fast computation onlytemporary, originally 100,000 % Max no. of draws from O(n) when checking empty identified set
opt.KTilde = 1000; %@YZ2025JIMF, No. of draws of Q used in Algorithm 1, original=10,000
opt.dispIter = 1; % Print number of draws remaining every dispIter draws
opt.gridLength = 1000; % Size of grid used when computing credible intervals

%rng(126453); % Set seed for random number generator

%% @YZ2025JIMF   Zero restrictions
% This part of codes is implanted from "Robust Bayesian inference for 
% set-identified models," Econometrica, 89(4), 1519-155 by R. Giancomini 
% and T. Kitagawa

% Each row of eqRestr contains a quartet (i,j,l,t) representing a 
% particular equality restriction, where t is the type of restriction:
% t = 1: the (ij)th element of A0 is zero
% t = 2: the (ij)th element of A0^(-1) is zero
% t = 3: the (ij)th element of A_l is zero
% t = 4: the (ij)th element of the long-run cumulative IR (LRCIR) is zero 
% The value of l does not matter except when t = 3.

%@YZ2025JIMF, Note that the order of variables are changed in line 55 of this code!!

restr.eqRestr = [1 5 0 1; % A0_15 (shock, variable)= 0
                2 5 0 1;
                3 5 0 1;
                4 5 0 1;
%                4 5 0 4;% LRCIR(vairable,shock) = 0
%                4 3 0 4;
%                4 2 0 4;
%                4 1 0 4;
                5 3 0 4;
                5 2 0 4;
                5 1 0 4;
                5 4 0 4]; 

%% Input identifying (traditional sign and narrative sign restrictions)

% Each row of shockSignRestr contains a vector (i,s,t) representing the
% shock-sign restriction that the ith shock in period t is >= 0 (s = 1) or
% <= 0 (s = -1). t is equal to the index of the observation.
% restr.shockSignRestr = []; % No shock-sign restrictions
%restr.shockSignRestr = [1 1 find(dates == datenum(yyyy,mm,dd))]; 

if nr_yz==1 
%@YZ2025JIMF Introduction of QE in March 2001 is negative.
restr.shockSignRestr = [2 -1 find(dates == datenum(2001,03,01))];
elseif nr_yz==2
%@YZ2025JIMF MP response to the Leman shock in Nov 2008 is negative.
restr.shockSignRestr = [2 -1 find(dates == datenum(2008,11,01))];
elseif nr_yz==3
%@YZ2025JIMF WEP(Transitory Global shock) response to the Leman shock in Nov 2008 is negative.
restr.shockSignRestr = [5 -1 find(dates == datenum(2008,11,01))];
elseif nr_yz==4
%@YZ2025JIMF The East Japan Great Earthquake supply shock in Mar 2011 is
%negative.
restr.shockSignRestr = [4 -1 find(dates == datenum(2011,03,01))];
elseif nr_yz==5
%@YZ2025JIMF Abe adiministration induced negative (depreciation) shock in
%Jan 2013.
restr.shockSignRestr = [3 -1 find(dates == datenum(2013,01,01))];
elseif nr_yz==6
%@YZ2025JIMF this case has one lag, 
restr.shockSignRestr = [5 1 find(dates == datenum(2022,03,01))];
elseif nr_yz==7
%@YZ2025JIMF two simultaneous restrictions
restr.shockSignRestr = [4 -1 find(dates == datenum(2011,03,01));
                        3 -1 find(dates == datenum(2013,01,01))];
end

% Each row of hdSignRestr contains a vector (i,j,t,h,k,s) representing a
% narrative restriction on the historical decomposition. 
% k represents the class of restriction (Type A or Type B) and s represents
% the type of restriction within the class.
% k = 1: Type A restriction
%  - s = 1:  The jth shock is the 'most important contributor' to the 
%   change in the ith variable between periods t and  t+h.
%  - s = -1: The jth shock is the 'least important contributor' to the
%   change in the ith variable between periods t and t+h.
% k = 2: Type B restriction
%  - s = 1: The jth shock is the 'overwhelming contributor' to the change
%   in the ith variable between periods t and t+h.
% - s = -1: The jth shock is a 'negligible contributor' to the change in
%   the ith variable between periods t and t+h.
% restr.hdSignRestr = []; % No restrictions on historical decomposition

if nr_yz==1
%@YZ2025JIMF Introduction of QE in March 2001 on ER is overwhelming.
restr.hdSignRestr = [3 2 find(dates == datenum(2001,03,01)) 0 2 1];
elseif nr_yz==2
%@YZ2025JIMF MP shock in Nov 2008 on ER is overwhelming.
restr.hdSignRestr = [3 2 find(dates == datenum(2008,11,01)) 0 2 1];
elseif nr_yz==3
%@YZ2025JIMF Transitory global shock in Nov 2008 on ER is overwhelming.
restr.hdSignRestr = [3 5 find(dates == datenum(2008,11,01)) 0 2 1];
elseif nr_yz==4
%@YZ2025JIMF Supply shock in Mar 2011 on GDP is overwhelming.
restr.hdSignRestr = [4 4 find(dates == datenum(2011,03,01)) 0 2 1];
elseif nr_yz==5
%@YZ2025JIMF ER shock in Jan 2013 on ER is overwhelming.
restr.hdSignRestr = [3 3 find(dates == datenum(2013,01,01)) 0 2 1];
elseif nr_yz==6
%@YZ2025JIMF Trnasitory Global shock in Apr 2022 on import price is overwhelming. 
restr.hdSignRestr = [6 5 find(dates == datenum(2022,04,01)) 0 2 1];
elseif nr_yz==7
%@YZ2025JIMF two simultanous restrictions, overwhelming
restr.hdSignRestr = [4 4 find(dates == datenum(2011,03,01)) 0 2 1;
                     3 3 find(dates == datenum(2013,01,01)) 0 2 1];
end

% Each row of signRestr contains a vector (i,j,h,s,t) representing a
% 'traditional' sign restriction, where t is the type of restriction:
% t = 1: the impulse response of the ith variable to the jth shock at the 
% hth horizon is nonnegative (s = 1) or nonpositive (s = -1).
% t = 2: the (ij)th element of A0 is nonnegative (s = 1) or nonpositive 
% (s = -1). 
% signRestr = []; % No sign restrictions on impulse responses
restr.signRestr = ...
      [1 4 0 -1 1; % Response of CPI to supply shock on impact is nonpositive
%      1 4 1 -1 1; % As above after one month
%      1 4 2 -1 1; % As above after two months
%      1 4 3 -1 1; % As above after three months
      4 4 0 1 1; % Response of Ind prod to supply shock on impact is nonnegative
%      4 4 1 1 1; % As above after one month
%      4 4 2 1 1; % As above after two months
%      4 4 3 1 1; % As above after three months
       4 1 0 1 1; %
%       4 1 1 1 1; % as above after one month
       1 1 0 1 1; % response of CPI to demand shock on impact is nonnegative
%       1 1 1 1 1; % as above after one month
       2 1 0 1 1; %
%       2 1 1 1 1; % as above after one month
       3 1 0 1 1; %
%       3 1 1 1 1; % as above after one month
       4 2 0 -1 1; % 
       1 2 0 -1 1; % 
       2 2 0 1 1;
       3 2 0 1 1;
       1 3 0 -1 1;
       2 3 0 -1 1;
       3 3 0 1 1; % Response of c to monetary policy shock is nonpositive
       5 5 0 1 1; % As above after one month, etc
       6 6 0 1 1; 
       5 6 0 1 1]; % As above after one month, etc   

%% Conduct (robust) posterior inference.
mainfileYZ;

%@YZ2025JIMF
%@YZ2025JIMF calculate ERPT for a single prior post
%@YZ2025JIMF ERPT is defined as the ratio of cumulated IRs.
        A_yy1=ones(opt.H+1); %@YZ2025JIMF 
        A_yy2=triu(A_yy1); %@YZ2025JIMF a matrix to construct cumulted IRs.
    for j_shock=1:6
    rSinglePriorPostERPT(:,:,j_shock)=(rSinglePriorPost(:,:,1,j_shock)*A_yy2)./(rSinglePriorPost(:,:,3,j_shock)*A_yy2);
    end

    
% Compute posterior probabilities of ERPT falling after two years.
% Compute probability under single prior.

posteriorProb = mean(rSinglePriorPostERPT(:,25,1) <0);
fprintf('\nPosterior probability that demand shock ERPT falls below zero after two years is %1.4g\n',...
    posteriorProb);

% Compute lower posterior probability.
lowerProb = mean(rMaxERPT(:,25,1) < 0);
fprintf('\nLower posterior probability that demand shock ERPT falls after two years is %1.4g\n',...
    lowerProb);

%@YZ2025JIMF ... the output file name is manually modified
save('2NSR11Mar_13Jan_End23Jul_results.mat');

%% Plot impulse responses.
%@YZ2025JIMF
%@YZ2025JIMF for multi-kind shocks (originally only for a single kind of shock)
%{
curFolder = pwd;
cd('Figures');
for jj = 1:nj
for ii = 1:ni
    
    figure;
    h1 = plot(0:opt.H,meanlb(:,ii,jj),'color','k','LineStyle','--');
    hold on;
    plot(0:opt.H,meanub(:,ii,jj),'color','k','LineStyle','--');
    h2 = plot(0:opt.H,credlb(:,ii,jj),'-k','LineWidth',2);
    plot(0:opt.H,credub(:,ii,jj),'-k','LineWidth',2); % fixed on Nov16,2023
    h3 = plot(0:opt.H,rSinglePriorMean(:,ii,jj),'ko','MarkerFaceColor','k');
    h4 = plot(0:opt.H,hpdlb(:,ii,jj),'--k','LineWidth',2);
    plot(0:opt.H,hpdub(:,ii,jj),'--k','LineWidth',2);
    
    if irfUnits(opt.ivar(ii)) == 0
        ylabel('%');
    elseif irfUnits(opt.ivar(ii)) == 1
        ylabel('ppt');
    elseif irfUnits(opt.ivar(ii)) == 2
        ylabel('bps');
    end
    
    title(sprintf('%s',varnames{opt.ivar(ii)},shocknames{opt.jshock(jj)}));
    xlabel('Horizon (months)');
    xticks(0:12:60);
    
    for hh = 0:opt.H
        line([hh hh],[meanlb(hh+1,ii,jj) meanub(hh+1,ii,jj)],'color','k');
    end

    line([0 opt.H],[0 0],'color','black','LineStyle',':');
    
    Ax = gca;
    Ax.FontSize = 14;
    opt.g_title=strcat(varnames{opt.ivar(ii)},' by ',shocknames{opt.jshock(jj)});
    %print(sprintf('%s_Mar01',varnames{opt.ivar(ii)},shocknames{opt.jshock(jj)}),'-dpng');
    print(sprintf('%s_MAR11NSR(SS_IP)',opt.g_title),'-dpng');

end
end
cd(oldFolder);
%}
%% compute credible regions for ERPT
%@YZ2025JIMF
% Compute robustified credible regions.
[ERPTcredlb,ERPTcredub] = credibleRegionERPT(rMinERPT,rMaxERPT,opt);
% Compute highest posterior density (HPD) interval under single prior.
        A_yy1=ones(opt.H+1);
        A_yy2=triu(A_yy1); % matrix to construct cumulted IRs.
for j_shock=1:6
rSinglePriorPostERPT(:,:,j_shock)=(rSinglePriorPost(:,:,1,j_shock)*A_yy2)./(rSinglePriorPost(:,:,3,j_shock)*A_yy2);
end
[hpdlbERPT,hpdubERPT] = highestPosteriorDensityERPT(rSinglePriorPostERPT,opt);
%% plot cumulated ERPT impulse responses
%@YZ2025JIMF
curFolder = pwd;
cd('Figures');
        A_yy1=ones(opt.H+1);
        A_yy2=triu(A_yy1); % matrix to construct cumulted IRs.

for j_shock=1:6
figure;
%ERPTlb = meanlb(:,1,:)./meanlb(:,3,:);
%ERPTub = meanub(:,1,:)./meanub(:,3,:);
ERPT=(rSinglePriorPost(:,:,1,j_shock)*A_yy2)./(rSinglePriorPost(:,:,3,j_shock)*A_yy2);
    for ii=1:opt.phiDraws
    h1 = plot(0:opt.H,ERPT(ii,:,1,1),'color','k','LineStyle','--');
    hold on;
    end
    
    plot(0:opt.H,meanERPTub(:,j_shock),'color','b','LineStyle','-',LineWidth=2);
    plot(0:opt.H,meanERPTlb(:,j_shock),'color','r','LineStyle','-',LineWidth=2);

    title(sprintf('%s','ERPT',shocknames{opt.jshock(j_shock)}));
    xlabel('Horizon (months)');
    xticks(0:12:60);
% credible bounds, adjusted for multiple shocks    
    h2 = plot(0:opt.H,ERPTcredlb(:,j_shock),'-k','LineWidth',2);
    plot(0:opt.H,ERPTcredub(:,j_shock),'-k','LineWidth',2); % fixed on Nov16,2023
% higher posibility density  
    h4 = plot(0:opt.H,hpdlbERPT(:,j_shock),'-g','LineWidth',2);
    plot(0:opt.H,hpdubERPT(:,j_shock),'-g','LineWidth',2);


    for hh = 0:opt.H
        line([hh hh],[meanERPTlb(hh+1,j_shock) meanERPTub(hh+1,j_shock)],'color','k');
    end
    opt.g_title=strcat(string(datetime('today','Format','uuuu-MM-dd')),'ERPT by ',shocknames{opt.jshock(j_shock)});
    print(sprintf('%s_MAR11NSR(SS_IP)',opt.g_title),'-dpng');    
end
cd(curFolder);
