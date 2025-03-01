% This program is coded by Yushi Yoshida for re-producing the
% ERPT impulse reponse functions satisfying restrictions.
% December 25, 2023
% Figures without gunshots, modified on March 3, 2024
% Multiple database, added on March 8, 2024

%% Output selections

clear;
clc;

oldFolder = pwd;
cd ..
addpath([oldFolder,'/resultFiles']);
cd(oldFolder);

%%% Single or Comparison
out_type= 2; % =1 (single), =2 (multiple)

%%% Gunshot option
gun_type=0; % =0 (No), =1 (Yes)


%%% change the result file number
numFile(1)= 7; % this is necesaary all the time, base results
numFile(2)= 10; % this is used only if out_type=2, overlaying

%% data files
%%% No NSR restriction
resultFileName(1)="NSRNo_End23Jul_results.mat";
resultFileName(2)="NSRNo_End20Jan_results.mat";

%%% Single NSR restriction
resultFileName(3)="NSR11Mar_End23Jul_results.mat";
resultFileName(4)="NSR11Mar_End20Jan_results.mat";
resultFileName(5)="NSR13Jan_End23Jul_results.mat";
resultFileName(6)="NSR13Jan_End20Jan_results.mat";
resultFileName(7)="NSR22Apr_End23Jul_results.mat";

%%% Two NSR restrictions
%resultFileName(8)="2NSR11Mar_13Jan_End23Jul_results.mat";

resultFileName(9)="2NSR11Mar_13Jan_End20Jan_results.mat";
resultFileName(10)="2NSR11Mar_13Jan_End23Jul_resultsNew.mat";
%%%%%%%

for i = 1:out_type
s.d(i)=load(resultFileName(numFile(i)));
end

%% Set working directories.

oldFolder = pwd;
cd ..
addpath([oldFolder,'/auxFunctions']);
cd(oldFolder);


%% compute credible regions for ERPT
for i=1:out_type
rMinERPT=s.d(i).rMinERPT;
rMaxERPT=s.d(i).rMaxERPT;
A_yy1=s.d(i).A_yy1;
A_yy2=s.d(i).A_yy2;
opt=s.d(i).opt
rSinglePriorPost=s.d(i).rSinglePriorPost;
% Compute robustified credible regions.
[ERPTcredlb,ERPTcredub] = credibleRegionERPT(rMinERPT,rMaxERPT,opt);
s.d(i).ERPTcredlb=ERPTcredlb;
s.d(i).ERPTcredub=ERPTcredub;
% Compute highest posterior density (HPD) interval under single prior.
        A_yy1=ones(opt.H+1);
        A_yy2=triu(A_yy1); % matrix to construct cumulted IRs.
    for j_shock=1:6
    rSinglePriorPostERPT(:,:,j_shock)=(rSinglePriorPost(:,:,1,j_shock)*A_yy2)./(rSinglePriorPost(:,:,3,j_shock)*A_yy2);
    end
    [hpdlbERPT,hpdubERPT] = highestPosteriorDensityERPT(rSinglePriorPostERPT,opt);
s.d(i).hpdlbERPT=hpdlbERPT;
s.d(i).hpdubERPT=hpdubERPT;
end
%% plot cumulated ERPT impulse responses
% created by Yoshida on Nov 16, 2023
% modified on Mar 8, 2024
curFolder = pwd;
cd('Figures');
        A_yy1=ones(opt.H+1);
        A_yy2=triu(A_yy1); % matrix to construct cumulted IRs.

for j_shock=1:6
for i=1:out_type
rSinglePriorPost=s.d(i).rSinglePriorPost;
meanERPTub=s.d(i).meanERPTub;    
meanERPTlb=s.d(i).meanERPTlb;
ERPTcredub=s.d(i).ERPTcredub;
ERPTcredlb=s.d(i).ERPTcredlb;
hpdubERPT=s.d(i).hpdubERPT;
hpdlbERPT=s.d(i).hpdlbERPT;
shocknames=s.d(i).shocknames;
if i==1
figure;
end
%ERPTlb = meanlb(:,1,:)./meanlb(:,3,:);
%ERPTub = meanub(:,1,:)./meanub(:,3,:);
ERPT=(rSinglePriorPost(:,:,1,j_shock)*A_yy2)./(rSinglePriorPost(:,:,3,j_shock)*A_yy2);

%%% Initializing a figue with horizontal line at zero 
     h1=line([0 opt.H],[0 0],'color','k','LineStyle','-');
     hold on;

% Gunshot IRFs     
    if out_type==1 & gun_type==1
        for ii=1:opt.phiDraws
        h1a = plot(0:opt.H,ERPT(ii,:,1,1),'color','k','LineStyle','--');
        hold on;
        end
    end

% the set of means
    if i==1
        plot(0:opt.H,meanERPTub(:,j_shock),'color','b','LineStyle','-.',LineWidth=1.5);
        plot(0:opt.H,meanERPTlb(:,j_shock),'color','b','LineStyle','-.',LineWidth=1.5);
    else
        plot(0:opt.H,meanERPTub(:,j_shock),'color',[0,0,1,0.3],'LineStyle','-.',LineWidth=2.5);
        plot(0:opt.H,meanERPTlb(:,j_shock),'color',[0,0,1,0.3],'LineStyle','-.',LineWidth=2.5);
    end

    title(sprintf('%s',shocknames{opt.jshock(j_shock)},' ERPT on CPI'));
    xlabel('Horizon (months)');
    xticks(0:12:60);
% credible bounds, adjusted for multiple shocks    
    if i==1
        h2 = plot(0:opt.H,ERPTcredlb(:,j_shock),'-k','LineWidth',2);
        plot(0:opt.H,ERPTcredub(:,j_shock),'-k','LineWidth',2); % fixed on Nov16,2023
    else
        h6 = plot(0:opt.H,ERPTcredlb(:,j_shock),'LineStyle','-','color',[0,0,0,0.3],'LineWidth',4);
        plot(0:opt.H,ERPTcredub(:,j_shock),'LineStyle','-','color',[0,0,0,0.3],'LineWidth',4); % fixed on Nov16,2023
    end
    % higher posibility density  
    if i==1
        h4 = plot(0:opt.H,hpdlbERPT(:,j_shock),'--r','LineWidth',2);
        plot(0:opt.H,hpdubERPT(:,j_shock),'--r','LineWidth',2);
    else
        h8 = plot(0:opt.H,hpdlbERPT(:,j_shock),'LineStyle','--','color',[1,0,0,0.3],'LineWidth',4);
        plot(0:opt.H,hpdubERPT(:,j_shock),'LineStyle','--','color',[1,0,0,0.3],'LineWidth',4);
    end


    for hh = 0:opt.H
        if i==1
        line([hh hh],[meanERPTlb(hh+1,j_shock) meanERPTub(hh+1,j_shock)],'color','k');
        else
        line([hh hh],[meanERPTlb(hh+1,j_shock) meanERPTub(hh+1,j_shock)],'color',[0.5,0.5,0.5,0.2],'LineWidth',4);
        end
    end
%hold on;
end
    opt.g_title=strcat('_ERPT by ',shocknames{opt.jshock(j_shock)});

    if out_type==1 && i==1
        saveas(gcf,strcat(resultFileName(numFile(1)),opt.g_title,'.png'));
        close;
    elseif out_type==2 && i==2
        saveas(gcf,strcat(resultFileName(numFile(1)),resultFileName(numFile(2)),opt.g_title,'.png'));
        close;
    end

end
cd(curFolder);
