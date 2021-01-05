%% load data
metaSignal         = load_csv('sFile','metaSignal.csv');

sTable                  = 'emg';
idSignalList            = getMeta(metaSignal,qry('sTable',sTable));
idSignalList([13,14])   = [];
sSignalList             = getMeta(metaSignal,qry('idSignal',idSignalList),'sSignal');
idSubjectList           = 2:11; % all subjects 2:11
idTrialTypeList         = [4,5,6];
sTrialTypeList          = {'Control','Resistive','Assistive'};
cCol                    = 'krb';

sTitle = ['DynMEPgain_allTrialType_allSignal']; % only gain is analyzed
load([sTitle,'.mat'])
%% plotting Fig. 7A and 7B    
for iSignal = 1:length(idSignalList)
    sTitle1 = ['DistribMEPgain_allTrialType_allSubject_',sSignalList{iSignal}];
    sTitle2 = ['DiffMEPgain_allTrialType_allSubject_',sSignalList{iSignal}];
    [hFig, hPlot]  = setPlot('nRow',1,'nCol',3,'sAnnotation',sTitle1);
    [hFig2, hPlot2]  = setPlot('nRow',1,'nCol',10,'sAnnotation',sTitle2);
    for iTrialType = 1:numel(idTrialTypeList)
        for iSubject = 1:numel(idSubjectList)
            nData11 = nMedian1{iSubject,iSignal,iTrialType}; % early bins 1,2,3
            nData21 = nMedian2{iSubject,iSignal,iTrialType}; % late bins 5,6
            nH1    = Hall(iSubject,iSignal,iTrialType);
            if nH1
                scatter(hPlot(1,iTrialType),...
                    iSubject*ones(length(nData11),1)-0.1,nData11',20,[0,0,0],'filled')
                if iSubject == 1
                    hold(hPlot(1,iTrialType),'on')
                end
                scatter(hPlot(1,iTrialType),...
                    iSubject*ones(length(nData21),1)+0.1,nData21',20,[0,0,1],'filled')
            else
                scatter(hPlot(1,iTrialType),...
                    iSubject*ones(length(nData11),1)-0.1,nData11',20,[0,0,0])
                if iSubject == 1
                    hold(hPlot(1,iTrialType),'on')
                end
                scatter(hPlot(1,iTrialType),...
                    iSubject*ones(length(nData21),1)+0.1,nData21',20,[0,0,1])
            end
            nDiff = [];
            for iRep = 1:numel(nData11)
                nDiff = [nDiff, nData11(iRep)-nData21];
            end
            histogram(hPlot2(iSubject),nDiff,'Orientation','horizontal',...
                'FaceColor',cCol(iTrialType))
            if iTrialType==1
                hold(hPlot2(iSubject),'on')
            end
        end
        ylabel(hPlot(1,iTrialType),sTrialTypeList{iTrialType})
        set(hPlot(1,iTrialType),'XTick',1:10)
        xlim(hPlot(1,iTrialType),[0.5,10.5])
        xlabel(hPlot(1,iTrialType),'idSubject')
    end
    nYLim1 = ylim(hPlot(1));
    nYLim2 = ylim(hPlot(2));
    nYLim3 = ylim(hPlot(3));
    nYlim = [min([nYLim1(1),nYLim2(1),nYLim3(1)]),max([nYLim1(2),nYLim2(2),nYLim3(2)])];
    ylim(hPlot(1),nYlim)
    ylim(hPlot(2),nYlim)
    ylim(hPlot(3),nYlim)
end
%% modulation test across individuals: Fig. 7C
nMedianAll = [];
for iSignal = 1:length(idSignalList)
    for iTrialType = 1:numel(idTrialTypeList)
        for iSubject = 1:numel(idSubjectList)
            nMedianAll(iSubject,iTrialType,1,iSignal) = nanmedian(nMedian1{iSubject,iSignal,iTrialType}); % early bins 1,2,3
            nMedianAll(iSubject,iTrialType,2,iSignal) = nanmedian(nMedian2{iSubject,iSignal,iTrialType}); % late bins 5,6v
        end
    end
end
% Sidak correction: 1-(1-0.05)^(1/3)
sTitle2 = ['DiffMedianMEPgain_allTrialType_allSubject_allSignal'];
[hFig2, hPlot2]  = setPlot('nRow',3,'nCol',4,'sAnnotation',sTitle2);
for iSignal = 1:length(idSignalList)
    nPtt = [];
    for iTrialType = 1:numel(idTrialTypeList)
        nDiff1 = nMedianAll(:,iTrialType,1,iSignal);
        nDiff2 = nMedianAll(:,iTrialType,2,iSignal);
        [h,P,~,stat] = ttest((nDiff1-nDiff2)); % Sidak correction for 3 trial types, Alpha = 0.017
        nPtt(iTrialType) = stat.tstat;
        histogram(hPlot2(iSignal),nDiff1-nDiff2,'Orientation','horizontal',...
            'FaceColor',cCol(iTrialType))
        if iTrialType==1
            hold(hPlot2(iSignal),'on')
        end
    end
    title(hPlot2(iSignal),sSignalList{iSignal})
    display([sSignalList{iSignal},' pC',num2str(nPtt(1)),' pR',num2str(nPtt(2)),' pA',num2str(nPtt(3))])
end
%% comparing across tasks for Fig. 8 & 9
load('nAngleAll.mat')
nAngle23early   = (nAngleAll(:,2,:,1)-nAngleAll(:,3,:,1));
nAngle13early   = (nAngleAll(:,1,:,1)-nAngleAll(:,3,:,1));
nAngle12early   = (nAngleAll(:,1,:,1)-nAngleAll(:,2,:,1));
nAngle23late    = (nAngleAll(:,2,:,2)-nAngleAll(:,3,:,2));
nAngle13late    = (nAngleAll(:,1,:,2)-nAngleAll(:,3,:,2));
nAngle12late    = (nAngleAll(:,1,:,2)-nAngleAll(:,2,:,2));
load('nTorAll_GT.mat')
nTorAll_GT      = nTorAll;
nGT23early      = (nTorAll_GT(:,2,:,1)-nTorAll_GT(:,3,:,1));
nGT13early      = (nTorAll_GT(:,1,:,1)-nTorAll_GT(:,3,:,1));
nGT12early      = (nTorAll_GT(:,1,:,1)-nTorAll_GT(:,2,:,1));
nGT23late       = (nTorAll_GT(:,2,:,2)-nTorAll_GT(:,3,:,2));
nGT13late       = (nTorAll_GT(:,1,:,2)-nTorAll_GT(:,3,:,2));
nGT12late       = (nTorAll_GT(:,1,:,2)-nTorAll_GT(:,2,:,2));
load('nTorAll_DT.mat')
nTorAll_DT      = nTorAll;
nDT23early      = (nTorAll_DT(:,2,:,1)-nTorAll_DT(:,3,:,1));
nDT13early      = (nTorAll_DT(:,1,:,1)-nTorAll_DT(:,3,:,1));
nDT12early      = (nTorAll_DT(:,1,:,1)-nTorAll_DT(:,2,:,1));
nDT23late       = (nTorAll_DT(:,2,:,2)-nTorAll_DT(:,3,:,2));
nDT13late       = (nTorAll_DT(:,1,:,2)-nTorAll_DT(:,3,:,2));
nDT12late       = (nTorAll_DT(:,1,:,2)-nTorAll_DT(:,2,:,2));
load('nTorAll_MT.mat')
nTorAll_MT      = nTorAll;
nMT23early      = (nTorAll_MT(:,2,:,1)-nTorAll_MT(:,3,:,1));
nMT13early      = (nTorAll_MT(:,1,:,1)-nTorAll_MT(:,3,:,1));
nMT12early      = (nTorAll_MT(:,1,:,1)-nTorAll_MT(:,2,:,1));
nMT23late       = (nTorAll_MT(:,2,:,2)-nTorAll_MT(:,3,:,2));
nMT13late       = (nTorAll_MT(:,1,:,2)-nTorAll_MT(:,3,:,2));
nMT12late       = (nTorAll_MT(:,1,:,2)-nTorAll_MT(:,2,:,2));

nZ1=[]; nZ2=[];
for iSignal = 1:12
    if iSignal<6 || iSignal==7 % include BicL and TriLo into shoulder regr.
        iSignalKD = 1;  % shoulder
    elseif iSignal==6 || iSignal==8
        iSignalKD = 2;% elbow
    else
        iSignalKD = 3; % wrist, include Br into wrist regression
    end
%     if iSignal<5
%         iSignalKD = 1; % shoulder
%     elseif iSignal>9
%         iSignalKD = 3; % wrist
%     else
%         iSignalKD = 2; % elbow
%     end
    for iSubject = 1:numel(idSubjectList)
        
        nData11 = nMedian1{iSubject,iSignal,1}; % early bins 1,2,3
        nData12 = nMedian1{iSubject,iSignal,2};
        nData13 = nMedian1{iSubject,iSignal,3};
        
        nData21 = nMedian2{iSubject,iSignal,1}; % late bins 5,6,7
        nData22 = nMedian2{iSubject,iSignal,2};
        nData23 = nMedian2{iSubject,iSignal,3};
        
        nData1112 = nanmedian(nData11)-nanmedian(nData12); % early bins 1-2
        nData1113 = nanmedian(nData11)-nanmedian(nData13);
        nData1213 = nanmedian(nData12)-nanmedian(nData13);
        
        
        nData2122 = nanmedian(nData21)-nanmedian(nData22); % late bins 1-2
        nData2123 = nanmedian(nData21)-nanmedian(nData23);
        nData2223 = nanmedian(nData22)-nanmedian(nData23);
        
        % Control - Resistive
        [H1,P1,kStat1]  = kstest2(nData11',nData12','Alpha',0.0051);  % early bins 2,3
        nZ1             = [nZ1;[nData1112,nAngle12early(iSubject,1,iSignalKD),...
            nGT12early(iSubject,1,iSignalKD),nDT12early(iSubject,1,iSignalKD),...
            nMT12early(iSubject,1,iSignalKD),12,H1,iSignal]];
        [H2,P2,kStat2]  = kstest2(nData21',nData22','Alpha',0.0051); % late bins 5,6
        nZ2             = [nZ2;[nData2122,nAngle12late(iSubject,1,iSignalKD),...
            nGT12late(iSubject,1,iSignalKD),nDT12late(iSubject,1,iSignalKD),...
            nMT12late(iSubject,1,iSignalKD),12,H2,iSignal]];
        
        % Resistive - Assistive
        [H1,P1,kStat1]  = kstest2(nData12',nData13','Alpha',0.0051);  % early bins 2,3
        nZ1             = [nZ1;[nData1213,nAngle23early(iSubject,1,iSignalKD),...
            nGT23early(iSubject,1,iSignalKD),nDT23early(iSubject,1,iSignalKD),...
            nMT23early(iSubject,1,iSignalKD),23,H1,iSignal]];
        [H2,P2,kStat2]  = kstest2(nData22',nData23','Alpha',0.0051); % late bins 5,6
        nZ2             = [nZ2;[nData2223,nAngle23late(iSubject,1,iSignalKD),...
            nGT23late(iSubject,1,iSignalKD),nDT23late(iSubject,1,iSignalKD),...
            nMT23late(iSubject,1,iSignalKD),23,H2,iSignal]];
        
        % Control- Assistive
        [H1,P1,kStat1]  = kstest2(nData11',nData13','Alpha',0.0051);  % early bins 2,3
        nZ1             = [nZ1;[nData1113,nAngle13early(iSubject,1,iSignalKD),...
            nGT13early(iSubject,1,iSignalKD),nDT13early(iSubject,1,iSignalKD),...
            nMT13early(iSubject,1,iSignalKD),13,H1,iSignal]];
        [H2,P2,kStat2]  = kstest2(nData21',nData23','Alpha',0.0051); % late bins 5,6
        nZ2             = [nZ2;[nData2123,nAngle13late(iSubject,1,iSignalKD),...
            nGT13late(iSubject,1,iSignalKD),nDT13late(iSubject,1,iSignalKD),...
            nMT13late(iSubject,1,iSignalKD),13,H2,iSignal]];
        
    end
end
%%
cMap        = jet(12);
sTitle      = ['DiffBetweenTasksMEPgain_allTrialType_BicLTriLo_Sho_allSubject'];

% bSignal1    = nZ1(:,7)==1 & nZ1(:,8)==2; %AD
% bSignal2    = nZ2(:,7)==1 & nZ2(:,8)==2;
% iSignalList = 2;

% bSignal1    = nZ1(:,7)==1 & (nZ1(:,8)<5 & nZ1(:,8)~=2); %PecPDTM
% bSignal2    = nZ2(:,7)==1 & (nZ2(:,8)<5 & nZ2(:,8)~=2);
% iSignalList = [1,3:4];
% 
% bSignal1   = nZ1(:,7)==1 & nZ1(:,8)>4 & nZ1(:,8)<9; %BicTri
% bSignal2    = nZ2(:,7)==1 & nZ2(:,8)>4 & nZ2(:,8)<9;
% iSignalList = 5:8;

% bSignal1       = nZ1(:,7)==1 & nZ1(:,8)>8; %Wrist
% bSignal2       = nZ2(:,7)==1 & nZ2(:,8)>8;
% iSignalList    = 9:12;

% bSignal1        = nZ1(:,7)==1; % all muscles
% bSignal2        = nZ2(:,7)==1;
% iSignalList     = 1:12;

bSignal1    = nZ1(:,7)==1 & (nZ1(:,8)==5 | nZ1(:,8)==7); %BicLTriLo
bSignal2    = nZ2(:,7)==1 & (nZ2(:,8)==5 | nZ2(:,8)==7);
iSignalList = [5,7];


nMax1 = max(nZ1(bSignal1,:)) - min(nZ1(bSignal1,:));
nMax2 = max(nZ2(bSignal2,:)) - min(nZ2(bSignal2,:));

% repeated measures 1-(1-0.05)^(1/4) = 0.0127
[hFig, hPlot]  = setPlot('nRow',2,'nCol',4,'sAnnotation',sTitle);
[B,R,P,hTitle] = getRegress(nZ1(bSignal1,2)/nMax1(2),...
    nZ1(bSignal1,1)/nMax1(1),'bPlot' ,1,'hPlot' ,hPlot(1,1));
hold(hPlot(1,1),'on')
[B,R,P,hTitle] = getRegress(nZ1(bSignal1,5)/nMax1(5),...
    nZ1(bSignal1,1)/nMax1(1),'bPlot' ,1,'hPlot' ,hPlot(1,2));
hold(hPlot(1,2),'on')
%% testing power of regression
x = nZ1(bSignal1,5)/nMax1(5);
y = nZ1(bSignal1,1)/nMax1(1);
[B,bint,r,rint,stats] = regress(y,[ones(length(x),1),x]); % test model y=B1+B2*x+err
mu0 = B(2); % significant slope
x = nZ2(bSignal2,5)/nMax2(5);
y = nZ2(bSignal2,1)/nMax2(1);
[B,bint,r,rint,stats] = regress(y,[ones(length(x),1),x]); % test model y=B1+B2*x+err
Bsd = std(r); % residuals from insignificant slope
sig = Bsd;
N = numel(r);
alpha = 0.05;
conf = 1-alpha;
cutoff = norminv(conf, mu0, sig/sqrt(N));
x = [linspace(0,cutoff), linspace(cutoff,2)];
y = normpdf(x,mu0,sig/sqrt(N));
figure
h1 = plot(x,y);
xhi = [cutoff, x(x>=cutoff)];
yhi = [0, y(x>=cutoff)];
patch(xhi,yhi,'b');
title(['Distribution of sample mean, N=',num2str(N)]);
xlabel('Sample mean');
ylabel('Density');
text(0.8,2,sprintf('Reject if mean>%.4g\nProb = 0.05',cutoff),'Color','b');
mu1 = mu0+0.2637;
y2 = normpdf(x,mu1,sig/sqrt(N));
h2 = line(x,y2,'Color','r');
yhi = [0, y2(x>=cutoff)];
patch(xhi,yhi,'r','FaceAlpha',0.25);
P = 1 - normcdf(cutoff,mu1,sig/sqrt(N));
text(1.2,5,sprintf('Reject if T>%.4g\nProb = %.2g',cutoff,P),'Color',[1 0 0]);
legend([h1 h2],'Null hypothesis','Alternative hypothesis');
%%
getRegress(nZ1(bSignal1,3)/nMax1(3),...
    nZ1(bSignal1,1)/nMax1(1),'bPlot' ,1,'hPlot' ,hPlot(1,3));
hold(hPlot(1,3),'on')
getRegress(nZ1(bSignal1,4)/nMax1(4),...
    nZ1(bSignal1,1)/nMax1(1),'bPlot' ,1,'hPlot' ,hPlot(1,4));
hold(hPlot(1,4),'on')
% late
getRegress(nZ2(bSignal2,2)/nMax2(2),...
    nZ2(bSignal2,1)/nMax2(1),'bPlot' ,1,'hPlot' ,hPlot(2,1));
hold(hPlot(2,1),'on')
[B,R,P,hTitle] = getRegress(nZ2(bSignal2,5)/nMax2(5),...
    nZ2(bSignal2,1)/nMax2(1),'bPlot' ,1,'hPlot' ,hPlot(2,2));
hold(hPlot(2,2),'on')
getRegress(nZ2(bSignal2,3)/nMax2(3),...
    nZ2(bSignal2,1)/nMax2(1),'bPlot' ,1,'hPlot' ,hPlot(2,3));
hold(hPlot(2,3),'on')
getRegress(nZ2(bSignal2,4)/nMax2(4),...
    nZ2(bSignal2,1)/nMax2(1),'bPlot' ,1,'hPlot' ,hPlot(2,4));
hold(hPlot(2,4),'on')
sLegend = [{''},{''}];
for iSignal = iSignalList
    bSignal1 = nZ1(:,7)==1 & nZ1(:,8)==iSignal;
    bSignal2 = nZ2(:,7)==1 & nZ2(:,8)==iSignal;
    % early
    plot(hPlot(1,1),nZ1(bSignal1,2)/nMax1(2),nZ1(bSignal1,1)/nMax1(1),'o','MarkerFaceColor',cMap(iSignal,:))
    plot(hPlot(1,2),nZ1(bSignal1,5)/nMax1(5),nZ1(bSignal1,1)/nMax1(1),'o','MarkerFaceColor',cMap(iSignal,:))
    plot(hPlot(1,3),nZ1(bSignal1,3)/nMax1(3),nZ1(bSignal1,1)/nMax1(1),'o','MarkerFaceColor',cMap(iSignal,:))
    plot(hPlot(1,4),nZ1(bSignal1,4)/nMax1(4),nZ1(bSignal1,1)/nMax1(1),'o','MarkerFaceColor',cMap(iSignal,:))
    % late
    plot(hPlot(2,1),nZ2(bSignal2,2)/nMax2(2),nZ2(bSignal2,1)/nMax2(1),'o','MarkerFaceColor',cMap(iSignal,:))
    plot(hPlot(2,2),nZ2(bSignal2,5)/nMax2(5),nZ2(bSignal2,1)/nMax2(1),'o','MarkerFaceColor',cMap(iSignal,:))
    plot(hPlot(2,3),nZ2(bSignal2,3)/nMax2(3),nZ2(bSignal2,1)/nMax2(1),'o','MarkerFaceColor',cMap(iSignal,:))
    plot(hPlot(2,4),nZ2(bSignal2,4)/nMax2(4),nZ2(bSignal2,1)/nMax2(1),'o','MarkerFaceColor',cMap(iSignal,:))
    
    sLegend = [sLegend,sSignalList{iSignal}];
end
ylabel(hPlot(1,1),'nDiff')
xlabel(hPlot(2,1),'Angle late')
ylabel(hPlot(2,1),'nDiff')
xlabel(hPlot(2,2),'MT')
xlabel(hPlot(2,3),'GT')
xlabel(hPlot(2,4),'DT')

legend(hPlot(2,4),sLegend)