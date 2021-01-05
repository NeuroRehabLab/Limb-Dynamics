%% means for paper between MEP magnitude and EMG, Fig. 6C & 6D

idTrialTypeList = [4:6];
sPath           = 'RSquared_MEParea_iMuscle_iTrialType_allSubjects';
sY              = 'MEP_area';
sSignalEMGList  =  {'Pectoralis','AntDelt','PostDelt','TeresMaj',...
    'TriLong','TriShort','BicepLong','BicepShort','Brachiorad',...
    'FlexCarpiRad','FlexCarpiUln','ExtCarpiRad'};
sSignalKinList  = {'Shoulder_Right1','Shoulder_Right1','Shoulder_Right1','Shoulder_Right1',...
    'Elbow_Right1','Elbow_Right1','Elbow_Right1','Elbow_Right1','Elbow_Right1',...
    'Wrist_Right1','Wrist_Right1','Wrist_Right1'};
clear nData nDataRat nDataGain
for iSignal = 1:numel(sSignalEMGList)
    sFileMat = [sPath,filesep,'RSquared_',sSignalEMGList{iSignal},'_',sSignalKinList{iSignal},'_',sY,'_allTrialType_allSubject.mat'];
    load(sFileMat)
    nData(:,iSignal) = nanmedian(nRsqured,2);
    nDataRat(1,:,iSignal) = nanmedian(nanmedian(nRatAll,1),3);
    nDataGain(1,:,iSignal) = nanmedian(nanmedian(nGainAll,1),3);
    nDiff(:,:,iSignal) = nRsqured;
end
nDataMean = prctile(nData,[25 50 75]);
nX = [1:12]';
nPlotMean = [nDataMean(1,:);nDataMean(2,:) - nDataMean(1,:);nDataMean(3,:) - nDataMean(2,:)];
figure
bar(nX,nPlotMean','stacked')
set(gca,'XTickLabel',sSignalEMGList,'XTick',[1:12],'XTickLabelRotation',90)
ylabel('R2')

% data for Tabel 2
nDataRatM = reshape(nanmedian(nDataRat,1),10,12);
nDataGainM = reshape(nanmedian(nDataGain,1),10,12);
prctile(nDataRatM,[25 50 75])
prctile(nDataGainM,[25 50 75])
%Plotting & Sidak correction: 1-(1-0.05)^(1/3)
sTitle2 = ['DiffMedianMEPmag_allTrialType_allSubject_allSignal'];
cCol = 'mbr';
[hFig2, hPlot2]  = setPlot('nRow',3,'nCol',4,'sAnnotation',sTitle2);
nPtt = []; nTtt = [];
for iSignal = 1:length(idSignalList)
    nPtt = [];
    for iTrialType = 1:numel(idTrialTypeList)
        if iTrialType==1
            nDiff1 = nDiff(:,1,iSignal);
            nDiff2 = nDiff(:,2,iSignal);
        elseif iTrialType==2
            nDiff1 = nDiff(:,1,iSignal);
            nDiff2 = nDiff(:,3,iSignal);
        elseif iTrialType==3
            nDiff1 = nDiff(:,2,iSignal);
            nDiff2 = nDiff(:,3,iSignal);
        end
        [~,P,~,stats] = ttest(nDiff1,nDiff2,'Alpha',0.017); % Sidak correction for 3 trial types
        nPtt(iTrialType) = P;
        nTtt(iTrialType) = stats.tstat;
        histogram(hPlot2(iSignal),nDiff1-nDiff2,'Orientation','horizontal',...
            'FaceColor',cCol(iTrialType))
        if iTrialType==1
            hold(hPlot2(iSignal),'on')
        end
    end
    title(hPlot2(iSignal),sSignalList{iSignal})
    display([sSignalList{iSignal},' tC',num2str(nTtt(1)),' tR',num2str(nTtt(2)),' tA',num2str(nTtt(3))])
end
legend(hPlot2(iSignal),'12','13','23')