%% Filter EMG
clear
% close all
idSubject = 2:11;

%%
idTrialList = getMeta('metaTrial',qry('idSubject',idSubject,'bTrial',1));
sScript = ['nDataTrial = butterfilt(nDataTrial,nRate,10,''nOrder'',2,''sType'',''high'');'];
sTable          = 'emgraw';
sSignalList     = getMeta('metaSignal',qry('sTable',sTable),'sSignal');
sSignalList(16:end) = [];
idSignalList    = getMeta('metaSignal',qry('sTable',sTable),'idSignal');
idSignalList(16:end) = [];
nRate = getMeta('metaSignal',qry('idSignal',idSignalList(1)),'nRate');


%%
tic
inx = 1;
multiWaitbar( 'Filtering EMG', 'Reset');
for idTrial = idTrialList
    multiWaitbar( 'Filtering EMG', inx/numel(idTrialList));
    tSync   = getMeta('metaSync',qry('idTrial',idTrial),sTable);
    idSignalTrial    = getMeta('metaTrial',qry('idTrial',idTrial),'idSignal');
    idSession        = getMeta('metaTrial',qry('idTrial',idTrial),'idSession');
    idTrialType      = getMeta('metaTrial',qry('idTrial',idTrial),'idTrialType');
    for iSignal = 1:length(idSignalList)
        sSignal = sSignalList{iSignal};
        if ~strcmp(sSignal(1:3),'Tim') && ~strcmp(sSignal(1:3),'Pha') && ~strcmp(sSignal(1:3),'TMS')
            nData   = getSignal(idTrial,idSignalList(iSignal),[],...
                'sScript',sScript,'bPlot',0);
            %% set signal
            setSignal(nData,qry('sTable','emg',...
                'sSignal',sSignal),idTrial,...
                'nRate'      ,nRate,...
                'sUnit'      ,'V',...
                'idSession'  ,idSession,...
                'idTrialType',idTrialType,...
                'tSync'      ,tSync);
        end
        %         else
        %             disp(['Signal ',sSignalList(iSignal),' is missing in trial ',num2str(idTrial),])
        %         end
    end
    inx=inx+1;
    save_TrialData(idTrial,'bVerbose',0)
end
toc