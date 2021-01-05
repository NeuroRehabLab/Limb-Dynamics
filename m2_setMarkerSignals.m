%%% Set Marker Signal:  Take marker raw data and transform it to marker
%%% data.  This Script will filter and interpolate raw data
clear
idSubject       = 12;  %% 3,5,9,11,12,13
idSession       = idSubject;
nLim            = 0.2;
idTrialTypeList = [4:9];
sTableOld       = {'MarkerRaw'};
sTableNew       = {'Marker'};
sSignalList     = getMeta('metaSignal',qry('sTable',sTableOld),'sSignal');
for iSignal = 1:length(sSignalList)
    sName   = sSignalList{iSignal};
    bSignalList(iSignal)  = strcmp(sName(end),'A')||strcmp(sName(end),'L') ||(strcmp(sName(2),'Y') && strcmp(sName(3),'N'));
end
sSignalList(bSignalList) = [];
if idSubject==2
    sSignalList(end-2:end) = [];
end
idSignalList    = getMeta('metaSignal',qry('sTable',sTableOld,'sSignal',sSignalList));
nRate           = getMeta('metaSignal',qry('sTable',sTableOld,'sSignal',sSignalList{1}),'nRate');

for iTrialType= 1:numel(idTrialTypeList);
    idTrialType = idTrialTypeList(iTrialType);
    idTrialList = getMeta('metaTrial',qry('idSubject',idSubject,'idTrialType',idTrialType));
    for iTrial= 1:numel(idTrialList)
        tSync   = getMeta('metaSync',qry('idTrial',idTrialList(iTrial)),sTableOld);
        for iSignal= 1:numel(sSignalList)
            idSignal    = idSignalList(iSignal);
            nData = getSignal(idTrialList(iTrial),idSignal,[]);
            if ~isempty(nData) && sum(isnan(nData)) ~=length(nData);
                iNaNlast = find(isnan(nData),1,'last');
                iVallast= find(~isnan(nData),1,'last');
                iVal=nData(iVallast);
                if numel(iNaNlast)~=0 && iNaNlast==numel(nData)
                    iVallast = find(~isnan(nData),1,'last');
                    nData(iVallast+1:end) = iVal;
                end
                
                nDataF = butterfilt(nData,nRate,10,'nOrder',2);
                nDataIF = interpolate(nDataF,'nLim',nRate*nLim,'sMethod','spline');
                nDataIF = interpolate(nDataF,'nLim',nRate*nLim);

                setSignal(nDataIF,qry('sTable',sTableNew,'sSignal',sSignalList(iSignal)),idTrialList(iTrial),...
                    'nRate',nRate,...
                    'sUnit','mm',...
                    'idSession',idSession,...
                    'idTrialType',idTrialTypeList(iTrialType),...
                    'tSync',0);
            end
        end
    end
end

disp('Complete');





      

      
