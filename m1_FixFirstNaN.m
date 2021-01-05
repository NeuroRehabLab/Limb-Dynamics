% Did Subjects 1-6 by 13 Oct 2014
% Subjects 7 & 8 completed 25 Jun 2015
idSubject = 2:11;
idTrialType = [4,5,6,7,8,9,15,16,17];
idTrialList = getMeta('metaTrial',qry('idSubject',idSubject,'idTrialType',idTrialType,'bTrial',1));
sTable          = 'emgraw'; %table from which signals are used
sSignalList = getMeta('metaSignal',qry('sTable',sTable),'sSignal');
for idTrial = idTrialList
   for iSignal = 1:length(sSignalList)
       nData   = getSignal(idTrial,qry('sTable',sTable,'sSignal',sSignalList{iSignal}),[],'bPlot',0);
       if ~isempty(nData) && isnan(nData(1))
          iVal = find(~isnan(nData),1);
          nData(1:iVal-1) = nData(iVal);
          setSignal(nData,qry('sTable',sTable,'sSignal',sSignalList{iSignal}),idTrial);
       end
   end
end
