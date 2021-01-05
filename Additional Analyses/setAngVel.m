idSubjectList = [6:10];
idTrialTypeList = [4,5,6,7,8,9];
sTable = 'Euler';
idSignalList = getMeta('metaSignal',qry('sTable',sTable));
sSignalList = getMeta('metaSignal',qry('sTable',sTable),'sSignal');
nRate = getMeta('metaSignal',qry('sTable','Euler','idSignal','71'),'nRate');
tic
for idSubject = idSubjectList
    idTrialList = getMeta('metaTrial',qry('idSubject',idSubject,'idTrialType',idTrialTypeList,'bTrial',1));
    for idTrial = idTrialList
        idSession        = getMeta('metaTrial',qry('idTrial',idTrial),'idSession');
        idTrialType      = getMeta('metaTrial',qry('idTrial',idTrial),'idTrialType');
        tSync            = getMeta('metaSync',qry('idTrial',idTrial),sTable);
        for idSignal =  1:numel(idSignalList)
            sSignal = sSignalList(idSignal);
            nData = getSignal(idTrial,idSignalList(idSignal),[]);
            nVel = diff(nData)*nRate;
            nVelFilt = butterfilt(nVel,nRate,10);
            
            setSignal(nVelFilt,qry('sTable','AngVel','sSignal',sSignal),idTrial,...
                'nRate'      ,nRate,...
                'sUnit'      ,'V',...
                'idSession'  ,idSession,...
                'idTrialType',idTrialType,...
                'tSync'      ,tSync);
        end
    end
end
toc