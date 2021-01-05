clear

idSubjectList = [2:10];
sTable = {'Rotatum'};
% idSignalKin = getMeta('metaSignal',qry('sTable',sTable));
idSignalKin = [208, 211, 214];
nRateKin = getMeta('metaSignal',qry('idSignal',idSignalKin(1)),'nRate');
sSignalList     = getMeta('metaSignal',qry('idSignal',idSignalKin),'sSignal');
idTrialTypeList = [4,5,6,7,8,9];
idSignalTorque = getMeta('metaSignal',qry('sTable','Torque','sSignal',sSignalList));
nDirFactor = [1,-1,-1];

for idSubject = idSubjectList
    disp(['idSubject: ',num2str(idSubject)]);
    tic
    idTrialList = getMeta('metaTrial',qry('idSubject',idSubject,...
        'idTrialType',          idTrialTypeList,...
        'bTrial'              ,1));
    for idTrial = idTrialList
        
        tSync   = getMeta('metaSync',qry('idTrial',idTrial),sTable);
        idSignalTrial    = getMeta('metaTrial',qry('idTrial',idTrial),'idSignal');
        idSession        = getMeta('metaTrial',qry('idTrial',idTrial),'idSession');
        idTrialType      = getMeta('metaTrial',qry('idTrial',idTrial),'idTrialType');
        
        for iSignal = 1:numel(idSignalKin)

            sSignal = sSignalList(iSignal);
            nDataRot = nDirFactor(iSignal)*getSignal(idTrial,idSignalKin(iSignal),[]);
            nPosInx = find(nDataRot > 0);
            
            nFlexion = abs(nDataRot);
            nFlexion(nPosInx) = 0;
            
            setSignal(nFlexion,qry('sTable','Flexion_Torque',...
                'sSignal',sSignal),idTrial,...
                'nRate'      ,nRateKin(1),...
                'sUnit'      ,'V',...
                'idSession'  ,idSession,...
                'idTrialType',idTrialType,...
                'tSync'      ,tSync.(sTable{1}));
            
            nNegInx = find(nDataRot < 0);
            nExtension = abs(nDataRot);
            nExtension(nNegInx) = 0;
            
            setSignal(nExtension,qry('sTable','Extension_Torque',...
                'sSignal',sSignal),idTrial,...
                'nRate'      ,nRateKin(1),...
                'sUnit'      ,'V',...
                'idSession'  ,idSession,...
                'idTrialType',idTrialType,...
                'tSync'      ,tSync.(sTable{1}));
        end
        
    end
    toc
end
