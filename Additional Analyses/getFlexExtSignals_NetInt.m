
clear

idSubjectList = [2:10];
sTable = {'Torque'};
idSignal = [80,83,86];
nRate = getMeta('metaSignal',qry('idSignal',idSignal(1)),'nRate');
sSignalList     = getMeta('metaSignal',qry('idSignal',idSignal),'sSignal');
idTrialTypeList = [4,5,6,7,8,9];

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
        
        for iSignal = 1:numel(idSignal)

            sSignal = sSignalList(iSignal);
            
            nData = getSignal(idTrial,idSignal(iSignal),[]);
            
            nPosInx = find(nData > 0);            
            nFlexion = abs(nData);
            nFlexion(nPosInx) = 0;
            
            setSignal(nFlexion,qry('sTable','Flexion_Torque',...
                'sSignal',sSignal),idTrial,...
                'nRate'      ,nRate,...
                'sUnit'      ,'V',...
                'idSession'  ,idSession,...
                'idTrialType',idTrialType,...
                'tSync'      ,tSync.(sTable{1}));
            
            nNegInx = find(nData < 0);
            nExtension = abs(nData);
            nExtension(nNegInx) = 0;
            
            setSignal(nExtension,qry('sTable','Extension_Torque',...
                'sSignal',sSignal),idTrial,...
                'nRate'      ,nRate,...
                'sUnit'      ,'V',...
                'idSession'  ,idSession,...
                'idTrialType',idTrialType,...
                'tSync'      ,tSync.(sTable{1}));
        end
        
    end
    toc
end
