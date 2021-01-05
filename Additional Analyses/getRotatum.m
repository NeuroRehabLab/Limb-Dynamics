function getRotatum

% User Defined Parameters
idSubjectList = [11:12];
idTrialTypeList = [4,5,6,7,8,9];

% sTable Information
sTable = {'Torque'};
idSignalDyn = getMeta('metaSignal',qry('sTable',sTable));
nRateDyn = getMeta('metaSignal',qry('idSignal',idSignalDyn(1)),'nRate');
nRateDyn = nRateDyn(1);
sSignalList     = getMeta('metaSignal',qry('sTable',sTable),'sSignal');

for idSubject = idSubjectList
    
    % Show Progress to User
    disp(['idSubject: ',num2str(idSubject)]);
    tic
    
    % Pull Trial List per Subject
    idTrialList = getMeta('metaTrial',qry('idSubject',idSubject,...
        'idTrialType',          idTrialTypeList,...
        'bTrial'                ,1));
    
    for idTrial = idTrialList
        
        % Pull metaInformation for Signal Writing
        tSync   = getMeta('metaSync',qry('idTrial',idTrial),sTable);
        idSignalTrial    = getMeta('metaTrial',qry('idTrial',idTrial),'idSignal');
        idSession        = getMeta('metaTrial',qry('idTrial',idTrial),'idSession');
        idTrialType      = getMeta('metaTrial',qry('idTrial',idTrial),'idTrialType');
        
        for iSignal = 1:numel(idSignalDyn)
            
            % Signal Name
            sSignal = sSignalList(iSignal);
            
            % Pull Dynamic Signals
            nDataDyn = getSignal(idTrial,idSignalDyn(iSignal),[]);
            
            % Derivative of Torque Signals
            nRotatum = diff(nDataDyn)*nRateDyn;
            
            % Write Signal to BoxSci
            setSignal(nRotatum,qry('sTable','Rotatum',...
                'sSignal',sSignal),idTrial,...
                'nRate'      ,nRateDyn,...
                'sUnit'      ,'V',...
                'idSession'  ,idSession,...
                'idTrialType',idTrialType,...
                'tSync'      ,tSync.(sTable{1}));
        end 
    end
    toc
end
