%% Did Subjects 1-6 by 13 Oct 2014
% FIX: nan removal in joint angles where not necessary
idTrialType = [4,5,6,7,8,9];
idSubject = 12;
idTrialList = getMeta('metaTrial',qry('idTrialType',idTrialType,'idSubject',idSubject,'bTrial',1));
sTable          = 'MarkerRaw'; %table from which signals are used
nRate           = getMeta('metaSignal',qry('sTable',sTable,...
    'sSignal','C71'),'nRate'); %sampling frequency
sScript         = ('nDataTrial = interpolate(nDataTrial);');

sCurrentPath =pwd;
sPath = [sCurrentPath,'/Functions/'];

addpath(sPath);

%%
for idTrial = idTrialList
    clear nData* nPosRaw Hkin
    
    tSync          = getMeta('metaSync',qry('idTrial',idTrial),sTable);
    
    %% Pull Kinematic Data
    nDataWA(1,:)   = getSignal(idTrial,qry('sTable',sTable,'sSignal','RWT1'),[],'bPlot',0,'sScript',sScript);% Wrist radius (thumb)
    nDataWA(2,:)   = getSignal(idTrial,qry('sTable',sTable,'sSignal','RWT2'),[],'bPlot',0,'sScript',sScript);
    nDataWA(3,:)   = getSignal(idTrial,qry('sTable',sTable,'sSignal','RWT3'),[],'bPlot',0,'sScript',sScript);
    
    nDataWB(1,:)   = getSignal(idTrial,qry('sTable',sTable,'sSignal','RWP1'),[],'bPlot',0,'sScript',sScript);% Wrist ulna (pinky)
    nDataWB(2,:)   = getSignal(idTrial,qry('sTable',sTable,'sSignal','RWP2'),[],'bPlot',0,'sScript',sScript);
    nDataWB(3,:)   = getSignal(idTrial,qry('sTable',sTable,'sSignal','RWP3'),[],'bPlot',0,'sScript',sScript);
    
    nDataF(1,:)   = getSignal(idTrial,qry('sTable',sTable,'sSignal','HAND1'),[],'bPlot',0,'sScript',sScript);% MCP marker
    nDataF(2,:)   = getSignal(idTrial,qry('sTable',sTable,'sSignal','HAND2'),[],'bPlot',0,'sScript',sScript);
    nDataF(3,:)   = getSignal(idTrial,qry('sTable',sTable,'sSignal','HAND3'),[],'bPlot',0,'sScript',sScript);
    
    nDataE(1,:)   = getSignal(idTrial,qry('sTable',sTable,'sSignal','RELB1'),[],'bPlot',0,'sScript',sScript);% Elbow marker
    nDataE(2,:)   = getSignal(idTrial,qry('sTable',sTable,'sSignal','RELB2'),[],'bPlot',0,'sScript',sScript);
    nDataE(3,:)   = getSignal(idTrial,qry('sTable',sTable,'sSignal','RELB3'),[],'bPlot',0,'sScript',sScript);
    
    nDataSR(1,:)   = getSignal(idTrial,qry('sTable',sTable,'sSignal','RSHOR1'),[],'bPlot',0,'sScript',sScript);% Rear shoulder marker
    nDataSR(2,:)   = getSignal(idTrial,qry('sTable',sTable,'sSignal','RSHOR2'),[],'bPlot',0,'sScript',sScript);
    nDataSR(3,:)   = getSignal(idTrial,qry('sTable',sTable,'sSignal','RSHOR3'),[],'bPlot',0,'sScript',sScript);
    
    nDataSF(1,:)   = getSignal(idTrial,qry('sTable',sTable,'sSignal','RSHOF1'),[],'bPlot',0,'sScript',sScript);% Front shoulder marker
    nDataSF(2,:)   = getSignal(idTrial,qry('sTable',sTable,'sSignal','RSHOF2'),[],'bPlot',0,'sScript',sScript);
    nDataSF(3,:)   = getSignal(idTrial,qry('sTable',sTable,'sSignal','RSHOF3'),[],'bPlot',0,'sScript',sScript);
    
    
    %% Final coordinates to input into model
    
    % Averaged shoulder coordinates (Front and Rear)
    nPosRaw(1,:) = (nDataSR(1,:)+nDataSF(1,:))/2;
    nPosRaw(2,:) = (nDataSR(2,:)+nDataSF(2,:))/2;
    nPosRaw(3,:) = (nDataSR(3,:)+nDataSF(3,:))/2;
    
    % Elbow coordinates
    nPosRaw(4,:) = nDataE(1,:);
    nPosRaw(5,:) = nDataE(2,:);
    nPosRaw(6,:) = nDataE(3,:);
    
    % Averaged wrist coordinates (Thumb and Pinky sides)
    nPosRaw(7,:) = (nDataWB(1,:)+nDataWA(1,:))/2;
    nPosRaw(8,:) = (nDataWB(2,:)+nDataWA(2,:))/2;
    nPosRaw(9,:) = (nDataWB(3,:)+nDataWA(3,:))/2;
    
    % Hand coordinates
    nPosRaw(10,:) = nDataF(1,:);
    nPosRaw(11,:) = nDataF(2,:);
    nPosRaw(12,:) = nDataF(3,:);
    
    % Convert nPosRaw from millimeters to meters
    nPosRaw = nPosRaw/1000;
    
    
    nDataW = nDataWA;  %Wrist marker coordinates
    nDataW = [nDataW;nDataWB];
    
    [nSEul,nEEul,nWEul] = getEuler_ARM_140205(nPosRaw',nRate,nDataW','bPlot',0);
    
    %% Set Euler angles
    setSignal(nSEul(1,:)',qry('sTable','Euler','sSignal','Shoulder_Right1'),idTrial,...
        'tSync',tSync,'nRate',nRate,'sUnit','rad');
    setSignal(nSEul(2,:)',qry('sTable','Euler','sSignal','Shoulder_Right2'),idTrial,...
        'tSync',tSync,'nRate',nRate,'sUnit','rad');
    setSignal(nSEul(3,:)',qry('sTable','Euler','sSignal','Shoulder_Right3'),idTrial,...
        'tSync',tSync,'nRate',nRate,'sUnit','rad');
    setSignal(nEEul(1,:)',qry('sTable','Euler','sSignal','Elbow_Right1'),idTrial,...
        'tSync',tSync,'nRate',nRate,'sUnit','rad');
    setSignal(nWEul(1,:)',qry('sTable','Euler','sSignal','Wrist_Right1'),idTrial,...
        'tSync',tSync,'nRate',nRate,'sUnit','rad');
    setSignal(nWEul(2,:)',qry('sTable','Euler','sSignal','Wrist_Right2'),idTrial,...
        'tSync',tSync,'nRate',nRate,'sUnit','rad');
    setSignal(nWEul(3,:)',qry('sTable','Euler','sSignal','Wrist_Right3'),idTrial,...
        'tSync',tSync,'nRate',nRate,'sUnit','rad');
    %% Calculate dynamics
    nEuler = [nSEul',nEEul(1,:)',nWEul'];
    % Elbow coordinates
    Hkin(1,:) = nDataE(1,:)/1000- nPosRaw(1,:);
    Hkin(2,:) = nDataE(2,:)/1000- nPosRaw(2,:);
    Hkin(3,:) = nDataE(3,:)/1000- nPosRaw(3,:);
    
    % Averaged wrist coordinates (Thumb and Pinky sides)
    Hkin(4,:) =  nPosRaw(7,:)- nPosRaw(1,:);
    Hkin(5,:) =  nPosRaw(8,:)- nPosRaw(2,:);
    Hkin(6,:) =  nPosRaw(9,:)- nPosRaw(3,:);
    
    % Hand coordinates
    Hkin(7,:) = nDataF(1,:)/1000- nPosRaw(1,:);
    Hkin(8,:) = nDataF(2,:)/1000- nPosRaw(2,:);
    Hkin(9,:) = nDataF(3,:)/1000- nPosRaw(3,:);
    
    % segment lengths
    nLength(1) = nanmean(sqrt(Hkin(1,:).^2+Hkin(2,:).^2+Hkin(3,:).^2)); %humerus length
    nLength(2) = nanmean(sqrt((Hkin(4,:)-Hkin(1,:)).^2 + (Hkin(5,:)-Hkin(2,:)).^2 + (Hkin(6,:)-Hkin(3,:)).^2));%radius/ulna length
    nLength(3) = nanmean(sqrt((Hkin(4,:)-Hkin(7,:)).^2 + (Hkin(5,:)-Hkin(8,:)).^2 + (Hkin(6,:)-Hkin(9,:)).^2));%hand length
    
    nWeight = getMeta()
    [nOut] = MODEL_getTorques(double(nEuler),double(nPosRaw'),480,double(nLength),50, 'bPlot',0);
    
    %% Set Torques in Database %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Shoulder-Right
    setSignal(nOut{1},qry('sTable','Torque','sSignal','ShoR_Mx'),idTrial,...
        'tSync',tSync,'nRate',nRate,'sUnit','Nm/rad');
    setSignal(nOut{8},qry('sTable','Torque','sSignal','ShoR_Gx'),idTrial,...
        'tSync',tSync,'nRate',nRate,'sUnit','Nm/rad');
    setSignal(nOut{29},qry('sTable','Torque','sSignal','ShoR_NetIntx'),idTrial,...
        'tSync',tSync,'nRate',nRate,'sUnit','Nm/rad');
    
    % Elbow-Right
    setSignal(nOut{4},qry('sTable','Torque','sSignal','ElbowR_Mx'),idTrial,...
        'tSync',tSync,'nRate',nRate,'sUnit','Nm/rad');
    setSignal(nOut{11},qry('sTable','Torque','sSignal','ElbowR_Gx'),idTrial,...
        'tSync',tSync,'nRate',nRate,'sUnit','Nm/rad');
    setSignal(nOut{32},qry('sTable','Torque','sSignal','ElbowR_NetIntx'),idTrial,...
        'tSync',tSync,'nRate',nRate,'sUnit','Nm/rad');
    
    % Wrist-Right
    setSignal(nOut{5},qry('sTable','Torque','sSignal','WristR_Mx'),idTrial,...
        'tSync',tSync,'nRate',nRate,'sUnit','Nm/rad');
    setSignal(nOut{12},qry('sTable','Torque','sSignal','WristR_Gx'),idTrial,...
        'tSync',tSync,'nRate',nRate,'sUnit','Nm/rad');
    setSignal(nOut{33},qry('sTable','Torque','sSignal','WristR_NetIntx'),idTrial,...
        'tSync',tSync,'nRate',nRate,'sUnit','Nm/rad');
    
    %% Set Error in Simulation Per Trial
    
    if isempty(getMeta('metaDynSim'))
        setMeta('metaDynSim','new',...
            'idSim',0,...
            'idTrial',0,...
            'nRMS',0);
    end
    
    setMeta('metaDynSim',qry('idTrial',idTrial),...
        'idTrial',idTrial,...
        'nRMS',nOut{36});
    
end
