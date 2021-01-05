idTrialTypeList = getMeta('metaTrial',qry('idTrialType',[7,8,9],'idSubject',[12])); % Movement types
% % % % % % % % % nTrial = numel(idTrialTypeList); % Number of trials

sTableList = {'MarkerRaw','MarkerRaw','MarkerRaw'};
sSignalList = {'HAND3','RWP3','RELB3'}; % Three marker traces used for calculating movement onset/offset times

nRate       = getMeta('metaSignal', qry('sTable','MarkerRaw','sSignal','RWP3'),'nRate'); % Sampling rate of marker data
idSignal    = getMeta('metaSignal', qry('sTable',sTableList,'sSignal',sSignalList));

% Remove 'on' and 'off' event markers from HAND3 if they already exist
modEvent(idTrialTypeList,qry('idSignal',idSignal(1),'idEventType',1),'bEvent',0);
modEvent(idTrialTypeList,qry('idSignal',idSignal(1),'idEventType',2),'bEvent',0);

for iTrial = idTrialTypeList%([1:16 18:132 134:end]) <--Skip trials where on/off marking fails
    idTrialType = getMeta('metaTrial',qry('idTrial', iTrial),'idTrialType');
   
    nData = [];
    nData(:,1) = getSignal(iTrial,idSignal(1),[],'sScript','nDataTrial = interpolate(nDataTrial);'); % HAND3
    nData(:,2) = getSignal(iTrial,idSignal(2),[],'sScript','nDataTrial = interpolate(nDataTrial);'); % RWP3
    nData(:,3) = getSignal(iTrial,idSignal(3),[],'sScript','nDataTrial = interpolate(nDataTrial);'); % RELB3
    
    %%%%%%%%%%%%%%%%%%%%% Kinematic normalization %%%%%%%%%%%%%%%%%%%%%%%%  
    % The three extracted kinematic traces for HAND3, RWP3, and RELB3 will
    % be flipped if necessary such that all traces are positive going in
    % the z-direction.
    
    if ismember(idTrialType,[4,6]) % Movement 1, downward motion, 
        nData = -1.*nData;            % sign change implemented to make
    end                               % positive going
    
    for i = 1:3  % Normalize markers' traces to span values from 0 to 1.
        nData(:,i) = (nData(:,i) - min(nData(:,i)))./max(nData(:,i) - min(nData(:,i)));
    end
    
    nDataMean = nanmean(nData,2);          % Mean of the three traces
    
    % 2Hz Low-Pass filter
    % Naturally, this does smooth the mean trace quite a bit
    nDataMean = butterfilt(nDataMean,nRate,2,'sType','low'); % 2Hz Low-Pass filter
    
    %%%%%%%%%%%%%%%%%%%% Find/mark onsets and offsets %%%%%%%%%%%%%%%%%%%%
    % Marking onsets and offsets is accomplished by the finding maxima of 
    % the 3rd derivative (aka 'jerk') of the normalized, meaned, and 
    % low-pass filtered marker data in the z (vertical) dimension from above.
    
    nVel = diff(nDataMean(:,1))*nRate; % 1st derivative of position: velocity
    nAccel = diff(nVel(:,1))*nRate;    % 2nd derivative of positiosn: acceleration
    nJerk = diff(nAccel(:,1))*nRate;   % 3rd derivative of position: jerk
    
    % Peakfinding with a somewhat random threshold of 20. This threshold is
    % used to prevent mis-identification of jerk peaks that occur at
    % movement onset/offset.
    [peaks,locs] = findpeaks(nJerk(200:end),'minpeakheight',8); 
    locs = locs + 200;
    % Set 'on' and 'off' events in HAND3 marker trace using the
    % onset/offset jerk times.
    try
%     setEvent(locs(1)/nRate,iTrial,idSignal(1),'on');
%     setEvent(locs(2)/nRate,iTrial,idSignal(1),'off');
    catch me
        disp(['Error in Trial #: ',num2str(iTrial)]);
    end
    
end
  