idSubjectList = 3:11;

for idSubject = idSubjectList
    disp(['Running Simulations for Subject ',num2str(idSubject)]);
    
    cd('B:\BitBucket\TMS Project\2 - EMG and Kinematic Processing\Analysis Scripts');
    
   EulerTor(idSubject);
    save_MetaData('bVerbose',0)
    save_TrialData([],'bVerbose',0)
end