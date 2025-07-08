function processFolder(channelIndex, removeStimArtefact, troughsOnly)
 
    % Get folder containing input files
    folderPath = uigetdir(pwd, 'Select your ABF files folder');
    if folderPath == 0, return; end %if file selection cancelled, exit function

    %get a list of all ABF files in slected folder
    fileList = dir(fullfile(folderPath, '*.abf')); 

    % Loop over all files
    for fileIndex = 1:length(fileList) %loop through each ABF file in the folder

        filename = fileList(fileIndex).name;
        fullPath = fullfile(folderPath, filename);

        % Exception handling here specifically for file load error

        %load the data and sampling interval from ABF file
        try
            [data, samplingInterval] = abfload(fullPath); 
        catch ex
            warning('Failed to load data file %s: %s ', fullPath, ex.message);
            return;
        end

        [~,participantID,~] = fileparts(filename);
        tetanicResults = analysistetanicstimulation(participantID, data, samplingInterval, channelIndex,...
            removeStimArtefact, troughsOnly);
        genOutfile(participantID, tetanicResults, channelIndex, troughsOnly);
    end % for each file
end % function processFolder


  