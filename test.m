
filename = "22d13002.abf";

[data, samplingInterval] = abfload(filename); %load the data and sampling interval from ABF file

[~,participantID,~] = fileparts(filename);

channelIndex = 1;
removeStimArtefact = true;
troughsOnly = false;

tetanicResults = analysistetanicstimulation(participantID, data, samplingInterval, channelIndex, removeStimArtefact, troughsOnly);
genOutfile(participantID, tetanicResults, channelIndex, troughsOnly);

channelIndex = 2;
removeStimArtefact = false;
troughsOnly = true;
tetanicResults = analysistetanicstimulation(participantID, data, samplingInterval, channelIndex, removeStimArtefact, troughsOnly);
genOutfile(participantID, tetanicResults, channelIndex, troughsOnly);
