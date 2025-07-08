% Set up processing parameters 
channelIndex = 1;
removeStimArtefact = true;
troughsOnly = false; 

% Run it
processFolder(channelIndex, removeStimArtefact, troughsOnly);

% Change processing parameters 
channelIndex = 2;
removeStimArtefact = false;
troughsOnly = true; 

processFolder(channelIndex, removeStimArtefact, troughsOnly);

