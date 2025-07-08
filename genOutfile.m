function genOutfile(participantID, tetanicResults, channelIndex, troughsOnly)
    
    headers = genHeaders(troughsOnly);
    outfilename = strcat('ephys_results_tetanic_', participantID, '_ch_', int2str(channelIndex), '.xlsx');
    writecell([headers; tetanicResults], outfilename);

end