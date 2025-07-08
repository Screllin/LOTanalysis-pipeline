function headers = genHeaders(troughsOnly)
    headers = {'ParticipantID','Block','DepAmp','RepAmp','DepDur','RepDur',...
       'DepSlope','RepSlope','Last5Amp','Last5Dur','StimTime'};
    if (troughsOnly)
        headers = {'ParticipantID','Block','Trough', 'StimTime'};
    end
end