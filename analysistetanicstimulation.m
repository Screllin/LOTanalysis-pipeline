
% Refactored to process a single channel
function validResults =...
    analysistetanicstimulation(participantID, data, samplingInterval,...
    channelIndex, removeStimArtefact, troughsOnly)

    CONSTANT_RESULTS_ROWS = 20;
    CONSTANT_BLOCK_DURATION = 10; % Seconds- set a constant block duration to segment data 
    CONSTANT_PEAK_MOVEMENT_TOLERANCE= 0.2;    % ms tolerance for shifting of peaks/trough between blocks before visual verification is re-triggered
    % CONSTANT_PEAK_MOVEMENT_TOLERANCE= 1000; % for testing


    % Preallocate array for results, assuming max 0 blocks (over-ridden with overflow later keeps it efficient)
    resultsRows = CONSTANT_RESULTS_ROWS;
    outfileHeaders = genHeaders(troughsOnly);
    resultsCols = length(outfileHeaders);
    results = cell(resultsRows, resultsCols);

    resultIndex = 1; %initialize results index counter
    samplingRate = 1e6 / samplingInterval;  % Compute the sampling rate in Hz (1e6 converts microseconds to seconds)
    [numberSamples, numberChannels, numberSweeps] = size(data); % Extract size information from the data (samples, channels, sweeps)
    sweepsPerBlock = max(1, round(CONSTANT_BLOCK_DURATION / (numberSamples/samplingRate))); % Determine how many sweeps make up one block (based on block duration)
    numberBlock = floor(numberSweeps / sweepsPerBlock); % Calculate the number of blocks that can be formed from the data

    % Initialize reference structures to store previously detected features

    chReference =  struct('trough', NaN, 'firstPeak', NaN, 'secondPeak', NaN, 'delta', NaN);
    if (troughsOnly)
        chReference = struct('trough', NaN, 'delta', NaN);
    end

    %--------------------------------------------------------------
    % Process appropriate channel 
   
    % Loop over all blocks 
    for blockIndex = 1:numberBlock
        
        fprintf('\nProcessing %s - Block %d/%d\n Channel %d', ...
            participantID, blockIndex, numberBlock, channelIndex);  % Display progress to user

        % Determine the sweep indices for the current block% Determine the sweep indices for the current block
        startOfSweep = (blockIndex-1)*sweepsPerBlock + 1;
        endOfSweep = min(blockIndex*sweepsPerBlock, numberSweeps);

         % Average data across sweeps 
        averageData = mean(data(:, channelIndex, startOfSweep:endOfSweep), 3);

        % Generate a time vector in seconds
        timeVector = (0:numberSamples-1)/samplingRate;

        % If the includePeaks input argument is not set to true it will default to the troughs only approach for channel 2
        try
            if (~troughsOnly)
                 [parametersCh1, troughCh1, stimulationCh1, firstPeakCh1, secondPeakCh1] = ...  
                processBlockTroughsPeaks(averageData, timeVector, samplingRate, participantID, blockIndex,...
                chReference, CONSTANT_PEAK_MOVEMENT_TOLERANCE, removeStimArtefact);
            
                % Update reference features if all values are valid
                if ~isnan(troughCh1) && ~isnan(stimulationCh1) && ~isnan(firstPeakCh1) && ~isnan(secondPeakCh1)
                    chReference.trough = troughCh1;
                    chReference.firstPeak = firstPeakCh1;
                    chReference.secondPeak = secondPeakCh1;
                    chReference.delta = troughCh1 - stimulationCh1;
                end

                % Accrue the results
                if ~isempty(parametersCh1)
                    %automatically grows the array beyond 20 files if needed
                    results(resultIndex, :) = [{participantID}, {blockIndex}, parametersCh1{:}, {stimulationCh1}]; 
                    resultIndex = resultIndex + 1;
                end
            
            else
                [troughCh2, stim2] = ... 
                    processBlockTroughs(averageData, timeVector, samplingRate, participantID, blockIndex,...
                    chReference, CONSTANT_PEAK_MOVEMENT_TOLERANCE, removeStimArtefact);
                
                % Update reference features if values are valid
                    if ~isnan(troughCh2) && ~isnan(stim2) 
                        chReference.trough = troughCh2;
                        chReference.delta = troughCh2 - stim2;
                    end
                    
                % Accrue the results    
                if ~isnan(troughCh2)
                    results(resultIndex, :) = [{participantID}, {blockIndex}, {troughCh2}, {stim2}]; 
                    resultIndex = resultIndex + 1;
                end
            
            end
            
        catch ex % If any error occurs while processing a block, show a warning and exit
            warning('Block %d Channel %d processing failed: %s ', blockIndex, channelIndex, ex.message);
            return;
        end
    end % for all blocks
        
    % Initialise the return argument
    validResults = results(1:resultIndex-1, :);   % Remove any unused rows from the results matrix

end % end function

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% % This function processes a single block of averaged channel 1 data to identify 
% % the stimulus, trough, and two peaks, and then compute quantitative feature
% 
function [parameters, troughTime, stimulationTime, firstPeakTime, secondPeakTime] =...
    processBlockTroughsPeaks(trace, timeVector, samplingRate, fileName, blockNumber,...
    previousReference, peakTolerance, removeStimArtefact)
    
    parameters = [];
    troughTime = NaN; stimulationTime = NaN; firstPeakTime = NaN; secondPeakTime = NaN;  % Initialize outputs to empty or NaN in case of failure or early exit
    
    % Redoing FOC for clarity and to reduce code dup
    
    triggerManual = false;
    triggerMessage = "";
    
    try % Try to identify the stimulation time using a dedicated function
        stimulationSample = findStimulus(trace, samplingRate); 
        if isempty(stimulationSample) % If no stimulation is found, skip this block
            warning('No stimulus Ch1 Block%d', blockNumber);
            return;
        end
        stimulationTime = (stimulationSample / samplingRate) * 1000; % ms - convert sample index of stimulation to time in milliseconds

        % Determine where to begin looking for the action potential, optionally skipping stimulus artifact
        if (removeStimArtefact)
            actionPotentialStart = max(1, stimulationSample + round(0.0002*samplingRate));%where you change how much is cut off after stimulus initiation
        else
            preStimulusOffset = -round(0.001 * samplingRate); % Negative offset = earlier start  
            actionPotentialStart = max(1, stimulationSample + preStimulusOffset);     
        end

        actionPotentialEnd = min(length(trace), actionPotentialStart + round(0.02*samplingRate));  % Define the endpoint of the action potential search window, 10 ms after start
        if actionPotentialStart >= actionPotentialEnd    % If the start and end don't define a valid range, skip
            warning('Invalid window Ch1 Block%d', blockNumber);
            return;
        end

        % Extract the segment of data that likely contains the action potential
        actionPotentialSegment = trace(actionPotentialStart:actionPotentialEnd);
        actionPotentialTime = (0:length(actionPotentialSegment)-1)/samplingRate * 1000 + (actionPotentialStart/samplingRate*1000);  % Generate a time vector for this segment in m
        actionPotentialSmooth = smoothdata(actionPotentialSegment, 'gaussian', max(3, round(0.0002*samplingRate))); % Smooth the signal using a Gaussian window to reduce noise
        [~, troughLocation] = findpeaks(-actionPotentialSmooth); % Identify troughs (negative peaks) by searching for peaks in the inverted signal
        [~, peakLocation] = findpeaks(actionPotentialSmooth); % Identify positive peaks for depolarization/repolarization phases
        % --- Determine if this is the first block or if auto-tracking should be used --- %
        if blockNumber == 1 || isnan(previousReference.trough)

            triggerManual = true;
            triggerMessage = sprintf('%s Block %d (manual)', fileName, blockNumber);
        else

            % Introduce function here to improve readability
            [troughTime, firstPeakCandidates, secondPeakCandidates] = genPeakCandidates(previousReference, troughLocation, peakLocation,...
            actionPotentialTime, actionPotentialSmooth);

            % --- Determine if this is the first block or if auto-tracking should be used --- %
            if isempty(firstPeakCandidates) || isempty(secondPeakCandidates)

                triggerManual = true;
                triggerMessage = sprintf('%s Block %d (no auto peak)', fileName, blockNumber);
            else

                % Introduce function to remove code duplication

                % Pre-peak: closest to previous
                firstPeakTime = findClosestPeakTime(previousReference.firstPeak, firstPeakCandidates,...
                    actionPotentialTime, actionPotentialSmooth);

                % Post-peak: closest to previous
                secondPeakTime = findClosestPeakTime(previousReference.secondPeak, secondPeakCandidates,...
                    actionPotentialTime, actionPotentialSmooth);


                % Only trigger manual if any deviate >0.2 ms
                deviationTrough = abs(troughTime - previousReference.trough);
                deviationFirstPeak = abs(firstPeakTime - previousReference.firstPeak);
                deviationSecondPeak = abs(secondPeakTime - previousReference.secondPeak);
                if any([deviationTrough, deviationFirstPeak, deviationSecondPeak] > peakTolerance)

                        triggerManual = true;
                        triggerMessage = sprintf('%s Block %d (Δtrough/peak)', fileName, blockNumber);
                end
            end
        end

        if (triggerManual)
            [troughTime, firstPeakTime, secondPeakTime] =...
                manualTroughAndPeaks(actionPotentialSmooth, actionPotentialTime, troughLocation, peakLocation, triggerMessage);
        end


         % If all three points are successfully found, calculate quantitative parameters
        if ~isnan(troughTime) && ~isnan(firstPeakTime) && ~isnan(secondPeakTime)
            parameters = calculateParametersWithPeaks(actionPotentialSmooth, actionPotentialTime, troughTime, firstPeakTime, secondPeakTime);
        end

    catch ME
        %If an error occurs in processing this block, output a warning
        warning('Ch1 Block%d error: %s', blockNumber, ME.message);
    end
end


function [troughTime, stimulationTime] = processBlockTroughs(trace, timeVector, samplingRate, fileName, blockNumber, previousReference, troughTolerance, removeStimArtefact)
    troughTime = NaN; stimulationTime = NaN;
    try
        stimulationSample = findStimulus(trace, samplingRate); % Attempt to locate the stimulus event in the signal
        if isempty(stimulationSample)  % If the stimulus is not found, issue warning and skip processing
            warning('No stimulus Ch2 Block%d', blockNumber);
            return;
        end
        stimulationTime = (stimulationSample / samplingRate) * 1000; % ms % Convert stimulus location from samples to milliseconds
                % Start the analysis window slightly before the stimulus to include the artifact

        % Determine where to begin looking for the action potential, optionally skipping stimulus artifact
        if (removeStimArtefact)
            actionPotentialStart = max(1, stimulationSample + round(0.0002*samplingRate));%where you change how much is cut off after stimulus initiation
        else
            preStimulusOffset = -round(0.001 * samplingRate); % Negative offset = earlier start  
            actionPotentialStart = max(1, stimulationSample + preStimulusOffset);     
        end

        actionPotentialEnd = min(length(trace), actionPotentialStart + round(0.01*samplingRate));
        % Sanity check for valid window
        if actionPotentialStart >= actionPotentialEnd
            warning('Invalid window Ch2 Block%d', blockNumber);
            return;
        end
         % Extract the trace segment for analysis
        actionPotentialSegment = trace(actionPotentialStart:actionPotentialEnd);
        % Generate corresponding time vector (in ms)
        actionPotentialTime = (0:length(actionPotentialSegment)-1)/samplingRate * 1000 + (actionPotentialStart/samplingRate*1000); % Generate corresponding time vector (in ms)
        actionPotentialSmooth = smoothdata(actionPotentialSegment, 'gaussian', max(3, round(0.0002*samplingRate))); % Smooth the data to reduce noise and facilitate peak detection
        [~, troughLocation] = findpeaks(-actionPotentialSmooth); % Detect troughs (most negative deflections)

        if blockNumber == 1 || isnan(previousReference.trough)
            troughTime = manualTroughSelection(actionPotentialSmooth, actionPotentialTime, troughLocation, ...
                sprintf('%s Ch2 Block%d (manual)', fileName, blockNumber));
        else % Find the detected trough closest to the previous block's trough
            troughDifference = abs(actionPotentialTime(troughLocation) - previousReference.trough);
            [minimumDifference, troughIndexRel] = min(troughDifference);
            troughIndex = troughLocation(troughIndexRel);
            troughTime = parabolicPeak(actionPotentialTime, actionPotentialSmooth, troughIndex);
            % If the new trough location is too far from the reference, redo manually
            if minimumDifference > troughTolerance
                troughTime = manualTroughSelection(actionPotentialSmooth, actionPotentialTime, troughLocation, ...
                    sprintf('%s Ch2 Block%d (Δtrough)', fileName, blockNumber));
            end
        end
    catch ME % Report any error encountered during block processing
        warning('Ch2 Block%d error: %s', blockNumber, ME.message);
    end
end

function [troughTime, firstPeakTime, secondPeakTime] = manualTroughAndPeaks(actionPotentialSmooth, actionPotentialTime, troughLocation, peakLocation, TitleString)
    troughTime = NaN; firstPeakTime = NaN; secondPeakTime = NaN; % Initialize all outputs as NaN in case the user cancels the selection
    accepted = false; % Loop until the user accepts the selected points
    while ~accepted
         % Create a new figure window for manual selection
        fig = figure('Name', TitleString, 'NumberTitle','off');
        clf(fig); % Clear figure window
        plot(actionPotentialTime, actionPotentialSmooth, 'k-'); hold on;
        plot(actionPotentialTime(troughLocation), actionPotentialSmooth(troughLocation), 'ro', 'MarkerSize',8); % Mark detected troughs with red circles
        plot(actionPotentialTime(peakLocation), actionPotentialSmooth(peakLocation), 'go', 'MarkerSize',8); % Mark detected peaks with green circles
        legend('Trace','Troughs','Peaks');

         % Prompt user to click the trough
        title({'Select the TROUGH (red circles).'});
        [x1,~] = ginput(1);
        [~, troughIndex] = min(abs(actionPotentialTime(troughLocation) - x1));
        troughX = actionPotentialTime(troughLocation(troughIndex));
        troughY = actionPotentialSmooth(troughLocation(troughIndex));
        troughTime = parabolicPeak(actionPotentialTime, actionPotentialSmooth, troughLocation(troughIndex));
        plot(troughX, troughY, 'ko', 'MarkerSize',12, 'MarkerFaceColor','k');

        % Prompt user to click the pre-trough peak
        title({'Select the FIRST-PEAK (pink circles BEFORE trough)'});
        firstPeaks = peakLocation(peakLocation < troughLocation(troughIndex));
        if isempty(firstPeaks)
            firstPeakTime = NaN;
        else
            plot(actionPotentialTime(firstPeaks), actionPotentialSmooth(firstPeaks), 'mo', 'MarkerSize',8);
            [x2,~] = ginput(1);
            [~, firstIndex] = min(abs(actionPotentialTime(firstPeaks) - x2));
            firstPeakX = actionPotentialTime(firstPeaks(firstIndex));
            firstPeaky = actionPotentialSmooth(firstPeaks(firstIndex));
            firstPeakTime = parabolicPeak(actionPotentialTime, actionPotentialSmooth, firstPeaks(firstIndex));
            plot(firstPeakX, firstPeaky, 'ko', 'MarkerSize',12, 'MarkerFaceColor','k');
        end

        % Prompt user to click the post-trough peak
        title({'Select the SECOND-PEAK (blue circles AFTER trough)'});
        secondPeaks = peakLocation(peakLocation > troughLocation(troughIndex));
        if isempty(secondPeaks)
            secondPeakTime = NaN;
        else
            plot(actionPotentialTime(secondPeaks), actionPotentialSmooth(secondPeaks), 'co', 'MarkerSize',8);
            [x3,~] = ginput(1);
            [~, secondIndex] = min(abs(actionPotentialTime(secondPeaks) - x3));
            secondPeakX = actionPotentialTime(secondPeaks(secondIndex));
            secondPeaky = actionPotentialSmooth(secondPeaks(secondIndex));
            secondPeakTime = parabolicPeak(actionPotentialTime, actionPotentialSmooth, secondPeaks(secondIndex));
            plot(secondPeakX, secondPeaky, 'ko', 'MarkerSize',12, 'MarkerFaceColor','k');
        end
% Final instructions for the user
        title({'Selected points: black filled circles.', ...
               'SPACE=accept, n=redo/back'});
        hold off;
        % Wait for user input
        waitforbuttonpress;
        key = get(fig, 'CurrentKey');
        if strcmpi(key, 'space')
            accepted = true;
        elseif strcmpi(key, 'n')
            close(fig);
            continue; % redo
        end
        close(fig);
    end
end
% Function for manually selecting a trough from an electrophysiological trace
function troughTime = manualTroughSelection(actionPotentialSmooth, actionPotentialTime, troughLocation, TitleString)
    troughTime = NaN; % Initialize output as NaN in case selection fails
    accepted = false; % Control flag to repeat selection until accepted
    while ~accepted % Loop until user accepts the selected point
        fig = figure('Name', TitleString, 'NumberTitle','off'); % Create figure window with a custom title
        clf(fig); % Clear figure content
        plot(actionPotentialTime, actionPotentialSmooth, 'k-'); hold on; % Plot the signal in black
        plot(actionPotentialTime(troughLocation), actionPotentialSmooth(troughLocation), 'ro', 'MarkerSize',8); % Overlay red circles at detected troughs
        title({'Select the TROUGH (red circles).'});  % Prompt user
        [x,~] = ginput(1); % Get user input (x = time) from a mouse click
        [~, index] = min(abs(actionPotentialTime(troughLocation) - x)); % Find the closest detected trough to the click
        selectLocation = troughLocation(index);  % Get index in the trace
        troughTime = parabolicPeak(actionPotentialTime, actionPotentialSmooth, selectLocation); % Interpolate for sub-sample timing
        plot(actionPotentialTime(selectLocation), actionPotentialSmooth(selectLocation), 'ko', 'MarkerSize',12, 'MarkerFaceColor','k');  % Show selected point
        % Prompt instructions for confirmation or redo
        title({'Selected trough: black filled circle.', ...
               'SPACE=accept, n=redo/back'});
        hold off; % Finalize plot
        waitforbuttonpress; % Wait for keyboard input
        key = get(fig, 'CurrentKey'); % Capture the pressed key
        if strcmpi(key, 'space')
            accepted = true; % Accept selection
        elseif strcmpi(key, 'n') % Close figure and retry
            close(fig); % Close figure after decision
            continue; % redo
        end
        close(fig);
    end
end
% Function to detect stimulation time based on signal derivative
function stimulationSample = findStimulus(trace, samplingRate)
    dTrace = diff(trace); % First derivative of the trace
     % Threshold is 10× standard deviation of the first 0.5s of derivative signal
    thresh = 10 * std(dTrace(1:min(round(0.5*samplingRate), end)));
    stimulationSamples = find(abs(dTrace) > thresh); % Find large jumps in derivative
    minimumInterval = round(0.02 * samplingRate); % Require 20 ms spacing to avoid counting same event multiple times
    validStimulus = [true; diff(stimulationSamples) > minimumInterval];  % Mask to remove closely spaced duplicates
    stimulationSample = stimulationSamples(validStimulus); % Apply mask
    if ~isempty(stimulationSample)
        stimulationSample = stimulationSample(1); % Return first detected stimulus
    else
        stimulationSample = []; % Return empty if none found
    end
end
% Function to calculate waveform metrics based on identified peaks
function parameters = calculateParametersWithPeaks(actionPotentialSmooth, actionPotentialTime, troughTime, firstPeakTime, secondPeakTime)
    % Convert times to nearest indices in the data arrays
    troughIndex = find(abs(actionPotentialTime - troughTime) == min(abs(actionPotentialTime - troughTime)), 1);
    firstIndex = find(abs(actionPotentialTime - firstPeakTime) == min(abs(actionPotentialTime - firstPeakTime)), 1);
    secondIndex = find(abs(actionPotentialTime - secondPeakTime) == min(abs(actionPotentialTime - secondPeakTime)), 1);
% Depolarization: amplitude and duration
    depolarisationAmplitude = actionPotentialSmooth(firstIndex) - actionPotentialSmooth(troughIndex);
    depolarisationDuration = troughTime - firstPeakTime;
    [depolarisationSlope, ~] = calculateSlope(actionPotentialSmooth, actionPotentialTime, firstIndex, troughIndex);
% Repolarization: amplitude and duration
    repolarisationAmplitude = actionPotentialSmooth(secondIndex) - actionPotentialSmooth(troughIndex);
    repolarisationDuration = secondPeakTime - troughTime;
    [repolarisationSlope, ~] = calculateSlope(actionPotentialSmooth, actionPotentialTime, troughIndex, secondIndex);
% Final phase: last 5% of repolarization
    [last5Amplitude, last5Duration] = calculatelast5percent(actionPotentialSmooth, actionPotentialTime, troughIndex, secondIndex);
   % Combine all calculated parameters into a cell array
    parameters = {depolarisationAmplitude, repolarisationAmplitude, depolarisationDuration, repolarisationDuration, abs(depolarisationSlope), repolarisationSlope, last5Amplitude, last5Duration};
end

% Function to estimate the amplitude and duration of the final 5% of repolarization
function [last5Amplitude, last5Duration] = calculatelast5percent(actionPotentialSmooth, actionPotentialTime, troughIndex, secondIndex)
    % Only consider repolarisation segment
    segment = actionPotentialSmooth(troughIndex:secondIndex);
    totalSegment = actionPotentialTime(troughIndex:secondIndex);

    % Calculate repolarisation amplitude (second peak - trough)
    repolAmp = actionPotentialSmooth(secondIndex) - actionPotentialSmooth(troughIndex);

    % 95% repolarisation threshold value
    threshold95 = actionPotentialSmooth(troughIndex) + 0.95 * repolAmp;

    % Find the first index where the segment crosses the 95% threshold
    crossIndex = find(segment >= threshold95, 1, 'first');

    if isempty(crossIndex) || crossIndex == 1
        last5Duration = NaN;
        last5Amplitude = NaN;
        return;
    end

    % Parabolic interpolation for sub-sample timing
    if crossIndex < numel(segment)
        t = totalSegment(crossIndex-1:crossIndex+1);
        y = segment(crossIndex-1:crossIndex+1);
        exactTime = interpolateCrossing(t, y, threshold95);
    else
        exactTime = totalSegment(crossIndex);
    end

    % Duration is from 95% crossing to second peak (end of repolarisation)
    last5Duration = totalSegment(end) - exactTime;

    % Last 5% amplitude is always 5% of repolarisation amplitude
    last5Amplitude = repolAmp * 0.05;
end

function exactTime = interpolateCrossing(t, y, target)
    % Parabolic fit and solve for y = target
    A = [t(:).^2, t(:), ones(3,1)];
    coeffs = A\y(:);
    a = coeffs(1); b = coeffs(2); c = coeffs(3);
    discriminant = b^2 - 4*a*(c - target);
    if discriminant >= 0 && a ~= 0
        t1 = (-b + sqrt(discriminant)) / (2*a);
        t2 = (-b - sqrt(discriminant)) / (2*a);
        validRoots = [t1, t2];
        validRoots = validRoots(validRoots >= min(t) & validRoots <= max(t));
        if ~isempty(validRoots)
            exactTime = validRoots(1);
        else
            exactTime = t(2); % fallback
        end
    else
        exactTime = t(2); % fallback
    end
end

function [troughTime, firstPeakCandidates, secondPeakCandidates] =...
    genPeakCandidates(previousReference, troughLocation, peakLocation,...
    actionPotentialTime, actionPotentialSmooth)

 % Find trough closest to previous
    troughDifference = abs(actionPotentialTime(troughLocation) - previousReference.trough);
    [~, troughIndexRel] = min(troughDifference);
    troughIndex = troughLocation(troughIndexRel);
    troughTime = parabolicPeak(actionPotentialTime, actionPotentialSmooth, troughIndex);

    % Peaks: highest within 4ms before/after, closest to previous
    firstPeakCandidates = peakLocation(actionPotentialTime(peakLocation) < actionPotentialTime(troughIndex) & actionPotentialTime(peakLocation) >= actionPotentialTime(troughIndex) - 4); % duration before trough to find first peak
    secondPeakCandidates = peakLocation(actionPotentialTime(peakLocation) > actionPotentialTime(troughIndex) & actionPotentialTime(peakLocation) <= actionPotentialTime(troughIndex) + 6); % duration to find second peak make longer if needed
end


function peakTime = findClosestPeakTime(referencePeak, peakCandidates, actionPotentialTime, actionPotentialSmooth)
    peakDifference = abs(actionPotentialTime(peakCandidates) - referencePeak);
    [~, peakIndexRel] = min(peakDifference);
    peakIndex = peakCandidates(peakIndexRel);
    peakTime = parabolicPeak(actionPotentialTime, actionPotentialSmooth, peakIndex);
end
