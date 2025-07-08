% Function to calculate the slope between two waveform points
function [slope, maximumSlope] = calculateSlope(actionPotentialSmooth, actionPotentialTime, startIndex, endIndex)
    segment = actionPotentialSmooth(startIndex:endIndex); % Voltage segment
    totalSegment = actionPotentialTime(startIndex:endIndex); % Corresponding time
    normTrace = (segment - min(segment))/(max(segment) - min(segment)); % Normalize to [0, 1]
    validIndex = find(normTrace >= 0.05 & normTrace <= 0.95); % Focus on central 90% of waveform
    if length(validIndex) > 2
        p = polyfit(totalSegment(validIndex), segment(validIndex), 1); % Linear fit
        slope = p(1); % Slope from linear fit
        maximumSlope = max(diff(segment)./diff(totalSegment)); % Instantaneous max slope
    else
        slope = NaN;
        maximumSlope = NaN;
    end
end