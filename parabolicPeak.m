% Function for sub-sample peak timing using parabolic interpolation
function peakTime = parabolicPeak(t, y, index)
    if index <= 1 || index >= length(y)
        peakTime = t(index);% Edge case: no interpolation
        return;
    end
    x = [-1 0 1]; % Relative sample positions
    yValues = y(index-1:index+1); % 3 neighboring amplitude values
    p = polyfit(x, yValues, 2); % Fit parabola y = ax² + bx + c
    peakOffset = -p(2)/(2*p(1)); % Vertex of parabola (x = -b/2a)
    peakTime = t(index) + peakOffset*(t(2)-t(1)); % Convert to real time using sample spacing
end