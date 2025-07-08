% --- Select folder with *_ch_1.xlsx files for channel 1 ---
folderPath = uigetdir(pwd, 'Select folder with *_ch_1.xlsx files');
if isequal(folderPath, 0)
    error('No folder selected.');
end

% --- Find all *_ch_1.xlsx files ---
files = dir(fullfile(folderPath, '*_ch_1.xlsx')); %change file ending if different
if isempty(files)
    error('No matching files found.');
end

% --- Initialize storage ---
fileTables = cell(1, numel(files));
fileNames = cell(1, numel(files));
varSet = [];

% --- Read all files, store cleaned names and valid variables ---
for i = 1:numel(files)
    fullPath = fullfile(folderPath, files(i).name);
    T = readtable(fullPath);

    % Clean name
    [~, base, ~] = fileparts(files(i).name);
    cleanName = erase(base, '_ch_1');
    fileNames{i} = cleanName;

    % Store table
    fileTables{i} = T;

    % Get valid (numeric) variable names only
    numericVars = T.Properties.VariableNames(varfun(@isnumeric, T, 'OutputFormat', 'uniform'));
    varSet = union(varSet, numericVars);  % Unique union of numeric variables
end

% --- Output Excel path ---
outputFile = fullfile(folderPath, 'Compiled_Results.xlsx'); %change output file name if needed

% --- Loop through each variable ---
for v = 1:length(varSet)
    varName = varSet{v};
    maxLen = 0;

    % Find max column length
    for i = 1:length(fileTables)
        T = fileTables{i};
        if ismember(varName, T.Properties.VariableNames)
            maxLen = max(maxLen, height(T));
        end
    end

    % Create output sheet cell array
    outputSheet = cell(maxLen + 1, numel(files));  % +1 for header row

    for i = 1:numel(files)
        T = fileTables{i};
        outputSheet{1, i} = fileNames{i};

        if ismember(varName, T.Properties.VariableNames)
            col = T.(varName);
            col = col(:);
            colCells = num2cell(col);
            outputSheet(2:1+length(colCells), i) = colCells;
        end
    end

    % Truncate sheet name if needed
    sheetName = matlab.lang.makeValidName(varName);
    if strlength(sheetName) > 31
        sheetName = extractBefore(sheetName, 32);
    end

    % Write sheet
    writecell(outputSheet, outputFile, 'Sheet', sheetName);
end

disp("Successfully compiled output files");
