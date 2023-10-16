
function desiredPaths = extract_paths(file_list)

% Preallocate the output cell array for efficiency
desiredPaths = cell(size(file_list, 1), 1);

% Loop through each path in the list
for i = 1:size(file_list, 1)
    % Split the path into its components
    parts = strsplit(strtrim(file_list(i, :)), '\');
    
    % Extract the desired component
    desiredPath = parts{end-2}; % Based on your example, it's the third-to-last component
    
    % Store in the output cell array
    desiredPaths{i} = ['\', desiredPath];
end

% % Convert cell array back to char array if needed
% desiredPaths_char = char(desiredPaths);


end 
