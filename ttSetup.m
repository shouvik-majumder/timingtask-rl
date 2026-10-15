function root = ttSetup()
%TTSETUP  Add the library folders to the MATLAB path and return the repository root.
%
%   root = ttSetup()
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src', 'data'), fullfile(root, 'src', 'dist'), fullfile(root, 'src', 'policy'), ...
    fullfile(root, 'src', 'groups'), fullfile(root, 'src', 'plots'));
end
