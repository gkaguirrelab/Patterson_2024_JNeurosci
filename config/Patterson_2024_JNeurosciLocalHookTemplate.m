function Patterson_2024_JNeurosciLocalHook
% Patterson_2024_JNeurosciLocalHook
%
% For use with the ToolboxToolbox.  If you copy this into your
% ToolboxToolbox localToolboxHooks directory (by default,
% ~/localToolboxHooks) and delete "LocalHooksTemplate" from the filename,
% this will get run when you execute tbUse({'Patterson_2024_JNeurosciConfig'}) to set up for
% this project.  You then edit your local copy to match your local machine.
%
% The main thing that this does is define Matlab preferences that specify input and output
% directories.
%
% You will need to edit the project location and i/o directory locations
% to match what is true on your computer.

 
%% Define project
projectName = 'Patterson_2024_JNeurosci';
 

%% Clear out old preferences
if (ispref(projectName))
    rmpref(projectName);
end


%% Paths to find DropBox directory
if ismac
    [~, userName] = system('whoami');
    userName = strtrim(userName);
elseif isunix
    userName = getenv('USER');
elseif ispc
    userName = getenv('username');
else
    disp('What are you using?')
end

dropboxBaseDir = ...
    fullfile('/Users', userName, ...
    'Dropbox (Aguirre-Brainard Lab)');

% Set pref to dropboxBaseDir
setpref(projectName,'dropboxBaseDir',dropboxBaseDir);

% Set pref to username
setpref(projectName,'userName',userName);


%% Specify and save project location
projectBaseDir = tbLocateProject(projectName);
setpref(projectName,'projectBaseDir',projectBaseDir);

%% Flywheel key preferences 
setpref(projectName, 'flywheelAPIkey', 'Enter your key here');

%% Check for required Matlab toolboxes
% The set of Matlab add-on toolboxes being used can be determined by
% running the routines in the projecy, followed by the license function.
%{
    license('inuse')
%}
% This provides a list of toolbox license names. In the following
% assignment, the license name is given in the comment string after the
% matching version name for each toolbox.
requiredAddOns = {...
    'Parallel Computing Toolbox',...              % distrib_computing_toolbox
    'Statistics and Machine Learning Toolbox',...   % statistics_toolbox
    'Symbolic Math Toolbox',...                 % symbolic_toolbox
    };
% Given this hard-coded list of add-on toolboxes, we then check for the
% presence of each and issue a warning if absent.
V = ver;
VName = {V.Name};
warnState = warning();
warning off backtrace
for ii=1:length(requiredAddOns)
    if ~any(strcmp(VName, requiredAddOns{ii}))
        warnString = ['The Matlab ' requiredAddOns{ii} ' is missing. ' toolboxName ' may not function properly.'];
        warning('localHook:requiredMatlabToolboxCheck',warnString);
    end
end
warning(warnState);


end
