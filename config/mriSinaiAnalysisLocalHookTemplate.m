function mriSinaiAnalysisLocalHook
% mriSinaiAnalysisLocalHook
%
% For use with the ToolboxToolbox.  If you copy this into your
% ToolboxToolbox localToolboxHooks directory (by default,
% ~/localToolboxHooks) and delete "LocalHooksTemplate" from the filename,
% this will get run when you execute tbUse({'mriSinaiAnalysisConfig'}) to set up for
% this project.  You then edit your local copy to match your local machine.
%
% The main thing that this does is define Matlab preferences that specify input and output
% directories.
%
% You will need to edit the project location and i/o directory locations
% to match what is true on your computer.

 
%% Define project
projectName = 'mriSinaiAnalysis';
 

%% Clear out old preferences
if (ispref(projectName))
    rmpref(projectName);
end


%% Specify and save project location
projectBaseDir = tbLocateProject(projectName);
setpref(projectName,'projectBaseDir',projectBaseDir);


%% Get the userID
[~, userID] = system('whoami');
userID = strtrim(userID);



end
