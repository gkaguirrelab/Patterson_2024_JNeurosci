function stimulusContrastScale = returnStimulusContrastScale(cellClass,stimulusDirection)

% The RGC data presented by Solomon, Lee, and colleagues was expressed as
% gain, in units of spikes / sec / % contrast. Percent contrast was defined
% in a relative sense
% as the percent of the maximum available contrast available in the
% particular stimulus direction. These were:
%   LMS: 100%
%   L-M: 12%
%     S: 80%
%
% So, an actual 50% S-directed modulation is characterized as a 50/86 =
% 58% relative contrast in this direction.
%
% The gain value was estimated from the linear portion of a Naka-Rushton
% function fit to the contrast response data. For the bistratified and
% midget RGCs, we will assume that the stimuli we are presenting are within
% the linear range of the response. Parasol cells, however, saturate at 25%
% actual / relative contrast; we account for that property here.


% Obtain the chromatic weights.
switch cellClass
    case 'midget'
        switch stimulusDirection
            case 'LminusM'
                stimulusContrastScale = 100*(8/12);
            case 'LMS'
                stimulusContrastScale = 90;
            case 'S'
                stimulusContrastScale = nan;
        end
    case 'parasol'
        switch stimulusDirection
            case 'LminusM'
                stimulusContrastScale = nan;
            case 'LMS'
                stimulusContrastScale = 25;
            case 'S'
                stimulusContrastScale = nan;
        end
    case 'bistratified'
        switch stimulusDirection
            case 'LminusM'
                stimulusContrastScale = nan;
            case 'LMS'
                stimulusContrastScale = nan;
            case 'S'
                stimulusContrastScale = 100*(50/86);
        end
end


end
