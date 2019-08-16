function GOP = GOPdate()

% GOP = GOPdate()
%
% Return the current date/time in GOP format.
%
% "GOP" stands for the "Grand Observing Program" of the UCSD IPS Lab
% (Interplanetary Scintillation Radio Observatory.)
%
% These strings sort the same numerically, alphabetically, and temporally.
%
% In case you need it in your OS, here is how to get these strings 
% using the Posix date command...
% Local time:
% date "+%Y%m%d_%H%M%Z"  # Leave off the %Z to drop the TZ.
% UTC:
% date -u "+%Y%m%d_%H%M%S%Z"
% 
% Johanan L. Codona: 20050721

GOP = datestr(now,'yyyymmdd_HHMMSS');
