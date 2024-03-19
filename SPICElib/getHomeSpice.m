function HOMESPICE = getHomeSpice
%Change this according to your system
% IMPORTANT:
% 1-The folder 'imgo' should be created inside the path returned by this
% function
% e.g. '/Volumes/Extreme SSD/SPICE/imgo'
% 2-The data files with image information have to be COPIED inside the 
% path returned by this function
% e.g. '/Volumes/Extreme SSD/SPICE/CASSINI'

[~, name] = system('hostname');
name=strtrim(name);
%HOMESPICE='/Volumes/Extreme SSD/SPICE/';
if strcmp(name,'MacBook-Pro-de-Paula-2.local') || ...
        strcmp(name,'10-192-55-128client.eduroam.upc.edu') ||  ...
        strcmp(name,'MBP-de-Paula-2.homenet.telecomitalia.it')
    HOMESPICE='/Users/paulabetriu/Desktop/MASTER/MUEA/TFM/SPICE/MICE/mice/';
else
    error("user not identified. Please update the SPICE MICE path in getHomeSpice")
end

end

