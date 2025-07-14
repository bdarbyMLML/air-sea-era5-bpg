% Calculating  east-west (u) and north-south (v) wind components 
% from wind speed and direction.
% Input wind direction is between [0 360], denoting direction from which
% wind is coming.
%
% Q. Wang, 8/1/2011
%
function [u,v]=wswd2uv(ws,wd)
DD=(wd+180.)/180*pi;   
u=ws.*sin(DD);
v=ws.*cos(DD);
return