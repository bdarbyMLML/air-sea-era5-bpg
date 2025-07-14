% Calculating wind speed and direction from east-west (u) and north-south
% (v) wind components
% Output wind direction is between [0 360], denoting direction from which
% wind is coming.
%
% Q. Wang, 8/1/2011
%
function [ws,wd]=uv2wswd(u,v)
wd=180.-atan2(-u,v)*180./pi;
id=find(wd<0);
wd(id)=wd(id)+360;
id=find(wd>360);
wd(id)=wd(id)-360;
ws=sqrt(u.^2+v.^2); 
return