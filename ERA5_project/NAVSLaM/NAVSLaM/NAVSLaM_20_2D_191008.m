function [ustar,tstar,qstar,l,edh,mmin,zprof,u,t,q,p,m,ctsq,cqsq,ctq,cnsq] = ...
    NAVSLaM_20_2D_191008(lambda,ws,tair,tsea,h,hflag,pr,s,lat,az, ...
    zu,zt,zh,zp,zinc,zmax)

% Navy Atmospheric Vertical Surface Layer Model (NAVSLaM) Version 2.0
% shell for 2-dimensional gridded data model runs.
%
% Release of 8 November 2019
%
% Input 2-dimensional arrays should have the dimensions (nlat,nlon), where
% nlat is the number of latitude (or y coordinate) points in the grid and
% nlon is the number of longitude (or x coordinate) points in the grid.
%
% The output 2-dimensional arrays can be contour plotted using a MATLAB
% statement of the form: contourf(lat,lon,cnsq(:,:,nht)), where lat is the 
% 1-D array of latitude values, lon is the 1-D array of longitude values 
% and nht is the array element for the height of interest to be plotted.
%
%  Input parameters:
%    lambda - Wavelength (microns) or radio frequency flag (lambda = 0)
%             Can be scalar or array for optical wavelength applications
%             Must be a scalar value of 0 for radio frequency applications
%    ws     - Wind speed (m/s)
%    tair   - Air temperature (C)
%    tsea   - Sea temperature (C)
%    h      - Humidity parameter, either relative or specific humidity
%    hflag  - Humidity flag, indicates humidity parameter for h input: 
%             1 = Relative humidity (%), 2 = Specific humidity (kg/kg)
%             Must be a scalar value
%    pr     - Atmospheric pressure (hPa) - default is 1013 hPa
%    s      - Ocean salinity (PSU)
%             Can be scalar or array input, default value is 35 PSU
%    lat    - Latitude in decimal degrees (positive north, negative south)
%             Can be scalar or array input, default value is 40 degrees
%    az     - Azimuth of propagation in degrees clockwise from north
%             Can be scalar or array input, default value is 45 degrees
%    zu     - Wind speed measurement/forecast height (m)
%             Can be scalar or array input
%    zt     - Air temperature measurement/forecast height (m)
%             Can be scalar or array input
%    zh     - Humidity measurement/forecast height (m)
%             Can be scalar or array input
%    zp     - Atmospheric pressure measurement/forecast height (m)
%             Can be scalar or array input
%    zinc   - Height array increment (scalar) for computing profiles (m)
%             or a pre-defined height array for output vertical profiles
%    zmax   - Maximum height for computing profiles (m) - scalar
%             Not used if zinc is a pre-defined height array
%
%  Output parameters:
%    ustar - Wind speed surface layer scaling parameter (m/s)
%    tstar - Potential temperature surface layer scaling parameter (C)
%    qstar - Specific humidity surface layer scaling parameter (kg/kg)
%    l     - Obukhov length scale (m)
%    edh   - Evaporation duct height (m)
%    mmin  - Modified refractivity value at the EDH (M-units)
%    zprof - Height above the surface array (m) corresponding to the
%            u, t, q, p, m, ctsq, cqsq, ctq and cnsq profile output arrays
%    u     - Vertical wind speed profile (m/s)
%    t     - Vertical air temperature profile (C)
%    q     - Vertical specific humidity profile (kg/kg)
%    p     - Vertical atmospheric pressure profile (hPa)
%    m     - Vertical refractivity profile for optical wavelength cases, or
%            vertical modified refractivity profile for radio frequency 
%            cases (no units)
%    ctsq  - Vertical temperature structure parameter profile (K^2 m^-2/3)
%    cqsq  - Vertical specific humidity structure parameter profile 
%            ((kg/kg)^2 m^-2/3)
%    ctq   - Vertical temperature-specific humidity cross structure 
%            parameter profile (K (kg/kg) m^-2/3)
%    cnsq  - Vertical refractive index structure parameter profile (m^-2/3)
%
% Developed and written by Paul Frederickson, Naval Postgraduate School

% Find size and length of input arrays: 
asize = size(ws) ; nlat = asize(1) ; nlon = asize(2) ; nlen = nlat*nlon ;

% Reshape original 2-D arrays into 1-D column vectors for input to NAVSLaM:
ws = reshape(ws,[nlen,1]) ;
tair = reshape(tair,[nlen,1]) ;
tsea = reshape(tsea,[nlen,1]) ;
h = reshape(h,[nlen,1]) ;
pr = reshape(pr,[nlen,1]) ;
if length(lambda) > 1 ; lambda = reshape(lambda,[nlen,1]) ; end
if length(s) > 1 ; s = reshape(s,[nlen,1]) ; end
if length(lat) > 1 ; lat = reshape(lat,[nlen,1]) ; end
if length(az) > 1 ; az = reshape(az,[nlen,1]) ; end
if length(zu) > 1 ; zu = reshape(zu,[nlen,1]) ; end
if length(zt) > 1 ; zt = reshape(zt,[nlen,1]) ; end
if length(zh) > 1 ; zh = reshape(zh,[nlen,1]) ; end
if length(zp) > 1 ; zp = reshape(zp,[nlen,1]) ; end

% Call NAVSLaM Version 2.0 with 1-D array inputs:
[ustar,tstar,qstar,l,edh,mmin,zprof,u,t,q,p,m,ctsq,cqsq,ctq,cnsq] = ...
    NAVSLaM_20_191008(lambda,ws,tair,tsea,h,hflag,pr,s,lat,az, ...
    zu,zt,zh,zp,zinc,zmax) ;

% Reshape NAVSLaM output 1-D column vectors back into 2-D lat-lon arrays:
n = length(zprof) ;
ustar = reshape(ustar,[nlat,nlon]) ; tstar = reshape(tstar,[nlat,nlon]) ;
qstar = reshape(qstar,[nlat,nlon]) ; l = reshape(l,[nlat,nlon]) ;
edh = reshape(edh,[nlat,nlon]) ; mmin = reshape(mmin,[nlat,nlon]) ;
m = reshape(m,[nlat,nlon,n]) ; t = reshape(t,[nlat,nlon,n]) ;
q = reshape(q,[nlat,nlon,n]) ; p = reshape(p,[nlat,nlon,n]) ;
u = reshape(u,[nlat,nlon,n]) ;
ctsq = reshape(ctsq,[nlat,nlon,n]) ; cqsq = reshape(cqsq,[nlat,nlon,n]) ;
ctq = reshape(ctq,[nlat,nlon,n]) ; cnsq = reshape(cnsq,[nlat,nlon,n]) ;

end % End of function NAVSLaM_20_2D_191008