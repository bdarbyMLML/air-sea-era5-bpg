function [ustar,tstar,qstar,l,edh,mmin,zprof,u,t,q,p,m,ctsq,cqsq,ctq,cnsq] = ...
    NAVSLaM_20_191008(lambda,ws,tair,tsea,h,hflag,pr,s,lat,az,...
    zu,zt,zh,zp,zinc,zmax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Navy Atmospheric Vertical Surface Layer Model (NAVSLaM), Version 2.0   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Release of 8 November 2019
%
%    NAVSLaM is a bulk model for computing surface-layer scaling 
%  parameters (u*, T*, q*, L), the evaporation duct height (for RF cases) 
%  and the vertical profiles of refractivity, air temperature, specific 
%  humidity, pressure, wind speed, and various structure parameters (CT2, 
%  Cq2, CTq, Cn2) from mean values of meteorological quantities measured or 
%  forecasted near the ocean surface and the sea surface temperature and 
%  salinity.  The input wind speed, air temperature and humidity input 
%  values must be measured or forecasted at heights between 0.5 and 50 m 
%  above the ocean surface.  NAVSLaM v2.0 is valid for radio-frequency and 
%  optical wavelength (0.3 to 14 micron) applications.  Note that this 
%  model is only valid for applications above the ocean or other large 
%  bodies of water.
%
%    The meteorological inputs (u, tair, tsea, h, pr) can be either single 
%  point (scalar) or 1D array inputs.  Other inputs (s, lat, az, zu, zt,
%  zh, zp) can be scalar or array inputs, regardless of whether the 
%  meteorological inputs are scalars or arrays.  For example, if all
%  wind speed data are for the same height, then zu can be input as a 
%  scalar value, rather than as an array of values.  If using input arrays, 
%  all array dimensions and shapes must be the same.
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
%    All output parameters are set to -99 if any input data are outside of
%  their allowed limits.
%    All output parameters  are set to -98 for cases when the iteration 
%  does not converge (which should not occur).
%    edh is set to -97 if the air temperature and/or specific humidity
%  profiles reach unrealistic values (i.e. T(edh) > 55, q(edh) <= 0).
%    edh is set to -96 if the evaporation duct height is undefined
%  (i.e. there is no local minimum m value below the top height zmax).         
%
%  NAVSLaM Version 2.0 was developed and written by:
%    Paul A. Frederickson
%    Department of Meteorology
%    Naval Postgraduate School
%    589 Dyer Rd, Room 254
%    Monterey, CA 93943 
%    Phone:  (831) 521-8670
%    E-mail:  pafreder@nps.edu


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Section to Define Input Dimensions and Check Intput Data Quality       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find size and length of input arrays: 
asize = size(ws) ; iscolumnu = iscolumn(ws) ; nlen = length(ws) ;

% Initially set all output parameters to -99:
bad = -99*ones(asize) ;
ustar = bad ; tstar = bad ; qstar = bad ; l = bad ; edh = bad ;
t = bad ; q = bad ; p = bad ; m = bad ; cnsq = bad ; zprof = bad ; 
u = bad ; ctsq = bad ; cqsq = bad ; ctq = bad ; mmin = bad ;
clear bad

% Check input lambda parameter:
if length(lambda) == 1
    if (lambda < 0 || (lambda > 0 && lambda < 0.3) || lambda > 14)
        fprintf(['Error in lambda input: lambda = ' ...
            num2str(lambda) ' microns\n'])
        return
    end
elseif length(lambda) > 1
    ibad = find(lambda < 0 | (lambda > 0 & lambda < 0.3) | lambda > 14);
    if ~isempty(ibad)
        lambda(ibad) = NaN ; clear ibad
        fprintf('Warning! Some lambda inputs are bad, set to -99...\n')
    end ; clear ibad
    if length(lambda(lambda == 0)) == length(lambda)
        clear lambda
        lambda = 0 ;
    end
end

% Check input hflag humidity indicator:
if length(hflag) > 1 || isempty(hflag)
    fprintf(['Error! Input hflag value is invalid.  ' ...
        'Must be a scalar value of 1 or 2...\n'])
    return
else
    if hflag < 1 || hflag > 2
        fprintf(['Error! Input hflag value is invalid.  ' ...
            'Must be a scalar value of 1 or 2...\n'])
        return
    end
end

% Section to define and check the vertical profile height array:
if length(zinc) == 1 % If zinc is a scalar value
% Use the user-defined zinc and zmax values to compute height array
    if length(zmax) > 1 || isempty(zmax)
        fprintf('Error! zmax input must be a scalar value...\n')
        return
    end
    if isnan(zinc) || zinc <= 0 || isnan(zmax) || zmax <= 0 || ...
        zinc > zmax
        fprintf(['Error in zinc and/or zmax inputs: zinc = ' ...
            num2str(zinc) ' m, zmax = ' num2str(zmax) ' m\n'])
        return
    end
    if length(lambda) == 1 && lambda == 0 && zinc > 0.5
        fprintf(['Warning! zinc value is too large to properly resolve' ...
            ' edh. zinc should be <= 0.5 m.\n'])
    end
% Maximum height cannot be greater than 100 m
    if zmax > 100
        fprintf('Maximum height greater than 100 m, reseting to 100 m.\n')
        zmax = 100 ;
    end
% Define height arrays:
    zprof = 0 : zinc : zmax ; nprof = length(zprof) ;
else % If zinc is an array value
% Use the user-defined input height array, but check for validity
    zprof = zinc ; nprof = length(zprof) ; zmax = zprof(nprof) ;
    if zprof(1) ~= 0
        fprintf('Input height array does not begin at the surface...\n')
        return
    end
    ibad = find(zprof < 0,1) ;
    if ~isempty(ibad)
        fprintf('Input height array cannot have negative values...\n')
        return
    end ; clear ibad
    ibad = find(zprof > 100,1) ;
    if ~isempty(ibad)
        fprintf('Input height array cannot have values above 100 m...\n')
        return
    end ; clear ibad
    iinc = find(zprof(2:nprof)-zprof(1:nprof-1) <= 0,1) ;
    if ~isempty(iinc)
        fprintf('Input height array does not monotonically increase...\n')
        return
    end ; clear iinc
end
if iscolumn(zprof) ; zprof = zprof' ; end
if nlen > 1 ; z = repmat(zprof,nlen,1) ; else ; z = zprof ; end

% Set output profile arrays to -99:
if iscolumnu || nlen == 1
    m = -99*ones(nlen,nprof) ;
else
    m = -99*ones(nprof,nlen) ;
end
u = m ; t = m ; q = m ; p = m ; cnsq = m ; ctsq = m ; cqsq = m ; ctq = m ;

% Check that dimensions of input arrays match, if not then exit:
if ~((asize == size(tair) & asize == size(tsea) & asize == size(h) & ...
      asize == size(pr) & ...
    ((asize == size(lambda) & length(lambda) > 1) | length(lambda) == 1) & ...
    ((asize == size(s) & length(s) > 1) | length(s) == 1) & ...
    ((asize == size(lat) & length(lat) > 1) | length(lat) == 1) & ...
    ((asize == size(az) & length(az) > 1) | length(az) == 1) & ...
    ((asize == size(zu) & length(zu) > 1) | length(zu) == 1) & ...
    ((asize == size(zt) & length(zt) > 1) | length(zt) == 1) & ...
    ((asize == size(zh) & length(zh) > 1) | length(zh) == 1) & ...
    ((asize == size(zp) & length(zp) > 1) | length(zp) == 1)))
    fprintf('Error!  Input array sizes do not match...\n')
    return
end

% If input arrays are row vectors, reshape to column vectors:
if ~iscolumnu
    fprintf('Inputs are row vectors, reshaping to column vectors...\n')
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
end
asize2 = size(ws) ;

% Check input parameter values:
% Set out-of-bounds secondary input parameters to default values:
ibad = find(pr < 920 | pr > 1060 | isnan(pr)) ;
if ~isempty(ibad) ; pr(ibad) = 1013.25 ; 
    fprintf('Warning! Some pressure inputs are bad, set to 1013 hPa...\n')
end ; clear ibad
ibad = find(zp < 0 | zp > 50 | isnan(zp)) ;
if ~isempty(ibad) ; zp(ibad) = 0 ; 
    fprintf('Warning! Some zp inputs are bad, set to 0 m...\n')
end ; clear ibad
ibad = find(s <= 0 | s > 60) ;
if ~isempty(ibad) ; s(ibad) = 35 ; 
    fprintf('Warning! Some salinity inputs are bad, set to 35 PSU...\n')
end ; clear ibad
ibad = find(abs(lat) > 90) ;
if ~isempty(ibad) ; lat(ibad) = 40 ; 
    fprintf('Warning! Some latitude inputs are bad, set to 40 deg...\n')
end ; clear ibad
ibad = find(az < 0 | az > 360) ;
if ~isempty(ibad) ; az(ibad) = 45 ; 
    fprintf('Warning! Some azimuth inputs are bad, set to 45 deg...\n')
end ; clear ibad
% Set out-of-bounds primary input parameters to NaN:
ws(ws >= 0 & ws < 0.1) = 0.1 ;
ibad = find(ws < 0 | ws > 30) ;
if ~isempty(ibad) ; ws(ibad) = NaN ; 
    fprintf('Warning! Some wind speed inputs are bad, set to -99...\n')
end ; clear ibad
ibad = find(abs(tair) > 50) ;
if ~isempty(ibad) ; tair(ibad) = NaN ; 
    fprintf('Warning! Some air temp inputs are bad, set to -99...\n')
end ; clear ibad
ibad = find(tsea < -2.6 | tsea > 36) ;
if ~isempty(ibad) ; tsea(ibad) = NaN ; 
    fprintf('Warning! Some sea temp inputs are bad, set to -99...\n')
end ; clear ibad
if length(zu) == 1
    if (zu < 0.5 || zu > 50)
        fprintf('Error! Wind speed height input is bad...\n')
        return
    end
else
    ibad = find(zu < 0.5 | zu > 50) ;
    if ~isempty(ibad) ; zu(ibad) = NaN ; ws(ibad) = NaN ; 
        fprintf('Warning! Some zu inputs are bad, set to -99...\n')
    end ; clear ibad
end
if length(zt) == 1
    if (zt < 0.5 || zt > 50)
        fprintf('Error! Air temperature height input is bad...\n')
        return
    end
else
    ibad = find(zt < 0.5 | zt > 50) ;
    if ~isempty(ibad) ; zt(ibad) = NaN ; tair(ibad) = NaN ; 
        fprintf('Warning! Some zt inputs are bad, set to -99...\n')
    end ; clear ibad
end
if length(zh) == 1
    if (zh < 0.5 || zh > 50)
        fprintf('Error! Humidity height input is bad...\n')
        return
    end
else
    ibad = find(zh < 0.5 | zh > 50) ;
    if ~isempty(ibad) ; zh(ibad) = NaN ; h(ibad) = NaN ; 
        fprintf('Warning! Some zh inputs are bad, set to -99...\n')
    end ; clear ibad
end
% Check humidity values:
if hflag == 1 % Relative humidity (%)
    ibad = find(h < 10 | h > 100) ;
    if ~isempty(ibad) ; h(ibad) = NaN ;
        fprintf('Warning! Some humidity inputs are bad, set to -99...\n')
    end ; clear ibad
elseif hflag == 2 % Specific humidity (kg/kg)
    ibad = find(h <= 0 | h > 0.1) ;
    if ~isempty(ibad) ; h(ibad) = NaN ; 
        fprintf('Warning! Some humidity inputs are bad, set to -99...\n')
    end ; clear ibad
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Section to Define Required Physical Parameters                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Atmospheric carbon dioxide concentration (ppm)
CO2 = 400 ;        
% Molar mass of dry air (kg/mol)
Md = 0.02896546 + 0.012011*(CO2*1e-6 - 0.0004) ;
% Molar mass of water vapor (kg/mol)
Mw = 0.018015268 ;  
% Ratio of molar mass of water vapor to dry air
eps = Mw./Md ;
% Virtual temperature equation constant
gam = (1 - eps)./eps ;
% Universal ideal gas constant (J/mol/K)
R = 8.3144598 ;     
% Ideal gas constant for dry air (J/kg/K)
Rd = R/Md ;
% Specific heat of dry air at constant pressure (J/kg/K):
cp = 1005.60 + 0.017211*tair + 0.000392*tair.*tair ;
% Kinematic viscocity of dry air (m^2/s) (Andreas 2005):
nu = 1.326e-5*(1 + 6.542e-3*tair + 8.301e-6*tair.^2 - 4.84e-9*tair.^3) ;
% Latent heat of vaporization (J/Kg):
%Lv = (25.0 - 0.02274*tair)*1e5 ;
% Celsius-Kelvin conversion factor (K)
T0 = 273.15 ;       
% von Karman constant
k = 0.4 ;           
% Charnock's constant, 0.011 for open ocean (Smith 1988)
alpha = 0.011 ;     
% Free convective parameter
beta = 1.25 ;       
% Estimated inversion height (m)
zi = 600 ;          
% Define acceleration due to gravity (m/s/s)
g = gravity(lat) ;
if length(g) == 1 && nlen > 1 ; g = g*ones(asize2) ; end

% Compute specific humidity, if input is relative humidity:
% Estimate pressure at specific humidity measurement height
p_zh = pr.*exp(g.*(zp-zh)./Rd./(tair + T0)) ;
if hflag == 1
    qair = spec_hum(tair,h,p_zh,eps) ; clear p_zh
elseif hflag == 2
    qair = h ; 
    qs = spec_hum(tair,100,p_zh,eps) ; clear p_zh
    ibad = find(qair > qs) ;
    if ~isempty(ibad) ; qair(ibad) = NaN ;
        fprintf(['Warning! Some humidity inputs are greater than ' ...
            'saturation values, set to -99...\n'])
    end ; clear ibad
end

% Estimate surface pressure:
tv = (tair + T0).*(1 + gam.*qair) ; % Compute virtual temperature
psfc = pr.*exp(g.*zp./Rd./tv) ;     % From hypsometric equation
clear tv

% Estimate sea surface specific humidity:
rh_sea = (1 - 0.000537*s)*100;    % Takes into account ocean salinity
qsea = spec_hum(tsea,rh_sea,psfc,eps);

% Compute temperatures:
tabs = tair + T0;                 % Absolute air temperature (K)
tsabs = tsea + T0;                % Absolute sea temperature (K)
theta = tabs + (g./cp).*zt;       % Potential temperature at height zt (K)
thetav = theta.*(1 + gam.*qair);  % Virtual potential temperature at zt (K)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Section to Iteratively Determine Surface Layer Scaling Parameters      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define first-guess values for iteration and pre-allocate arrays:
wstar = zeros(asize2);
us = sqrt(ws.*ws + 0.01) ;
ustaro = us.*sqrt((0.6 + 0.07*us)/1000) ;
tstaro = 0.0011*(us./ustaro).*(theta - tsabs) ;
qstaro = 0.0011*(us./ustaro).*(qair - qsea) ; clear us
tstarv = tstaro.*(1 + gam.*qair) + gam.*tabs.*qstaro ;
lo = thetav.*ustaro.*ustaro./(k*g.*tstarv) ;

% Iteration loop to determine the surface layer scaling parameters 
% (ustar, tstar, qstar and l)
% Iterate up to 25 times, if necessary
counter = 0; % Counter for number of iterations
aconv = 1 ; conv = 0.001 ;
while (counter < 25 && aconv > conv)
%while (counter < 25)
    counter = counter + 1;
    zo = alpha*ustaro.*ustaro./g + 0.11*nu./ustaro;
%   zot = 5.5e-5*(zo.*ustaro./nu).^(-0.63) ;             % COARE 2.5 eqn
    zot = min(1.15e-4,5.5e-5*(zo.*ustaro./nu).^(-0.6)) ; % COARE 3.0 eqn
    if length(g) > 1
        wstar(tstarv < 0) = ...
            (-zi.*g(tstarv < 0).*ustaro(tstarv < 0).*tstarv(tstarv < 0)./ ...
            thetav(tstarv < 0)).^(1/3);
    else
        wstar(tstarv < 0) = ...
            (-zi.*g.*ustaro(tstarv < 0).*tstarv(tstarv < 0)./ ...
            thetav(tstarv < 0)).^(1/3);
    end
    us = sqrt(ws.*ws + (beta*wstar).*(beta*wstar));
    ustar = us*k./(log(zu./zo) - psi_m(zu./lo));
    tstar = (theta - tsabs).*k./(log(zt./zot) - psi_t(zt./lo));
    qstar = (qair - qsea).*k./(log(zh./zot) - psi_t(zh./lo));
    tstarv = tstar.*(1 + gam.*qair) + gam.*tabs.*qstar;
    l = thetav.*ustar.*ustar./(k*g.*tstarv);
    checku = abs((ustar - ustaro)./ustar) ;
    checkt = abs((tstar - tstaro)./tstar) ;
    checkq = abs((qstar - qstaro)./qstar) ;
    checkl = abs((l - lo)./l) ;
    ustaro = ustar; tstaro = tstar ; qstaro = qstar ; lo = l;
    aconv = max([max(checkl),max(checku),max(checkt),max(checkq)]) ;
    wstar = zeros(asize2);
end


% Update roughness lengths with final ustar values:
zo = alpha*ustar.*ustar./g + 0.11*nu./ustar;
zot = min(1.15e-4,5.5e-5*(zo.*ustar./nu).^(-0.6)) ; % COARE 3.0 eqn

clear nu ustaro tstaro qstaro lo us tstarv aconv


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Section to Compute Vertical Profiles and Evaporation Duct Height       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nlen > 1
    tstarm = repmat(tstar,1,nprof) ; qstarm = repmat(qstar,1,nprof) ;
else
    tstarm = tstar ; qstarm = qstar ; 
end

% Pre-allocate vertical profile arrays:
t = zeros(nlen,nprof) ; q = zeros(nlen,nprof) ; p = zeros(nlen,nprof) ;
u = zeros(nlen,nprof) ; gt = zeros(nlen,nprof) ;
% Define surface values:
u(:,1) = 0 ;          % Surface temperature
t(:,1) = tsabs ;      % Surface temperature
q(:,1) = qsea ;       % Surface specific humidity
p(:,1) = psfc ;       % Surface pressure
gt(:,1) = NaN*ones(size(tsabs)) ; % Dim temp structure parameter function

% Loop to compute wind speed, temperature, humidity and pressure profiles:
for i = 2 : nprof
    zz = zprof(i) ;
    u(:,i) = ws + ustar.*(log(zz./zu) - psi_m(zz./l) + psi_m(zu./l))/k;
    br = (log(zz./zot) - psi_t(zz./l))/k ;
    t(:,i) = tsabs + tstar.*br - g.*zz./cp;
    q(:,i) = qsea + qstar.*br;
    tv_avg = (t(:,i).*(1 + gam*q(:,i)) + t(:,i-1).*(1 + gam*q(:,i-1)))/2 ;
    zincr = zprof(i) - zprof(i-1) ;
    p(:,i) = p(:,i-1).*exp(-g.*zincr/Rd./tv_avg) ;
    gt(:,i) = dim_ctsq(zz./l) ; % Dimensionless temp structure parameter
end

% Compute vapor pressure profile (hPa):
qbr = (eps + (1-eps)*q) ;
e = q.*p./qbr ;

% Compute refractivity or modified refractivity profiles and
% the evaporation duct height for radio frequencies:
if length(lambda) == 1 && lambda == 0 % Radio frequencies
    
    % Radio refractivity equation empirical constants:
    %k1 = 77.6 ; k2 = -5.6 ; k3 = 375000 ; Re_inv = 0.157 ; 
    %k1 = 77.6904 ; k2 = -6.39522 ; k3 = 375463 ; Re_inv = 0.157 ;
    c1 = 77.6681 + CO2*1e-6*(133.48 - 77.6681) ; % K/hPa
    c2 = 71.2952 - c1 ; % K/hPa
    c3 = 3.75463e5 ; % K^2/hPa

    % Compute radius of curvature of the earth
    rcurv = radius_curv(lat,az) ;
    if length(lat) > 1 || length(az) > 1 
        rcurv = repmat(rcurv,1,nprof) ;
    end
    
    % Compute vertical modified refractivity profile:
    m = c1.*p./t + c2.*e./t + c3.*e./t./t + 1e6*z./rcurv ;

    %  Find the evaporation duct height:
    [mmin,iduct] = min(m,[],2) ;
    edh = zprof(iduct) ;

    % If the evaporation duct height is not defined, set = -96:
    mmin(edh == zmax) = -96 ; edh(edh == zmax) = -96 ; 

else % Optical wavelengths

    % Define optical wavelength refractivity empirical constants:
    [m1, m2] = Refractivity_wavelength_coeffs(lambda, CO2) ;
    c1 = m1 ; c2 = m2 - m1 ; c3 = 0 ;
    if length(lambda) > 1
        c1 = repmat(c1,1,nprof) ; c2 = repmat(c2,1,nprof) ; 
        c3 = zeros(size(c1)) ;
    end
    
    % Compute vertical refractivity profile:
    m = c1.*p./t + c2.*e./t ;
    
    % Set EDH and mmin to missing:
    edh = -99*ones(flip(asize)) ; mmin = -99*ones(flip(asize)) ;
    
end

% Compute refractive index derivatives:
A = -1e-6*p./t./t.*(c1 + c2.*q./qbr + 2*c3.*q./qbr./t) ; % A = dn/dT
B = 1e-6*p./t./qbr.*(1 - (1-eps).*q./qbr).*(c2 + c3./t) ; % B = dn/dq
%C = 1e-6*(c1 + c2*q./qbr + c3*q./t./qbr)./t ; % C = dn/dP

% Convert output temperature profile values from Kelvin to degrees Celsius:
t = t - T0 ;

% Define the temperature-humidity fluctuation correlation:
%rtq = 0.75*ones(size(tstarm)) ;
%rtq(tstarm./qstarm < 0) = -0.5 ;
cp = 1005.60 + 0.017211*t + 0.000392*t.*t ; % Specific heat of dry air
Lv = (25.0 - 0.02274*t)*1e5 ; % Latent heat of vaporization
Bo = cp.*tstarm./Lv./qstarm ; % Compute the Bowen ratio
rtq = 0.75*ones(size(t)) ;
rtq(Bo < -0.75) = -0.55 ;
rtq(Bo >= -0.75 & Bo < 0) = Bo(Bo >= -0.75 & Bo < 0)*0.55/0.75 ;
rtq(Bo >= 0 & Bo < 0.5) = Bo(Bo >= 0 & Bo < 0.5)*0.75/0.5 ;

% Compute the dimensionless temperature structure parameter function:
%gt = dim_ctsq(z./lm) ; 

% Compute the vertical structure parameter profiles:
ctsq = gt.*z.^(-2/3).*tstarm.*tstarm ;
cqsq = gt.*z.^(-2/3).*qstarm.*qstarm ;
ctq =  rtq.*gt.*z.^(-2/3).*tstarm.*qstarm ;
cnsq = A.*A.*ctsq + 2*A.*B.*ctq + B.*B.*cqsq ;
%cnsq = gt.*z.^(-2/3).*(A.*A.*tstarm.*tstarm + ...
%    2*rtq.*A.*B.*tstarm.*qstarm + B.*B.*qstarm.*qstarm);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Section to Check Output Array Dimensions and Check for Bad Outputs     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If optical wavelength is bad, set m & cnsq to -99:
m(isnan(m)) = -99 ; mmin(isnan(m)) = -99 ; cnsq(isnan(cnsq)) = -99 ;

% Set surface values of structure parameters to -99:
if zprof(1) == 0
    ctsq(:,1) = -99 ; cqsq(:,1) = -99 ; ctq(:,1) = -99 ; cnsq(:,1) = -99 ;
end

% Check vertical profiles for bad values & set bad values to -97:
qs = spec_hum(t,100,p,eps) ; % Saturation specific humidity
ibad = find(t > 50 | q < 0 | q > qs) ;
if ~isempty(ibad)
    t(ibad) = -97 ; q(ibad) = -97 ; p(ibad) = -97 ; m(ibad) = -97 ;
    u(ibad) = -97 ; 
    cnsq(ibad) = -97 ; ctsq(ibad) = -97 ; cqsq(ibad) = -97 ; 
    ctq(ibad) = -97 ;
end ; clear ibad
% For radio-frequency cases only:
if length(lambda) == 1 && lambda == 0
% If any profile values from surface to edh = -97, then set edh = -97
    for i = 1 : nlen
        if edh(i) >= 0
            if min(t(i,1:iduct(i)+1)) == -97
                edh(i) = -97 ; mmin(i) = -97 ; 
            end
        else
            if min(t(i,1:iduct(i))) == -97
                edh(i) = -97 ; mmin(i) = -97 ;
            end
        end
    end
end

% Set output for cases with bad input data = -99:
ibad = isnan(l) ;
if ~isempty(ibad)
    ustar(ibad) = -99 ; tstar(ibad) = -99 ; qstar(ibad) = -99 ;
    l(ibad) = -99 ; edh(ibad) = -99 ; u(ibad,:) = -99 ; 
    t(ibad,:) = -99 ; q(ibad,:) = -99 ; p(ibad,:) = -99 ; 
    m(ibad,:) = -99 ; cnsq(ibad,:) = -99 ; ctsq(ibad,:) = -99 ;
    cqsq(ibad,:) = -99 ; ctq(ibad,:) = -99 ; mmin(ibad,:) = -99 ;
end ; clear ibad

% Set output for cases which did not converge = -98:
ibad = find(checku > conv | checkt > conv | checkq > conv | checkl > conv);
if ~isempty(ibad)
    ustar(ibad) = -98 ; tstar(ibad) = -98 ; qstar(ibad) = -98 ; 
    l(ibad) = -98 ; edh(ibad) = -98 ; u(ibad,:) = -98 ;
    t(ibad,:) = -98 ; q(ibad,:) = -98 ; p(ibad,:) = -98 ;
    m(ibad,:) = -98 ; cnsq(ibad,:) = -98 ; ctsq(ibad,:) = -98 ;
    cqsq(ibad,:) = -98 ; ctq(ibad,:) = -98 ; mmin(ibad,:) = -98 ;
end ; clear ibad

% Convert output vectors back to original input vector shape, if necessary:
if ~iscolumnu
    ustar = reshape(ustar,[1,nlen]) ; tstar = reshape(tstar,[1,nlen]) ;
    qstar = reshape(qstar,[1,nlen]) ; l = reshape(l,[1,nlen]) ;
    edh = reshape(edh,[1,nlen]) ; u = reshape(u',[nprof,nlen]) ;
    m = reshape(m',[nprof,nlen]) ; t = reshape(t',[nprof,nlen]) ;
    q = reshape(q',[nprof,nlen]) ; p = reshape(p',[nprof,nlen]) ;
    cnsq = reshape(cnsq',[nprof,nlen]) ; ctsq = reshape(ctsq',[nprof,nlen]) ;
    cqsq = reshape(cqsq',[nprof,nlen]) ; ctq = reshape(ctq',[nprof,nlen]) ;
    zprof = zprof' ; mmin = reshape(mmin,[1,nlen]) ;
else
    edh = edh' ; mmin = mmin' ;
end

end % End of function NAVSLaM_20_191008


function [psim] = psi_m(zeta)

% Computes the dimensionless profile function for wind speed
% Input: zeta = z/L, where z is the height and L is the Obukhov length

% Pre-allocate array to neutral value:
psim = ones(size(zeta)) ;

%  Unstable case (zeta < 0):
x = (1 - 16*zeta(zeta < 0)).^(1/4);
psim(zeta < 0) = 2*log((1 + x)/2) + log((1 + x.*x)/2) - 2*atan(x) + pi/2;

%  Stable cases (zeta > 0):
y = (1 + zeta(zeta > 0)).^(1/3) ;
am = 5 ; bm = am/6.5 ; Bm = ((1 - bm)/bm)^(1/3) ; sr3 = sqrt(3) ;
psim(zeta > 0) = -3*am*(y-1)/bm + (am*Bm/2/bm)*(2*log((y+Bm)/(1+Bm)) - ...
    log((y.*y - y*Bm + Bm^2)/(1 - Bm + Bm^2)) + ...
    2*sr3*(atan((2*y-Bm)/sr3/Bm) - atan((2-Bm)/sr3/Bm))) ;

end % End of function psi_m


function [psit] = psi_t(zeta)

% Computes the dimensionless profile function for potential temperature
% Input: zeta = z/L, where z is the height and L is the Obukhov length

% Pre-allocate array to neutral value:
psit = ones(size(zeta)) ;

%  Unstable case (zeta < 0):
x = (1 - 16*zeta(zeta < 0)).^(1/2);
psit(zeta < 0) = 2*log((1 + x)/2);

%  Stable cases (zeta > 0):
ah = 10 ; bh = 4.5 ; ch = 3 ; Bh = sqrt(5) ;
%ah = 5 ; bh = 5 ;
y = zeta(zeta > 0) ;
psit(zeta > 0) = -bh/2*log(1 + ch*y + y.*y) + ...
    (-ah/Bh + bh*ch/2/Bh)*(log((2*y+ch-Bh)./(2*y+ch+Bh)) - ...
    log((ch-Bh)/(ch+Bh))) ;

end % End of function psi_t


function gt = dim_ctsq(zeta)

% Computes the dimensionless temperature structure parameter function.
% Input: zeta = z/L, where z is the height and L is the Obukhov length

ct = 5 ;
gt = ct*ones(size(zeta)) ;

% Unstable case (zeta < 0):
gt(zeta < 0) = ct*(1 - 8*zeta(zeta < 0)).^(-2/3);

% Stable case (zeta > 0):
gt(zeta > 0) = ct*(1 + 0.2*zeta(zeta > 0).^(2/3));

end % End of function gt


function [q] = spec_hum(t,rh,p,eps)

% Computes the specific humidity as a function of temperature, relative
% humidity and pressure.
%
% Inputs:
%   t  - Air temperature in C
%   rh - Relative humidity in %
%   pr - Atmospheric pressure in hPa
%
% Output:
%   q  - Specific humidity in kg/kg

%Md = 28.9655;     % Molar mass of dry air with CO2 = 400 ppm (g/mol)
%Mw = 18.015268;   % Molar mass of water vapor (g/mol)
%eps = Mw/Md;      % Ratio of molar mass of water vapor to dry air
% Compute the saturation vapor pressure using the Buck (1981) equations:
es = 6.1121*(1.0007 + 3.46e-6*p) ;
%es = es.*exp(17.502*t./(240.97+t));
if min(t) >= 0
    es = es.*exp(17.368*t./(238.88+t));
else
    es(t >= 0) = es(t >= 0).*exp(17.368*t(t >= 0)./(238.88+t(t >= 0)));
    es(t < 0) = es(t < 0).*exp(17.966*t(t < 0)./(247.15+t(t < 0)));
end
% Compute the vapor pressure:
e = rh.*es/100;
% Compute the specific humidity:
q = eps.*e./(p - (1-eps).*e);

end % End of function spec_hum


function [g] = gravity(lat)

% Computes the acceleration due to gravity as a function of latitude.
%
% Input:  lat - Latitude in degrees
% Output: g   - Acceleration due to gravity in m/s^2

ge = 9.7803253359 ; % Gravity at the equator in m/s^2
x2 = sind(lat).*sind(lat) ;
g = ge*(1+0.00193185265241*x2)./sqrt(1-0.00669437999013*x2) ;
g(abs(lat) > 90) = 9.80665 ;

end % End of function gravity


function [rcurv] = radius_curv(lat,az)
    
% Computes the radius of curvature of the earth in m as a function of 
% latitude and azimuth.    
%
% Inputs: lat   - Latitude in degrees
%         az    - Azimuth of propagation in degrees
% Output: rcurv - Radius of curvature of the earth in m

a = 6378137.0 ; % Earth's equatorial radius (m)
b = 6356752.3 ; % Earth's polar radius (m)
e2 = 1 - b*b/a/a ;
br = 1 - e2*sind(lat).^2 ;
N = a./sqrt(br) ;
M = N.*(1 - e2)./br ;
rcurv = 1./(cosd(az).^2./M + sind(az).^2./N) ;
rcurv(abs(lat) > 90 | az < 0 | az > 360) = 6374371 ; % lat = 40, az = 45

end % End of function radius_curv


function [m1, m2] = Refractivity_wavelength_coeffs(wl, CO2)

% Function to define wavelenth-dependence functions for the refractivity
% of dry air and water vapor for wavelengths between 0.3 and 14 microns, 
% excluding strong H2O and CO2 absorption bands.
%
% Input parameters: 
%   wl  - Wavelength (microns)
%   CO2 - Concentration of atmospheric carbon dioxide (ppm)
% Output parameters: 
%   m1 - Wavelength function for the refractivity of dry air (K/hPa)
%   m2 - Wavelength function for the refractivity of water vapor (K/hPa)
%
%   The refractivity of moist air (N) is given by the equation:
%
%     N = (n - 1) x 10^6 = m1*Ptot/T + (m2 - m1)*e/T
%
% where
%
%     Ptot - Total atmospheric pressure (hPa) = Pdry + e
%     e    - Water vapor partial pressure (hPa)
%     T    - Absolute air temperature (K)
%     n    - Index of refraction of air
%
%   m1 and m2 can be converted to 'specific refractivities' for dry air and
% water vapor by dividing by the ideal gas constants for dry air and water
% vapor, respectively.
%
%   Uses Ciddor (1996), Appl. Opt., vol. 35, 1566-1573 refractivity equations 
% for dry air for 0.3 <= wl <= 2.5 and for water vapor for 0.3 <= wl <= 0.7, 
% and the CO2 adjustment factor equation.  Ignores compressibility of air.  
%   Developed new dry air and water vapor refractivity wavelength-dependence 
% equations for other wavelengths based on polynomial fits to data derived 
% from Colavita et al. (2004), Pub. Astronom. Soc. Pac., vol. 116, 876-885.
%   The CO2 concentration equation is from Picard et al. (2008), Metrologia,
% vol. 45, 149-155.
%
% Developed and written by:
%   Paul A. Frederickson
%   Department of Meteorology
%   Naval Postgraduate School
%   Monterey, CA 93943-5114
%   Phone: (831) 521-8670
%   Email: pafreder@nps.edu
%   April 2015

% Initialize m1 and m2 output arrays and set equal to NaN:
m1 = NaN*ones(size(wl)) ; m2 = NaN*ones(size(wl)) ;

% Define refractivity wavelength-dependence functions for 0.3 to 0.7 microns
i1 = find(wl >= 0.3 & wl <= 0.7) ; 
if ~isempty(i1)
    s = wl(i1).^-2 ;
    m1(i1) = (288.15/101325)*(5792105./(238.0185 - s) + 167917./(57.362 - s)) ;
    m2(i1) = (293.15/1333)*1.022*(295.235 + 2.6422*s - 0.032380*s.^2 + 0.004028*s.^3) ;
    clear s
end

% Define refractivity wavelength-dependence functions for 0.7 to 2.5 microns
i2 = find(wl > 0.7 & wl <= 2.5) ;
if ~isempty(i2)
    s = wl(i2).^-2 ;
    m1(i2) = (288.15/101325)*(5792105./(238.0185 - s) + 167917./(57.362 - s)) ;
    clear s
    s = wl(i2) ;
    p2 = [-1.13859e+00 8.40095e+00 -2.48600e+01 3.70196e+01 -2.89942e+01 7.640266e+01] ;
    m2(i2) = p2(1)*s.^5 + p2(2)*s.^4 + p2(3)*s.^3 + p2(4)*s.^2 + p2(5)*s + p2(6) ;
    clear s p2
end

% Define refractivity wavelength-dependence functions for 2.8 to 4.2 microns
i3 = find(wl >= 2.8 & wl <= 4.2) ;
if ~isempty(i3)
    s = wl(i3) ;
    p1 = [-3.06766e-02 4.56347e-01 -2.69712e+00 7.93057e+00 -1.16668e+01 8.451128e+01] ;
    m1(i3) = p1(1)*s.^5 + p1(2)*s.^4 + p1(3)*s.^3 + p1(4)*s.^2 + p1(5)*s + p1(6) ;
    p2 = [-9.43108e-02 1.83431e+00 -1.47156e+01 6.02712e+01 -1.27375e+02 1.768031e+02] ;
    m2(i3) = p2(1)*s.^5 + p2(2)*s.^4 + p2(3)*s.^3 + p2(4)*s.^2 + p2(5)*s + p2(6) ;
    clear s p1 p2
end

% Define refractivity wavelength-dependence functions for 4.4 to 5.2 microns
i4 = find(wl >= 4.4 & wl <= 5.2) ;
if ~isempty(i4)
    s = wl(i4) ;
    p1 = [-4.80460e-03 1.40384e-01 -1.63588e+00 9.50714e+00 -2.75767e+01 1.095189e+02] ;
    m1(i4) = p1(1)*s.^5 + p1(2)*s.^4 + p1(3)*s.^3 + p1(4)*s.^2 + p1(5)*s + p1(6) ;
    p2 = [-9.43108e-02 1.83431e+00 -1.47156e+01 6.02712e+01 -1.27375e+02 1.768031e+02] ;
    m2(i4) = p2(1)*s.^5 + p2(2)*s.^4 + p2(3)*s.^3 + p2(4)*s.^2 + p2(5)*s + p2(6) ;
    clear s p1 p2
end

% Define refractivity wavelength-dependence functions for 7.5 to 14 microns
i5 = find(wl >= 7.5 & wl <= 14) ;
if ~isempty(i5)
    s = wl(i5) ;
    p1 = [-1.55586e-05 7.63250e-04 -1.48806e-02 1.43954e-01 -6.92399e-01 7.886519e+01] ;
    m1(i5) = p1(1)*s.^5 + p1(2)*s.^4 + p1(3)*s.^3 + p1(4)*s.^2 + p1(5)*s + p1(6) ;
    p2 = [-3.75919e-03 2.11262e-01 -4.74435e+00 5.29907e+01 -2.98185e+02 7.396402e+02] ;
    m2(i5) = p2(1)*s.^5 + p2(2)*s.^4 + p2(3)*s.^3 + p2(4)*s.^2 + p2(5)*s + p2(6) ;
    clear s p1 p2
end

% Correct dry air refractivity for CO2 concentration deviation from 450 ppm:
i6 = find(CO2 > 0 & ~isnan(CO2)) ;
if ~isempty(i6)
    md_CO2 = 28.96546 + 12.0107*(CO2(i6)-400)/10^6 ; % Molar mass of dry air
    md_450 = 28.96546 + 12.0107*(450-400)/10^6 ;
    m1(i6) = m1(i6).*(md_CO2/md_450).*(1 + 0.534*(CO2(i6)-450)/10^6) ;
end

if ~isempty(find(isnan(m1))==1)
    fprintf('Warning! Optical wavelength is bad...\n')
end

end % End of function Refractivity_wavelength_coeffs