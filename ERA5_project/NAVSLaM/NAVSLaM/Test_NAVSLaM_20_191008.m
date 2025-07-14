function Test_NAVSLaM_20_191008

% Function to run the MATLAB NAVSLaM Version 2.0 model code with
% test case inputs.
%
% Written by Paul Frederickson, Naval Postgraduate School, August 2018

% Loop to read the four sample input data files and run NAVSLaM:
for ifile = 1 : 4

    % Open input and output data files:
    if ifile == 1 % Radio frequency & relative humidity cases:
        fin = fopen('NAVSLaM_20_test_input_rf_rh.txt','r') ;
        fout = fopen('NAVSLaM_20_test_output_rf_rh.txt','wt') ;
        fprintf('Radio frequency case with relative humidity inputs:\n')
        tt = '(RF, RH inputs)' ;
        hflag = 1 ;
    elseif ifile == 2 % Radio frequency & specific humidity cases:
        fin = fopen('NAVSLaM_20_test_input_rf_q.txt','r') ;
        fout = fopen('NAVSLaM_20_test_output_rf_q.txt','wt') ;
        fprintf('Radio frequency case with specific humidity inputs:\n')
        tt = '(RF, q inputs)' ;
        hflag = 2 ;
    elseif ifile == 3 % Optical wavelength & relative humidity cases:
        fin = fopen('NAVSLaM_20_test_input_op_rh.txt','r') ;
        fout = fopen('NAVSLaM_20_test_output_op_rh.txt','wt') ;
        fprintf('Optical wavelength case with relative humidity inputs:\n')
        tt = '(Optical, RH inputs)' ;
        hflag = 1 ;
    elseif ifile == 4 % Optical wavelength & specific humidity cases:
        fin = fopen('NAVSLaM_20_test_input_op_q.txt','r') ;
        fout = fopen('NAVSLaM_20_test_output_op_q.txt','wt') ;
        fprintf('Optical wavelength case with specific humidity inputs:\n')
        tt = '(Optical, q inputs)' ;
        hflag = 2 ;
    end

    % Read input data and define input data parameters:
    C = textscan(fin,'%d %n %n %n %n %n %n %n %d %n %n %n %n %n %n %n %n'); 
    fclose(fin) ;
    tcase(:,1) = C{1} ;
    lambda(:,1) = C{2} ;
    zmax(:,1) = C{3} ;
    zinc(:,1) = C{4} ;
    ws(:,1) = C{5} ;
    tair(:,1) = C{6} ;
    tsea(:,1) = C{7} ;
    h(:,1) = C{8} ;
    hind(:,1) = C{9} ;
    pr(:,1) = C{10} ;
    s(:,1) = C{11} ;
    zu(:,1) = C{12} ;
    zt(:,1) = C{13} ;
    zq(:,1) = C{14} ;
    zp(:,1) = C{15} ; 
    lat(:,1) = C{16} ;
    az(:,1) = C{17} ;
    clear C
    
%    if length(hind(hind == 1)) == length(hind) ; hflag = 1 ; end
%    if length(hind(hind == 2)) == length(hind) ; hflag = 2 ; end

    % Run NAVSLaM Version 2.0 with input data arrays:
    [ustar,tstar,qstar,l,edh,mmin,z,u,t,q,p,m,ctsq,cqsq,ctq,cnsq] = ...
        NAVSLaM_20_191008(lambda,ws,tair,tsea,h,hflag,pr,s,lat,az,...
        zu,zt,zq,zp,zinc(1),zmax(1)) ;
    
    iz = find(z==5) ;
    q(q>0) = q(q>0)*1000 ;
    fmt = ['%2d %5.1f %8.4f %12.5e %12.5e %12.5e %12.5e %8.4f %8.4f %8.4f ' ...
        '%7.2f %8.4f %12.5e\n'] ;
    
    % Loop to write output data to files:
    for i = 1 : length(tcase)
        fprintf(fmt, ...
            tcase(i),edh(i),mmin(i),ustar(i),tstar(i),qstar(i),l(i), ...
            u(i,iz),t(i,iz),q(i,iz),p(i,iz),m(i,iz),cnsq(i,iz)) ;
        fprintf(fout,fmt, ...
            tcase(i),edh(i),mmin(i),ustar(i),tstar(i),qstar(i),l(i), ...
            u(i,iz),t(i,iz),q(i,iz),p(i,iz),m(i,iz),cnsq(i,iz)) ;
    end
    fclose(fout) ;
    
    % Set bad outputs to NaN for plotting purposes:
    edh(edh == -96) = zmax(1) ;
    u(u<0) = NaN ; t(t<=-97) = NaN ; q(q<=-97) = NaN ; p(p<=-97) = NaN ;
    m(m<=-97) = NaN ; ctsq(ctsq<0) = NaN ; cqsq(cqsq<0) = NaN ; 
    ctq(ctq<0) = NaN ; cnsq(cnsq<0) = NaN ;
    
    % Plot vertical profiles for test case i:
    i = 2 ;
    cn2 = squeeze(cnsq(i,:)) ; cn2(cn2<=0) = NaN ;
    
%    ta = squeeze(t(i,:)) ; ta(ta<=-97) = NaN ;
%    qa = squeeze(q(i,:)) ; qa(qa<=-97) = NaN ;
%    pres = squeeze(p(i,:)) ; pres(pres<=-97) = NaN ;
%    ref = squeeze(m(i,:)) ; ref(ref<=-97) = NaN ;

    % Create figure:
    figure(ifile)
    set(gcf,'position',[10 50 1800 1000])
    fsize = 13 ;
    xpos = [0.03 0.275 0.520 0.765] ; ypos = [0.07 0.56] ; 
    xlen = 0.22 ; ylen = 0.41 ;
    
    % Temperature
    axes('position',[xpos(1) ypos(2) xlen ylen])
    plot(t(i,:),z,'b-','linewidth',2)
    hold on
    plot(tair(i),zt(i),'bo','markerfacecolor','b')
    plot(tsea(i),0,'bo','markerfacecolor','b')
    xlabel('Temperature (\circC)','fontweight','bold','fontsize',fsize)
    ylabel('Height (m)','fontweight','bold','fontsize',fsize)
    title(['File ' int2str(ifile) ' ' tt ', Test Case ' int2str(i)],...
        'fontsize',fsize+1)
    set(gca,'linewidth',2,'ylim',[0 zmax(i)],'fontweight','bold')
    
    % Specific Humidity
    axes('position',[xpos(2) ypos(2) xlen ylen])
    plot(q(i,:),z,'b-','linewidth',2)
    if hflag == 2
        hold on
        plot(h(i)*1000,zq(i),'bo','markerfacecolor','b')
    end
    xlabel('Specific Humidity (g/kg)','fontweight','bold','fontsize',fsize)
    set(gca,'linewidth',2,'ylim',[0 zmax(i)],'fontweight','bold')

    % Pressure
    axes('position',[xpos(3) ypos(2) xlen ylen])
    plot(p(i,:),z,'b-','linewidth',2)
    hold on
    plot(pr(i),zp(i),'bo','markerfacecolor','b')
    xlabel('Pressure (hPa)','fontweight','bold','fontsize',fsize)
    set(gca,'linewidth',2,'ylim',[0 zmax(i)],'fontweight','bold')

    % Wind speed
    axes('position',[xpos(4) ypos(2) xlen ylen])
    plot(u(i,:),z,'b-','linewidth',2)
    hold on
    plot(ws(i),zu(i),'bo','markerfacecolor','b')
    xlabel('Wind Speed (m/s)','fontweight','bold','fontsize',fsize)
    set(gca,'linewidth',2,'ylim',[0 zmax(i)],'fontweight','bold')
    
    % Refractivity
    axes('position',[xpos(1) ypos(1) xlen ylen])
    plot(m(i,:),z,'b-','linewidth',2)
    hold on 
    if edh(i) >= 0
       yp = plot(m(i,z==edh(i)),edh(i),'bo','markerfacecolor','b') ;
    end
    if ifile == 1 || ifile == 2
    xlabel('Modified Refractivity','fontweight','bold','fontsize',fsize)
    legend(yp,['EDH = ' num2str(edh(i),'%5.1f') ' m'])
    else
    xlabel('Refractivity','fontweight','bold','fontsize',fsize)
    end
    ylabel('Height (m)','fontweight','bold','fontsize',fsize)
    set(gca,'linewidth',2,'ylim',[0 zmax(i)],'fontweight','bold')
    
    % CT2
    axes('position',[xpos(2) ypos(1) xlen ylen])
    semilogx(ctsq(i,:),z,'b-','linewidth',2)
    xlabel('C_T^2 (K^2 m^-^2^/^3)','fontweight','bold','fontsize',fsize)
    set(gca,'linewidth',2,'ylim',[0 zmax(i)],'fontweight','bold')
    
    % Cq2
    axes('position',[xpos(3) ypos(1) xlen ylen])
    semilogx(cqsq(i,:),z,'b-','linewidth',2)
    xlabel('C_q^2 (m^-^2^/^3)','fontweight','bold','fontsize',fsize)
    set(gca,'linewidth',2,'ylim',[0 zmax(i)],'fontweight','bold')
    
    % Cn2
    axes('position',[xpos(4) ypos(1) xlen ylen])
    semilogx(cnsq(i,:),z,'b-','linewidth',2)
    xlabel('C_n^2 (m^-^2^/^3)','fontweight','bold','fontsize',fsize)
    set(gca,'linewidth',2,'ylim',[0 zmax(i)],'fontweight','bold')

    set(gcf,'paperpositionmode','auto')
    print(['NAVSLaM_20_Test_File' int2str(ifile) '_Case' int2str(i) ...
        '.png'],'-dpng','-r300')
    
    clear tcase lambda zinc zmax ws tair tsea h hflag hind pr s
    clear zu zt zq zp lat az ustar tstar qstar l edh z u t q p m cnsq
    clear ctsq ctq cqsq
    
end % end of ifile loop


% Example of running NAVSLaM with a pre-defined height array with
% logarithmic vertical spacing for optical wavelength cases

% Pre-define height array:
z(1) = 0 ; 
z(2:12) = (10.^[0:0.1:1])/10 ;
z(13:32) = 10.^[0.1:0.1:2] ;
%z(13:52,1) = 10.^[0.05:0.05:2] ; % For higher vertical resolution
zinc = z ; nz = length(z) ;
zmax = NaN ; % Not used when zinc is a pre-defined height array

% Define input data:
lambda = 1.06 ;
u = [7 5 10]' ;
tair = [18 19 20]' ;
tsea = [20 19 19]' ;
h = [75 70 80]' ; hflag = 1 ;
pr = [1013 1010 1020]' ;
s = [34 35 33]' ;
lat = [38 45 30]' ;
az = [45 90 215]' ;
zu = 10 ;
zt = 2 ; zq = zt ;
zp = 0 ;

% Run NAVSLaM Version 2.0 with input data array:
[ustar,tstar,qstar,l,edh,mmin,z,u,t,q,p,m,ctsq,cqsq,ctq,cnsq] = ...
    NAVSLaM_20_191008(lambda,u,tair,tsea,h,hflag,pr,s,lat,az,...
    zu,zt,zq,zp,zinc,zmax) ;

% Plot Cn2 profiles for all three test cases:
figure(5)
set(gcf,'position',[10 50 1000 800])
fsize = 13 ;
% Refractivity
axes('position',[0.06 0.09 0.9 0.87])
col = {'r','g','b'} ;
for i = 1 : 3
    loglog(cnsq(i,2:nz),z(2:nz),'o','markersize',5,'color',col{i},...
        'markerfacecolor',col{i})
    hold on
    yp(i) = loglog(cnsq(i,2:nz),z(2:nz),'-','linewidth',2,'color',col{i}) ;
end
xlabel('C_n^2 (m^-^2^/^3)','fontweight','bold','fontsize',fsize)
ylabel('Height (m)','fontweight','bold','fontsize',fsize)
title('NAVSLaM Examples with Logarithmic Vertical Spacing', ...
    'fontsize',fsize+2)
legend(yp,{'Case 1','Case 2','Case 3'})
set(gca,'linewidth',2,'fontweight','bold')

set(gcf,'paperpositionmode','auto')
print('NAVSLaM_20_logz_example.png','-dpng','-r300')

end % End of function Test_NAVSLaM_20_191008