close all; clc; clear;

%% Load Data
dataFolder = '.\';    % Data's directory path
outputFolder = '.\output\';    % Output's directory path

warning off;
mkdir(outputFolder);

load ([dataFolder 'cmemsSSThMonthly19932019.mat']);
load ([dataFolder 'era5Wnd10Monthly19932019.mat']);

%% Derived parameters calc
% SST gradient, °C/km
[dxOce,dyOce] = cdtdim(latOce,lonOce,'km'); 
[sstX,sstY] = gradient(ssthetao);
dsstXdx = sstX./dxOce;
dsstYdy = sstY./dyOce;
sstGrad = hypot(dsstXdx,dsstYdy);
% sstGrad(sstGrad.*1e2<2) = NaN;    % Strong thermal front (>0.02 °C/km)

% Wind stress, N/m^2
[TauX,TauY] = windstress(u10,v10,'Cd',0.001,'rho',1.160);
Tau = hypot(TauX,TauY);

% Wind stress curl, N/m^3
[dxAtm,dyAtm] = cdtdim(latAtm,lonAtm,'m');
[dTauXx,dTauXy] = gradient(TauX);
[dTauYx,dTauYy] = gradient(TauY);
TauCurl = (dTauYx./dxAtm) - (dTauXy./dyAtm);

% Crosswind SST gradient
xsstGrad = dsstXdx(1:3:end,1:3:end,:).*TauY-dsstYdy(1:3:end,1:3:end,:).*TauX;

% Ekman Pumping (UE & VE, m^2/s; wE, m/s)
[UE,VE,wE] = ekman(latAtm,lonAtm,u10,v10,'Cd',0.001);

% Upwelling transport, 10^6 m^3/s = 1 Sv
upw = wE.*abs(dxAtm).*abs(dyAtm).*1e-6;
% upw (upw<=0) = NaN;    % Only upwelling (positive value)

%% Climatologies and anomalies calc
variable = [
    "ssthetao";
    "u10"; "v10";
    "sstGrad"; "xsstGrad";
    "TauX"; "TauY"; "Tau"; "TauCurl";
    "UE"; "VE"; "wE"; "upw"
    ];

for nv = 1:length(variable)
    eval(["["+variable(nv)+"Climat,"+variable(nv)+"Anom] = climatAnom ("+variable(nv)+",time(1:length(time)));"]);
end

%% Correlation
corr_sstwnd = []; 
for nlat = 1:length(latAtm(:,1))
    for nlon = 1:length(lonAtm(1,:))
        a = (squeeze(xsstGrad(nlat,nlon,:)));
        b = (squeeze(TauCurl(nlat,nlon,:)));
            
        %Calculate the Correlation
        [r,lag] = xcorr(a,b,'coeff');
        lag=permute(lag,[2 1]);
            
        %Find Correlation on Zero Lag
        c=find(lag(:,1)==0);
        corr_sstwnd(nlat,nlon)=(r(c));
    end
end

%% Visualization

n = 7;    % North
s = -13;      % South
w = 94;   % West
e = 142;    % East

% % Single time plot (annual mean)
figure ('Color','w');   % Set figure's window
m_proj('miller','lon',[w e],'lat',[s n]); % Map's projection

% Shaded plot
m_contourf(lonOce,latOce,nanmean(ssthetaoClimat,3),(26:0.25:32),'edgeColor','none'); hold on;   % Filled 2-D contour plot

% Contour plot
m_contour(lonOce,latOce,nanmean(sstGradClimat,3).*1e2,(0:0.75:0.75),'lineColor','k','lineStyle','-','lineWidth',1); hold on;

% Quiver or velocity plot

m_gshhs_i('patch',[.8 .8 .8],'edgeColor',[.2 .2 .2],'linewidth',.5); hold on;  % Map's land patch 
m_grid('linestyle','none','xtick',6,'ytick',4,'fontsize',16);  % Map's grid line

% title(["Mean Sea Surface Temperature"],'fontsize',16,'fontweight','b');   % Map's title

xlabel ('Longitude','fontsize',16); % Label x-axis
ylabel ('Latitude','fontsize',16);  % Label y-axis
cb = colorbar ('eastoutside','fontsize',16);    % Colorbar
ylabel(cb,'Annual mean of sea surface temperature (°C)','fontsize',16) % Colorbar's label
caxis([26 32]);  % Colorbar's limit
cmap = getPyPlot_cMap('RdYlBu_r', 24);    % Colormap
colormap (cmap);

set(gcf, 'Position', get(0, 'Screensize')); % Set figure's window fullscreen

export_fig([outputFolder ['meanSST']], '-png');   % Export figure

figure ('Color','w');   % Set figure's window
m_proj('miller','lon',[w e],'lat',[s n]); % Map's projection

% Shaded plot
% m_pcolor(lonOce,latOce,nanmean(ssthetaoClimat,3));    % Pseudocolor (checkerboard) plot
% shading flat; hold on;
m_contourf(lonAtm,latAtm,nanmean(wEClimat,3).*1e4,(-1:1/12:1),'edgeColor','none'); hold on;   % Filled 2-D contour plot

% Contour plot
% m_contour(lonAtm,latAtm,nanmean(wEClimat,3).*1e4,(0:1/12:1),'lineColor','k','lineStyle','-','lineWidth',1); hold on;

% Quiver or velocity plot
m_quiver(lonAtm(1:4:end,1:4:end),latAtm(1:4:end,1:4:end),nanmean(TauXClimat(1:4:end,1:4:end,:),3),nanmean(TauYClimat(1:4:end,1:4:end,:),3),'Color','k','lineWidth',1);
hold on;
m_quiver(lonAtm(1:4:end,1:4:end),latAtm(1:4:end,1:4:end),nanmean(UEClimat(1:4:end,1:4:end,:),3),nanmean(VEClimat(1:4:end,1:4:end,:),3),'Color','r','lineWidth',1);
hold on;

m_gshhs_i('patch',[.8 .8 .8],'edgeColor',[.2 .2 .2],'linewidth',.5); hold on;  % Map's land patch 
m_grid('linestyle','none','xtick',6,'ytick',4,'fontsize',16);  % Map's grid line

% title(["Mean Sea Surface Temperature"],'fontsize',16,'fontweight','b');   % Map's title

xlabel ('Longitude','fontsize',16); % Label x-axis
ylabel ('Latitude','fontsize',16);  % Label y-axis
cb = colorbar ('eastoutside','fontsize',16);    % Colorbar
ylabel(cb,'Annual mean of Ekman pumping ({\times}10^{-4} m/s)','fontsize',16) % Colorbar's label
caxis([-1 1]);  % Colorbar's limit
cmap = getPyPlot_cMap('BrBG', 24);    % Colormap
colormap (cmap);

set(gcf, 'Position', get(0, 'Screensize')); % Set figure's window fullscreen

export_fig([outputFolder ['meanwE']], '-png');   % Export figure

% % Multiple time plot (climatologies)
for nt = 1:12
    figure ('Color','w')
    m_proj('miller','lon',[w e],'lat',[s n]);

%     m_pcolor(lonOce,latOce,sstClimat(:,:,nt));
%     shading flat; hold on;
%     m_contourf(lonOce,latOce,sstClimat(:,:,nt),(24:0.25:32),'edgeColor','none'); hold on;
    
%     m_pcolor(lonAtm,latAtm,wEClimat(:,:,nt).*1e4);
%     shading flat; hold on;
    m_contourf(lonAtm,latAtm,wEClimat(:,:,nt).*1e4,(-2:2/12:2),'edgeColor','none'); hold on;

    m_contour(lonOce,latOce,sstGradClimat(:,:,nt).*1e2,(0:2:2),'lineColor','k','lineStyle','-','lineWidth',2); hold on;
    
%     m_quiver(lonOce(1:6:end,1:6:end),latOce(1:6:end,1:6:end),uoMonthly(1:6:end,1:6:end,nt),voMonthly(1:6:end,1:6:end,nt),'Color',[.4 .4 .4],'lineWidth',1.5);
%     hold on;
    m_quiver(lonAtm(1:4:end,1:4:end),latAtm(1:4:end,1:4:end),TauXClimat(1:4:end,1:4:end,nt),TauYClimat(1:4:end,1:4:end,nt),'Color',[.4 .4 .4],'lineWidth',1.5);
    hold on;

    m_gshhs_i('patch',[.8 .8 .8],'edgeColor',[.2 .2 .2],'linewidth',.5); hold on;   
    m_grid('linestyle','none','xtick',6,'ytick',4,'fontsize',16);

    title([datestr(datenum(1,nt,1),'mmmm')],'fontsize',16,'fontweight','b');

    xlabel ('Longitude','fontsize',16);
    ylabel ('Latitude','fontsize',16);
    cb = colorbar ('eastoutside','fontsize',16);
    xlabel(cb,'Ekman pumping ({\times}10^{-4} m/s)','fontsize',16)

    caxis([-2 2]);
    cmap = getPyPlot_cMap('RdBu', 24);
    colormap (cmap);
    
    set(gcf, 'Position', get(0, 'Screensize'));

    export_fig([outputFolder sprintf(['wE' datestr(datenum(time(nt,:)),'mm')])], '-png');

    close;
end

% % Multiple time plot (seasonal)
season = [
    "DJF" "MAM" "JJA" "SON";
    "Dec - Jan - Feb" "Mar - Apr - May" "Jun - Jul - Aug" "Sep - Oct - Nov";
    "[1:2,12]" "[3:5]" "[6:8]" "[9:11]";
    ];

for nt = 1:4
    figure ()
    m_proj('miller','lon',[w e],'lat',[s n]);
    
    % Shaded plot
    eval(["m_contourf(lonAtm,latAtm,nanmean(wEClimat(:,:,"+season(3,nt)+"),3).*1e4,(-2:2/12:2),'edgeColor','none');"]);
    hold on;   % Filled 2-D contour plot

    % Contour plot
    eval(["m_contour(lonOce,latOce,nanmean(sstGradClimat(:,:,"+season(3,nt)+"),3).*1e2,(0:2:2),'lineColor','k','lineStyle','-','lineWidth',2);"]);
    hold on;   % Filled 2-D contour plot

    % Quiver or velocity plot
    eval(["m_quiver(lonAtm(1:4:end,1:4:end),latAtm(1:4:end,1:4:end),nanmean(TauXClimat(1:4:end,1:4:end,"+season(3,nt)+"),3),nanmean(TauYClimat(1:4:end,1:4:end,"+season(3,nt)+"),3),'Color',[.4 .4 .4],'lineWidth',1.5);"]);
    hold on;
    
    m_gshhs_i('patch',[.8 .8 .8],'edgeColor',[.2 .2 .2],'linewidth',.5); hold on;  % Map's land patch 
    m_grid('linestyle','none','xtick',6,'ytick',4,'fontsize',16);  % Map's grid line

    title([season(2,nt)],'fontsize',16,'fontweight','b');   % Map's title

    xlabel ('Longitude','fontsize',16); % Label x-axis
    ylabel ('Latitude','fontsize',16);  % Label y-axis
    cb = colorbar ('eastoutside','fontsize',16);    % Colorbar
    xlabel(cb,'Ekman pumping ({\times}10^{-4} m/s)','fontsize',16); % Colorbar's label
    caxis([-2 2]);  % Colorbar's limit
    cmap = getPyPlot_cMap('RdBu_r', 24);    % Colormap
    colormap (cmap);

    set(gcf, 'Position', get(0, 'Screensize')); % Set figure's window fullscreen

    export_fig([outputFolder ['meanwE' char(season(1,nt))]], '-png', '-transparent');   % Export figure
end

% % Spatial correlation plot
figure ('Color','w');   % Set figure's window
m_proj('miller','lon',[w e],'lat',[s n]); % Map's projection

% Shaded plot
m_pcolor(lonAtm,latAtm,corr_sstwnd);    % Pseudocolor (checkerboard) plot
shading flat; hold on;
m_contourf(lonAtm,latAtm,corr_sstwnd,(-1:0.1:1),'edgeColor','none'); hold on;   % Filled 2-D contour plot

% Contour plot
m_contour(lonAtm,latAtm,corr_sstwnd,-1:1.2:0.2,'lineColor','k','lineStyle','--','lineWidth',1); hold on;

m_gshhs_i('patch',[.8 .8 .8],'edgeColor',[.2 .2 .2],'linewidth',.5); hold on;  % Map's land patch 
m_grid('linestyle','none','xtick',6,'ytick',4,'fontsize',16);  % Map's grid line

% title(["Spatial Correlation of Crosswind-SST Gradient and Wind's Curl"],'fontsize',16,'fontweight','b');   % Map's title

xlabel ('Longitude','fontsize',16); % Label x-axis
ylabel ('Latitude','fontsize',16);  % Label y-axis
cb = colorbar ('eastoutside','fontsize',16);    % Colorbar
ylabel(cb,'Correlation coefficient','fontsize',16) % Colorbar's label
caxis([-1 1]);  % Colorbar's limit
cmap = getPyPlot_cMap('RdBu_r', 20);    % Colormap
colormap (cmap);

set(gcf, 'Position', get(0, 'Screensize')); % Set figure's window fullscreen

export_fig([outputFolder ['corrcoef2']], '-png', '-transparent');   % Export figure

%% Subset
[~,latOceMin] = min(abs(latOce(:,1)-n));
[~,latOceMax] = min(abs(latOce(:,1)-s));
[~,lonOceMin] = min(abs(lonOce(1,:)-w));
[~,lonOceMax] = min(abs(lonOce(1,:)-e));

[~,latAtmMin] = min(abs(latAtm(:,1)-n));
[~,latAtmMax] = min(abs(latAtm(:,1)-s));
[~,lonAtmMin] = min(abs(lonAtm(1,:)-w));
[~,lonAtmMax] = min(abs(lonAtm(1,:)-e));

sstAnomSS = squeeze(nanmean(nanmean(sstAnom(latOceMin:latOceMax,lonOceMin:lonOceMax,:),1),2));
% mask = nanvar(sstGradMonthly(latMin:latMax,lonMin:lonMax,:),[],3).*1e6;
% mask (mask<20) = NaN;
% mask (mask>=20) = 1;

Subset = [
    "maskJava" "JS";
    "maskNorthJava" "NJ";
    "maskEastSumatra" "ES";
    "maskWestBorneo" "WB";
    "maskEastBorneo" "WE";
    ];

for nv = [1:3,6]
    eval([""+variable(nv)+"SS = [];"]);
    eval([""+variable(nv)+"MonthlySS = [];"]);
    eval([""+variable(nv)+"AnomSS = [];"]);
        
    for ns = 1:5
        eval(["load "+Subset(ns,1)+".txt;"]);
        eval([""+variable(nv)+""+Subset(ns,2)+" = squeeze(nanmean(nanmean("+variable(nv)+"(latOceMin:latOceMax,lonOceMin:lonOceMax,:).*"+Subset(ns,1)+",1),2));"]); 
        eval([""+variable(nv)+"SS = ["+variable(nv)+"SS "+variable(nv)+""+Subset(ns,2)+"];"]);
        eval([""+variable(nv)+"Monthly"+Subset(ns,2)+" = squeeze(nanmean(nanmean("+variable(nv)+"Monthly(latOceMin:latOceMax,lonOceMin:lonOceMax,:).*"+Subset(ns,1)+",1),2));"]); 
        eval([""+variable(nv)+"MonthlySS = ["+variable(nv)+"MonthlySS "+variable(nv)+"Monthly"+Subset(ns,2)+"];"]);
        eval([""+variable(nv)+"Anom"+Subset(ns,2)+" = squeeze(nanmean(nanmean("+variable(nv)+"Anom(latOceMin:latOceMax,lonOceMin:lonOceMax,:).*"+Subset(ns,1)+",1),2));"]); 
        eval([""+variable(nv)+"AnomSS = ["+variable(nv)+"AnomSS "+variable(nv)+"Anom"+Subset(ns,2)+"];"]);
    end
end

for nv = [4:5,7:15]
    eval([""+variable(nv)+"SS = [];"]);
    eval([""+variable(nv)+"MonthlySS = [];"]);
    eval([""+variable(nv)+"AnomSS = [];"]);
        
    for ns = 1:5
        eval(["load "+Subset(ns,1)+".txt;"]);
        eval([""+variable(nv)+""+Subset(ns,2)+" = squeeze(nanmean(nanmean("+variable(nv)+"(latAtmMin:latAtmMax,lonAtmMin:lonAtmMax,:).*"+Subset(ns,1)+"(1:3:end,1:3:end),1),2));"]); 
        eval([""+variable(nv)+"SS = ["+variable(nv)+"SS "+variable(nv)+""+Subset(ns,2)+"];"]);
        eval([""+variable(nv)+"Monthly"+Subset(ns,2)+" = squeeze(nanmean(nanmean("+variable(nv)+"Monthly(latAtmMin:latAtmMax,lonAtmMin:lonAtmMax,:).*"+Subset(ns,1)+"(1:3:end,1:3:end),1),2));"]); 
        eval([""+variable(nv)+"MonthlySS = ["+variable(nv)+"MonthlySS "+variable(nv)+"Monthly"+Subset(ns,2)+"];"]);
        eval([""+variable(nv)+"Anom"+Subset(ns,2)+" = squeeze(nanmean(nanmean("+variable(nv)+"Anom(latAtmMin:latAtmMax,lonAtmMin:lonAtmMax,:).*"+Subset(ns,1)+"(1:3:end,1:3:end),1),2));"]); 
        eval([""+variable(nv)+"AnomSS = ["+variable(nv)+"AnomSS "+variable(nv)+"Anom"+Subset(ns,2)+"];"]);
    end
end
