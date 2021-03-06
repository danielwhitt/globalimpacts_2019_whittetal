% heat map
% need to run fig 1 fric vel and 7 and save end workspace
% run fig 7 first
restoredefaultpath
clear all
close all
addpath ./utility
addpath ./cmocean_v1
addpath ./export_fig-master/

%load('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/MLD_data.mat','diffmo','varoutoldclmo','epsilon1');
load('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/N2barXMXLp2monmean.mat','N2barXMXLp2monmean016','N2barXMXLp2monmean036');
load('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholsonCarranza_public/fig_4_bflux_CI.mat','varoutoldmo');
Bfluxoldclmo=squeeze(nanmean(varoutoldmo,4));
load('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholsonCarranza_public/fig_4_bflux_CI.mat','varoutmo');
Bfluxclmo=squeeze(nanmean(varoutmo,4));
clear varoutmo
load('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholsonCarranza_public/fig_1_fricvel_fig_9_load.mat','varoutoldclmo','tlat','tlon');
ustaroldclmo=varoutoldclmo;
load('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholsonCarranza_public/fig_1_fricvel_fig_9_load.mat','varoutclmo','tlat','tlon');
ustarclmo=varoutclmo;

load('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholsonCarranza_public/fig_1_fricvel_fig_9_load.mat','diffmo','varoutoldclmo','epsilon1');
clear tempvar tempmeanvar
tempvar=diffmo;
tempvar=reshape(tempvar,[size(tempvar,1)*size(tempvar,2)*size(tempvar,3) size(tempvar,4)]);
tempmeanvar=remove_insignificant_points(tempvar,0);
diffmomeanyrs=reshape(tempmeanvar,[size(diffmo,1) size(diffmo,2) size(diffmo,3)]);
deltaustar=diffmomeanyrs./(abs(varoutoldclmo));
clear diffmo varoutoldclmo
clear varoutoldclmo varoutclmo
clear varoutoldmo
load('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholsonCarranza_public/fig_7_XMXL_fig_9_load.mat','diffmo','varoutoldclmo','varoutclmo','epsilon1');
N2barXMXLp2monmean016=N2barXMXLp2monmean016(1:6:end,1:6:end,:);
N2barXMXLp2monmean036=N2barXMXLp2monmean036(1:6:end,1:6:end,:);

for flag = 5:6
coords=load('/glade/work/dwhitt/roms.1dmodel/coords_boxes.txt');
coords(abs(coords)>500)=nan;
CASENM=strcat('FILE4',num2str(flag))
POPHISTPATH='/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/archive/g.e20b07.2000_DATM%NYF_SLND_CICE_POP2_DROF%NYF_SGLC_SWAV.T62_t13.hybrid.016/ocn/hist/'
POPHISTFILE='g.e20b07.2000_DATM%NYF_SLND_CICE_POP2_DROF%NYF_SGLC_SWAV.T62_t13.hybrid.016.pop.h.0011-02.nc'
fn=strcat(POPHISTPATH,POPHISTFILE);
thresh=0.3;
lonmin=coords(flag,1)
lonmax=coords(flag,3)+coords(flag,1)
latmin=coords(flag,2)
latmax=coords(flag,2)+coords(flag,4)
addpath /glade/work/dwhitt/roms.1dmodel
[stidx{flag},lenidx{flag}]=identify_pop_idx(fn,thresh,lonmin,lonmax,latmin,latmax)
clear lonmax lonmin latmax latmin thresh fn coords
end

coriolis = 2.*7.292E-5.*sind(tlat);

figure;
for it = 1:4
subplot(2,2,it),...
load coast
axm=axesm('MapProjection','robinson');
surfm(tlat,tlon,real(log10(sqrt(N2barXMXLp2monmean016(:,:,3*it-2))./abs(coriolis)))); % this makes a color-filled plot
patchesm(lat,long,0*[1 1 1]);    % this makes grey shading over land
framem;
gridm('on')
gridm(':')
gridm('GLineWidth',1)
axis off
h = colorbar('h','location','eastoutside')
title(['P, log_{10} N/f, Month: ',num2str(3*it-2)]);
caxis([1 3.5])
addpath ~/cmocean_v1
%cmap=cmocean('thermal');
colormap(gca,jet)
set(gcf,'color','w')
end
export_fig('fig_Pr_4panel.png','-r300')

figure;
for it = 1:4
subplot(2,2,it),...
load coast
axm=axesm('MapProjection','robinson');
surfm(tlat,tlon,real((Bfluxoldclmo(:,:,3*it-2))./(abs(coriolis).*(ustaroldclmo(:,:,3*it-2).^2)))); % this makes a color-filled plot
patchesm(lat,long,0*[1 1 1]);    % this makes grey shading over land
framem;
gridm('on')
gridm(':')
gridm('GLineWidth',1)
axis off
h = colorbar('h','location','eastoutside')
title(['Z=B/u_*^2f, Month: ',num2str(3*it-2)]);
caxis([-30 30])
addpath ~/cmocean_v1
cmap=cmocean('balance');
colormap(gca,cmap)
set(gcf,'color','w')
end
export_fig('fig_Z_4panel.png','-r300')
Zmoold=(Bfluxoldclmo)./(repmat(abs(coriolis),[1 1 12]).*(ustaroldclmo.^2));
Pmoold=real(log10(sqrt(N2barXMXLp2monmean016)./repmat(abs(coriolis),[1 1 12])));
Zmo=Zmoold;
Pmo=Pmoold;
ZmoLP=(Bfluxclmo)./(repmat(abs(coriolis),[1 1 12]).*(ustarclmo.^2));
PmoLP=real(log10(sqrt(N2barXMXLp2monmean036)./repmat(abs(coriolis),[1 1 12])));


clear diffmomeanyrs diffmostdyrs y3
clear tempvar tempmeanvar
tempvar=diffmo;
tempvar=reshape(tempvar,[size(tempvar,1)*size(tempvar,2)*size(tempvar,3) size(tempvar,4)]);
tempmeanvar=remove_insignificant_points(tempvar,0);
diffmomeanyrs=reshape(tempmeanvar,[size(diffmo,1) size(diffmo,2) size(diffmo,3)]);
deltaMLD=diffmomeanyrs./(abs(varoutoldclmo));


%%% PZ circles
figure;
subplot(4,2,1),...
scatter(squeeze(nanmean(nanmean(Zmoold(round(stidx{5}(1)/6):(round((stidx{5}(1)+lenidx{5}(1))/6)),round(stidx{5}(2)/6):(round((stidx{5}(2)+lenidx{5}(2))/6)),:),1),2)),...
    squeeze(nanmean(nanmean(Pmoold(round(stidx{5}(1)/6):(round((stidx{5}(1)+lenidx{5}(1))/6)),round(stidx{5}(2)/6):(round((stidx{5}(2)+lenidx{5}(2))/6)),:),1),2)),...
    squeeze(nanmean(nanmean(varoutoldclmo(round(stidx{5}(1)/6):(round((stidx{5}(1)+lenidx{5}(1))/6)),round(stidx{5}(2)/6):(round((stidx{5}(2)+lenidx{5}(2))/6)),:),1),2)),1:12,'filled');
hold on;
scatter(squeeze(nanmean(nanmean(ZmoLP(round(stidx{5}(1)/6):(round((stidx{5}(1)+lenidx{5}(1))/6)),round(stidx{5}(2)/6):(round((stidx{5}(2)+lenidx{5}(2))/6)),:),1),2)),...
    squeeze(nanmean(nanmean(PmoLP(round(stidx{5}(1)/6):(round((stidx{5}(1)+lenidx{5}(1))/6)),round(stidx{5}(2)/6):(round((stidx{5}(2)+lenidx{5}(2))/6)),:),1),2)),...
    squeeze(nanmean(nanmean(varoutclmo(round(stidx{5}(1)/6):(round((stidx{5}(1)+lenidx{5}(1))/6)),round(stidx{5}(2)/6):(round((stidx{5}(2)+lenidx{5}(2))/6)),:),1),2)),1:12,'filled','marker','square');
grid on
colormap(gca,jet);

subplot(4,2,2),...
scatter(squeeze(nanmean(nanmean(Zmoold(round(stidx{6}(1)/6):(round((stidx{6}(1)+lenidx{6}(1))/6)),round(stidx{6}(2)/6):(round((stidx{6}(2)+lenidx{6}(2))/6)),:),1),2)),...
    squeeze(nanmean(nanmean(Pmoold(round(stidx{6}(1)/6):(round((stidx{6}(1)+lenidx{6}(1))/6)),round(stidx{6}(2)/6):(round((stidx{6}(2)+lenidx{6}(2))/6)),:),1),2)),...
    squeeze(nanmean(nanmean(varoutoldclmo(round(stidx{6}(1)/6):(round((stidx{6}(1)+lenidx{6}(1))/6)),round(stidx{6}(2)/6):(round((stidx{6}(2)+lenidx{6}(2))/6)),:),1),2)),1:12,'filled');
hold on;
scatter(squeeze(nanmean(nanmean(ZmoLP(round(stidx{6}(1)/6):(round((stidx{6}(1)+lenidx{6}(1))/6)),round(stidx{6}(2)/6):(round((stidx{6}(2)+lenidx{6}(2))/6)),:),1),2)),...
    squeeze(nanmean(nanmean(PmoLP(round(stidx{6}(1)/6):(round((stidx{6}(1)+lenidx{6}(1))/6)),round(stidx{6}(2)/6):(round((stidx{6}(2)+lenidx{6}(2))/6)),:),1),2)),...
    squeeze(nanmean(nanmean(varoutclmo(round(stidx{6}(1)/6):(round((stidx{6}(1)+lenidx{6}(1))/6)),round(stidx{6}(2)/6):(round((stidx{6}(2)+lenidx{6}(2))/6)),:),1),2)),1:12,'filled','marker','square');
grid on
colormap(gca,jet);

deltaOUT=(deltaMLD-deltaustar);
%deltaOUT=deltaMLD;
subplot(4,2,3),...
scatter(squeeze(nanmean(nanmean(Zmoold(round(stidx{5}(1)/6):(round((stidx{5}(1)+lenidx{5}(1))/6)),round(stidx{5}(2)/6):(round((stidx{5}(2)+lenidx{5}(2))/6)),:),1),2)),...
    squeeze(nanmean(nanmean(Pmoold(round(stidx{5}(1)/6):(round((stidx{5}(1)+lenidx{5}(1))/6)),round(stidx{5}(2)/6):(round((stidx{5}(2)+lenidx{5}(2))/6)),:),1),2)),60,...
    squeeze(nanmean(nanmean(deltaOUT(round(stidx{5}(1)/6):(round((stidx{5}(1)+lenidx{5}(1))/6)),round(stidx{5}(2)/6):(round((stidx{5}(2)+lenidx{5}(2))/6)),:),1),2)),'filled');
colormap(gca,cmap);
caxis([-.5 .5])
grid on
subplot(4,2,4),...
scatter(squeeze(nanmean(nanmean(Zmoold(round(stidx{6}(1)/6):(round((stidx{6}(1)+lenidx{6}(1))/6)),round(stidx{6}(2)/6):(round((stidx{6}(2)+lenidx{6}(2))/6)),:),1),2)),...
    squeeze(nanmean(nanmean(Pmoold(round(stidx{6}(1)/6):(round((stidx{6}(1)+lenidx{6}(1))/6)),round(stidx{6}(2)/6):(round((stidx{6}(2)+lenidx{6}(2))/6)),:),1),2)),60,...
    squeeze(nanmean(nanmean(deltaOUT(round(stidx{6}(1)/6):(round((stidx{6}(1)+lenidx{6}(1))/6)),round(stidx{6}(2)/6):(round((stidx{6}(2)+lenidx{6}(2))/6)),:),1),2)),'filled');
grid on
colormap(gca,cmap);
caxis([-.5 .5])

%deltaOUT=(deltaMLD-deltaustar);
deltaOUT=deltaMLD;
subplot(4,2,5),...
scatter(squeeze(nanmean(nanmean(Zmoold(round(stidx{5}(1)/6):(round((stidx{5}(1)+lenidx{5}(1))/6)),round(stidx{5}(2)/6):(round((stidx{5}(2)+lenidx{5}(2))/6)),:),1),2)),...
    squeeze(nanmean(nanmean(Pmoold(round(stidx{5}(1)/6):(round((stidx{5}(1)+lenidx{5}(1))/6)),round(stidx{5}(2)/6):(round((stidx{5}(2)+lenidx{5}(2))/6)),:),1),2)),60,...
    squeeze(nanmean(nanmean(deltaOUT(round(stidx{5}(1)/6):(round((stidx{5}(1)+lenidx{5}(1))/6)),round(stidx{5}(2)/6):(round((stidx{5}(2)+lenidx{5}(2))/6)),:),1),2)),'filled');
colormap(gca,cmap);
caxis([-.5 .5])
grid on
subplot(4,2,6),...
scatter(squeeze(nanmean(nanmean(Zmoold(round(stidx{6}(1)/6):(round((stidx{6}(1)+lenidx{6}(1))/6)),round(stidx{6}(2)/6):(round((stidx{6}(2)+lenidx{6}(2))/6)),:),1),2)),...
    squeeze(nanmean(nanmean(Pmoold(round(stidx{6}(1)/6):(round((stidx{6}(1)+lenidx{6}(1))/6)),round(stidx{6}(2)/6):(round((stidx{6}(2)+lenidx{6}(2))/6)),:),1),2)),60,...
    squeeze(nanmean(nanmean(deltaOUT(round(stidx{6}(1)/6):(round((stidx{6}(1)+lenidx{6}(1))/6)),round(stidx{6}(2)/6):(round((stidx{6}(2)+lenidx{6}(2))/6)),:),1),2)),'filled');
grid on
colormap(gca,cmap);
caxis([-.5 .5])

deltaOUT=deltaustar;
subplot(4,2,7),...
scatter(squeeze(nanmean(nanmean(Zmoold(round(stidx{5}(1)/6):(round((stidx{5}(1)+lenidx{5}(1))/6)),round(stidx{5}(2)/6):(round((stidx{5}(2)+lenidx{5}(2))/6)),:),1),2)),...
    squeeze(nanmean(nanmean(Pmoold(round(stidx{5}(1)/6):(round((stidx{5}(1)+lenidx{5}(1))/6)),round(stidx{5}(2)/6):(round((stidx{5}(2)+lenidx{5}(2))/6)),:),1),2)),60,...
    squeeze(nanmean(nanmean(deltaOUT(round(stidx{5}(1)/6):(round((stidx{5}(1)+lenidx{5}(1))/6)),round(stidx{5}(2)/6):(round((stidx{5}(2)+lenidx{5}(2))/6)),:),1),2)),'filled');
colormap(gca,cmap);
caxis([-.5 .5])
grid on
subplot(4,2,8),...
scatter(squeeze(nanmean(nanmean(Zmoold(round(stidx{6}(1)/6):(round((stidx{6}(1)+lenidx{6}(1))/6)),round(stidx{6}(2)/6):(round((stidx{6}(2)+lenidx{6}(2))/6)),:),1),2)),...
    squeeze(nanmean(nanmean(Pmoold(round(stidx{6}(1)/6):(round((stidx{6}(1)+lenidx{6}(1))/6)),round(stidx{6}(2)/6):(round((stidx{6}(2)+lenidx{6}(2))/6)),:),1),2)),60,...
    squeeze(nanmean(nanmean(deltaOUT(round(stidx{6}(1)/6):(round((stidx{6}(1)+lenidx{6}(1))/6)),round(stidx{6}(2)/6):(round((stidx{6}(2)+lenidx{6}(2))/6)),:),1),2)),'filled');
grid on
colormap(gca,cmap);
caxis([-.5 .5])
deltaOUT=deltaMLD-deltaustar;


figure; 
OUT2=squeeze(nanmean(nanmean(deltaMLD(round(stidx{6}(1)/6):(round((stidx{6}(1)+lenidx{6}(1))/6)),round(stidx{6}(2)/6):(round((stidx{6}(2)+lenidx{6}(2))/6)),:),1),2));
OUT1=squeeze(nanmean(nanmean(deltaustar(round(stidx{6}(1)/6):(round((stidx{6}(1)+lenidx{6}(1))/6)),round(stidx{6}(2)/6):(round((stidx{6}(2)+lenidx{6}(2))/6)),:),1),2));
windowSize = 3; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;

subplot(2,2,3),...
    plot(0:13,filtfilt(b,a,[OUT1(end) OUT1' OUT1(1)]),'r-',...
    0:13,filtfilt(b,a,[OUT2(end)  OUT2'  OUT2(1)]),'b-','linewidth',2);
hold on
    %1:12,squeeze(nanmean(nanmean(deltaOUT(round(stidx{6}(1)/6):(round((stidx{6}(1)+lenidx{6}(1))/6)),round(stidx{6}(2)/6):(round((stidx{6}(2)+lenidx{6}(2))/6)),:),1),2))
xlim([1 12])
ylim([0 1])
grid on
OUT1=squeeze(nanmean(nanmean(deltaustar(round(stidx{5}(1)/6):(round((stidx{5}(1)+lenidx{5}(1))/6)),round(stidx{5}(2)/6):(round((stidx{5}(2)+lenidx{5}(2))/6)),:),1),2));
OUT2=squeeze(nanmean(nanmean(deltaMLD(round(stidx{5}(1)/6):(round((stidx{5}(1)+lenidx{5}(1))/6)),round(stidx{5}(2)/6):(round((stidx{5}(2)+lenidx{5}(2))/6)),:),1),2));
subplot(2,2,1),...
    plot(0:13,filtfilt(b,a,[OUT1(end) OUT1' OUT1(1)]),'r-',...
    0:13,filtfilt(b,a,[OUT2(end) OUT2' OUT2(1)]),'b-','linewidth',2);
hold on;
    %1:12,squeeze(nanmean(nanmean(deltaOUT(round(stidx{5}(1)/6):(round((stidx{5}(1)+lenidx{5}(1))/6)),round(stidx{5}(2)/6):(round((stidx{5}(2)+lenidx{5}(2))/6)),:),1),2)))
xlim([1 12])
ylim([0 1])
grid on
%tlon(round(stidx{6}(1)/6):(round((stidx{6}(1)+lenidx{6}(1))/6)),round(stidx{6}(2)/6):(round((stidx{6}(2)+lenidx{6}(2))/6)))
%tlat(round(stidx{6}(1)/6):(round((stidx{6}(1)+lenidx{6}(1))/6)),round(stidx{6}(2)/6):(round((stidx{6}(2)+lenidx{6}(2))/6)))

%%% load ROMS MLD
tidx1=(1781-59):(2145+60)
windowSize = 90; 
b1 = (1/windowSize)*ones(1,windowSize);
a = 1;
addpath ~/ML_diagnostics/matlabfunctions/seawater/
CASENM1={'FILE5','','FILE6'}
for itplt=[1,3]
   CASENM=CASENM1{itplt}
clear temp salt pden b dbmax idxdbmax dbdz mask z_g hsbl MXLROMS
%CASENM='FILE5'

temp=squeeze(nanmean(nanmean(ncread(strcat('/glade/work/dwhitt/roms.1dmodel/OUT/ocean_avg_roms.1dmodel_',CASENM,'.nc'),'temp'),1),2));
salt=squeeze(nanmean(nanmean(ncread(strcat('/glade/work/dwhitt/roms.1dmodel/OUT/ocean_avg_roms.1dmodel_',CASENM,'.nc'),'salt'),1),2));
salt=salt(:,tidx1);
temp=temp(:,tidx1);
addpath /glade/work/dwhitt/roms.1dmodel/INI
s_r=ncread(strcat('/glade/work/dwhitt/roms.1dmodel/OUT/ocean_avg_roms.1dmodel_',CASENM,'.nc'),'s_rho');
hc=ncread(strcat('/glade/work/dwhitt/roms.1dmodel/OUT/ocean_avg_roms.1dmodel_',CASENM,'.nc'),'hc');
h=nanmean(nanmean(ncread(strcat('/glade/work/dwhitt/roms.1dmodel/OUT/ocean_avg_roms.1dmodel_',CASENM,'.nc'),'h'),1),2);
hsbl=squeeze(nanmean(nanmean(ncread(strcat('/glade/work/dwhitt/roms.1dmodel/OUT/ocean_avg_roms.1dmodel_',CASENM,'.nc'),'Hsbl'),1),2));
theta_s=nanmean(nanmean(ncread(strcat('/glade/work/dwhitt/roms.1dmodel/OUT/ocean_avg_roms.1dmodel_',CASENM,'.nc'),'theta_s'),1),2);
theta_b=nanmean(nanmean(ncread(strcat('/glade/work/dwhitt/roms.1dmodel/OUT/ocean_avg_roms.1dmodel_',CASENM,'.nc'),'theta_b'),1),2);
z_r= zROMS1D(s_r,hc,h,theta_s,theta_b,0);
pden = sw_dens(salt,temp,zeros(size(temp)));
hsbl=squeeze(nanmean(nanmean(ncread(strcat('/glade/work/dwhitt/roms.1dmodel/OUT/ocean_avg_roms.1dmodel_',CASENM,'.nc'),'Hsbl'),1),2));
b=-9.81.*pden./1027;
[dbmax,idxdbmax]=max(-(b(1:end-5,:)-repmat((b(end,:)),[length(z_r)-5 1]))./abs(repmat(z_r(1:end-5),[1 size(b,2)])),[],1);
dbdz=(b(3:end,:)-b(1:end-2,:))./(repmat(z_r(3:end)-z_r(1:end-2),[1 size(b,2)]));
mask=dbdz<repmat(dbmax,[size(dbdz,1) 1]);
z_g=repmat(z_r(2:end-1),[1 size(b,2)]);
z_g(mask)=nan;
MXLROMSCTL=min(abs(z_g),[],1);
clear sustr svstr ustarroms
sustr=squeeze(nanmean(nanmean(ncread(strcat('/glade/work/dwhitt/roms.1dmodel/OUT/ocean_avg_roms.1dmodel_',CASENM,'.nc'),'sustr'),1),2));
svstr=squeeze(nanmean(nanmean(ncread(strcat('/glade/work/dwhitt/roms.1dmodel/OUT/ocean_avg_roms.1dmodel_',CASENM,'.nc'),'svstr'),1),2));
ustarroms=double(sqrt((sqrt(svstr.^2 + sustr.^2))./1025));
ustarromsCTL=ustarroms(tidx1);

clear temp salt pden b dbmax idxdbmax dbdz mask z_g hsbl MXLROMS
FILTNM='060'
temp=squeeze(nanmean(nanmean(ncread(strcat('/glade/work/dwhitt/roms.1dmodel/OUT/ocean_avg_roms.1dmodel_',CASENM,'_filt',FILTNM,'.nc'),'temp'),1),2));
salt=squeeze(nanmean(nanmean(ncread(strcat('/glade/work/dwhitt/roms.1dmodel/OUT/ocean_avg_roms.1dmodel_',CASENM,'_filt',FILTNM,'.nc'),'salt'),1),2));
salt=salt(:,tidx1);
temp=temp(:,tidx1);
pden = sw_dens(salt,temp,zeros(size(temp)));
hsbl=squeeze(nanmean(nanmean(ncread(strcat('/glade/work/dwhitt/roms.1dmodel/OUT/ocean_avg_roms.1dmodel_',CASENM,'_filt',FILTNM,'.nc'),'Hsbl'),1),2));
b=-9.81.*pden./1027;
[dbmax,idxdbmax]=max(-(b(1:end-5,:)-repmat((b(end,:)),[length(z_r)-5 1]))./abs(repmat(z_r(1:end-5),[1 size(b,2)])),[],1);
dbdz=(b(3:end,:)-b(1:end-2,:))./(repmat(z_r(3:end)-z_r(1:end-2),[1 size(b,2)]));
%dbdz(end-5:end,:)=0;
mask=dbdz<repmat(dbmax,[size(dbdz,1) 1]);
z_g=repmat(z_r(2:end-1),[1 size(b,2)]);
z_g(mask)=nan;
MXLROMS=min(abs(z_g),[],1);
clear sustr svstr ustarroms
clear temp salt pden b dbmax idxdbmax dbdz mask z_g hsbl
FILTNM='060'
sustr=squeeze(nanmean(nanmean(ncread(strcat('/glade/work/dwhitt/roms.1dmodel/OUT/ocean_avg_roms.1dmodel_',CASENM,'_filt',FILTNM,'.nc'),'sustr'),1),2));
svstr=squeeze(nanmean(nanmean(ncread(strcat('/glade/work/dwhitt/roms.1dmodel/OUT/ocean_avg_roms.1dmodel_',CASENM,'_filt',FILTNM,'.nc'),'svstr'),1),2));
ustarroms=double(sqrt(sqrt(svstr.^2 + sustr.^2)./1026));
ustarroms=ustarroms(tidx1);

subplot(2,2,itplt),...
plot(15./30.15+(-59:(365+59))./(30+5./30),(filtfilt(b1,a,MXLROMSCTL)-filtfilt(b1,a,MXLROMS))./filtfilt(b1,a,MXLROMSCTL),'b--','linewidth',1);
%xlim([15 350]./(30+5./30))
xlim([1 12])
subplot(2,2,itplt),...
plot(15./30.15+(-59:(365+59))./(30+5./30),(filtfilt(b1,a,ustarromsCTL)-filtfilt(b1,a,ustarroms))./filtfilt(b1,a,ustarromsCTL),'r--','linewidth',1);
%xlim([15 350]./(30+5./30))
xlim([1 12])
%wts=[31 28 31 30 31 30 31 31 30 31 30 31]';
%set(gca,'XTick',[0.5:11.5]);
end


CASENM1={'','FILE5','','FILE6'}
for itplt=[2,4]
   CASENM=CASENM1{itplt}
filtname={'120','060','030','015','010','005','002','001'}
for it=1:8
FILTNM=filtname{it}
clear temp salt pden b dbmax idxdbmax dbdz mask z_g hsbl MXLROMS ustarroms
%CASENM='FILE5'
temp=squeeze(nanmean(nanmean(ncread(strcat('/glade/work/dwhitt/roms.1dmodel/OUT/ocean_avg_roms.1dmodel_',CASENM,'_filt',FILTNM,'.nc'),'temp'),1),2));
salt=squeeze(nanmean(nanmean(ncread(strcat('/glade/work/dwhitt/roms.1dmodel/OUT/ocean_avg_roms.1dmodel_',CASENM,'_filt',FILTNM,'.nc'),'salt'),1),2));
salt=salt(:,tidx1);
temp=temp(:,tidx1);
pden = sw_dens(salt,temp,zeros(size(temp)));
hsbl=squeeze(nanmean(nanmean(ncread(strcat('/glade/work/dwhitt/roms.1dmodel/OUT/ocean_avg_roms.1dmodel_',CASENM,'_filt',FILTNM,'.nc'),'Hsbl'),1),2));
b=-9.81.*pden./1027;
[dbmax,idxdbmax]=max(-(b(1:end-5,:)-repmat((b(end,:)),[length(z_r)-5 1]))./abs(repmat(z_r(1:end-5),[1 size(b,2)])),[],1);
dbdz=(b(3:end,:)-b(1:end-2,:))./(repmat(z_r(3:end)-z_r(1:end-2),[1 size(b,2)]));
%dbdz(end-5:end,:)=0;
mask=dbdz<repmat(dbmax,[size(dbdz,1) 1]);
z_g=repmat(z_r(2:end-1),[1 size(b,2)]);
z_g(mask)=nan;
MXLROMS=min(abs(z_g),[],1);
clear sustr svstr ustarroms
clear temp salt pden b dbmax idxdbmax dbdz mask z_g hsbl
sustr=squeeze(nanmean(nanmean(ncread(strcat('/glade/work/dwhitt/roms.1dmodel/OUT/ocean_avg_roms.1dmodel_',CASENM,'_filt',FILTNM,'.nc'),'sustr'),1),2));
svstr=squeeze(nanmean(nanmean(ncread(strcat('/glade/work/dwhitt/roms.1dmodel/OUT/ocean_avg_roms.1dmodel_',CASENM,'_filt',FILTNM,'.nc'),'svstr'),1),2));
ustarroms=double(sqrt(sqrt(svstr.^2 + sustr.^2)./1026));
ustarroms=ustarroms(tidx1);
cmap=colormap(cmocean('thermal',21));
subplot(2,2,itplt),...
plot(15./30.15+(-59:(365+59))./(30+5./30),((filtfilt(b1,a,MXLROMSCTL)-filtfilt(b1,a,MXLROMS))./filtfilt(b1,a,MXLROMSCTL)),...
'color',cmap(2*it,:),'linewidth',1);
hold on
plot(15./30.15+(-59:(365+59))./(30+5./30),((filtfilt(b1,a,ustarromsCTL')-filtfilt(b1,a,ustarroms'))./filtfilt(b1,a,ustarromsCTL')),'--','color',cmap(2*it,:),'linewidth',1);

hold on
%xlim([15 350]./(30+5./30))
xlim([1 12])
%wts=[31 28 31 30 31 30 31 31 30 31 30 31]';
%set(gca,'XTick',[0.5:11.5]);
end
end

pause




X=Zmo;
y=real(log10(sqrt(N2barXMXLp2monmean016)./repmat(abs(coriolis),[1 1 12])));
clear Z

figure;
Z = [X(~isnan(X+y)) y(~isnan(X+y))];
hist3(Z,'Ctrs',{(-30.5:.5:30.5) (1:.05:4)},'CdataMode','auto')
shading flat
%colormap(
xlabel('Z')
ylabel('P')
colorbar
view(2)
pause(1)
grid on
set(gca,'layer','top')
caxis([1 5E4])
cmapCMR=CMRmap(128);
colormap(flipud(cmapCMR))
set(gca,'colorscale','log')
    set(gca,'layer','top')
    set(gcf,'color','w')
    hold on
xlim([-30 30])
ylim([1 4])
title(['P-Z diagram from monthly climatologies in CTL'],'fontweight','normal','fontsize',13)
set(gca,'Fontsize',13)
export_fig('fig_PZmo_diagram.png','-r300')
clear X y
pause;

clear y3
diffmostdyrs=squeeze(nanstd(diffmo,0,4));
diffmotemp=diffmo;
diffmotemp(diffmotemp==0)=nan;
ZtestMLDdiff=(diffmotemp-repmat(diffmomeanyrs,[1 1 1 5]))./(repmat(diffmostdyrs,[1 1 1 5]));
y3=diffmomeanyrs./(abs(varoutoldclmo)+epsilon1);

%load('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/ustar_end_workspace_nobigvar.mat','diffmo','varoutoldclmo','epsilon1','tlon','tlat');
load('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholsonCarranza_public/fig_1_fricvel_fig_9_load.mat','diffmo','varoutoldclmo','epsilon1','tlon','tlat');
clear diffmomeanyrs diffmostdyrs X3
diffmomeanyrs=squeeze(nanmean(diffmo,4));
diffmostdyrs=squeeze(nanstd(diffmo,0,4));
diffmotemp=diffmo;
diffmotemp(diffmotemp==0)=nan;
Ztestfdiff=(diffmotemp-repmat(diffmomeanyrs,[1 1 1 5]))./(repmat(diffmostdyrs,[1 1 1 5]));

clear tempvar tempmeanvar
tempvar=diffmo;
tempvar=reshape(tempvar,[size(tempvar,1)*size(tempvar,2)*size(tempvar,3) size(tempvar,4)]);
tempmeanvar=remove_insignificant_points(tempvar,0);
diffmomeanyrs=reshape(tempmeanvar,[size(diffmo,1) size(diffmo,2) size(diffmo,3)]);
X3=diffmomeanyrs./(abs(varoutoldclmo)+epsilon1);


y4=(y3-X3);
Zmo(abs(Zmo)>10)=nan;
Pmo(abs(Pmo)>log10(316))=nan;
Zmo(abs(Pmo)>log10(316))=nan;
Pmo(abs(Zmo)>10)=nan;
figure('Position',[100 100 1000 400]);
subplot(1,2,1),...
clear Z
Z = [Pmo(~isnan(Pmo+y4)) y4(~isnan(Pmo+y4))];
%hist3(Z,'Ctrs',{(-50:.2:50) (-98.75:2.5:98.75)./100},'CdataMode','auto')
hist3(Z,'Ctrs',{(0:.02:4) (-98.75:2.5:98.75)./100},'CdataMode','auto')

shading flat
xlabel('Z')
ylabel('\Delta MLD')
colorbar
view(2)
pause(1)
grid on
set(gca,'layer','top')
caxis([1 0.5E4])
cmapCMR=CMRmap(128);
colormap(flipud(cmapCMR))
set(gca,'colorscale','log')
set(gca,'colorscale','linear')
    set(gca,'layer','top')
    set(gcf,'color','w')
    hold on
    %plot3(linspace(1,3.5,10),linspace(-.5,.5,10),5e4.*ones(10,1),'k-')
xlim([1 log10(316)])
ylabel('\Delta MLD - \Delta u_*')
xlabel('P')
ylim([-.5 .5])
[r,p]=corrcoef(Pmo(~isnan(Pmo+y4)),y4(~isnan(Pmo+y4)))
title(['(A) Relationship , r=',num2str(r(1,2),3)],'fontweight','normal','fontsize',10)
axis square
subplot(1,2,2),...
clear Z
Z = [Zmo(~isnan(Zmo+y4)) y4(~isnan(Zmo+y4))];
hist3(Z,'Ctrs',{(-50:.2:50) (-98.75:2.5:98.75)./100},'CdataMode','auto')
%hist3(Z,'Ctrs',{(0:.02:4) (-98.75:2.5:98.75)./100},'CdataMode','auto')

shading flat
xlabel('Z')
colorbar
view(2)
pause(1)
grid on
set(gca,'layer','top')
caxis([1 0.5E4])
cmapCMR=CMRmap(128);
colormap(flipud(cmapCMR))
set(gca,'colorscale','log')
set(gca,'colorscale','linear')
    set(gca,'layer','top')
    set(gcf,'color','w')
    hold on
  %  plot3(linspace(-5,5,10),linspace(-.5,.5,10),5e4.*ones(10,1),'k-')
xlim([-10 10])
ylabel('\Delta MLD - \Delta u_*')
ylim([-.5 .5])
[r,p]=corrcoef(Zmo(~isnan(Zmo+y4)),y4(~isnan(Zmo+y4)))
title(['(B) Relationship , r=',num2str(r(1,2),3)],'fontweight','normal','fontsize',10)
axis square
export_fig('fig_PmoandZmovsDeltaMLD.png','-r300')
pause


clear Z
figure;
subplot(2,1,1),...
Z = [X3(~isnan(X3+y3)) y3(~isnan(X3+y3))];
hist3(Z,'Ctrs',{(-98.75:2.5:98.75)./100 (-98.75:2.5:98.75)./100},'CdataMode','auto')
shading flat
%colormap(
xlabel('\Delta u_*')
ylabel('\Delta MLD')
colorbar
view(2)
pause(1)
grid on
set(gca,'layer','top')
caxis([1 3E4])
cmapCMR=CMRmap(128);
colormap(flipud(cmapCMR))
set(gca,'colorscale','log')
    set(gca,'layer','top')
    set(gcf,'color','w')
    hold on
    plot3(linspace(-.9875,.9875,10),linspace(-.9875,.9875,10),5e4.*ones(10,1),'k-')
xlim([-.9875 .9875])
ylim([-.9875 .9875])
[r,p]=corrcoef(X3(~isnan(X3+y3)),y3(~isnan(X3+y3)))
title(['(A) Normalized monthly mean sensitivity of MLD to friction velocity, r=',num2str(r(1,2),3)],'fontweight','normal','fontsize',10)



%load('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/MLD_data.mat','diffmomean','varoutoldclmomean','epsilon1');
load('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholsonCarranza_public/fig_7_XMXL_fig_9_load.mat','diffmomean','varoutoldannmean','epsilon1');
varoutoldclmomean=varoutoldannmean;
clear y3
diffmomean(diffmomean==0)=nan;
y3=diffmomean./(abs(varoutoldclmomean)+epsilon1);
%load('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/ustar_end_workspace_nobigvar.mat','diffmomean','varoutoldclmomean','epsilon1');
load('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholsonCarranza_public/fig_1_fricvel_fig_9_load.mat','diffmomean','varoutoldannmean','epsilon1');
varoutoldclmomean=varoutoldannmean;
clear X3
diffmomean(diffmomean==0)=nan;
X3=diffmomean./(abs(varoutoldclmomean)+epsilon1);
axis square

clear Z
subplot(2,1,2),...
Z = [X3(~isnan(X3+y3)) y3(~isnan(X3+y3))];
hist3(Z,'Ctrs',{(-98.75:2.5:98.75)./100 (-98.75:2.5:98.75)./100},'CdataMode','auto')
shading flat
%colormap(
xlabel('\Delta u_*')
ylabel('\Delta MLD')
colorbar
view(2)
pause(1)
grid on
set(gca,'layer','top')
caxis([1 3E3])
cmapCMR=CMRmap(128);
colormap(flipud(cmapCMR))
set(gca,'colorscale','log')
    set(gca,'layer','top')
    set(gcf,'color','w')
    hold on
    plot3(linspace(-.9875,.9875,10),linspace(-.9875,.9875,10),5e4.*ones(10,1),'k-')
xlim([-.9875 .9875])
ylim([-.9875 .9875])
clear r p
[r,p]=corrcoef(X3(~isnan(X3+y3)),y3(~isnan(X3+y3)))
title(['(B) Normalized annual mean sensitivity of MLD to friction velocity, r=',num2str(r(1,2),3)],'fontweight','normal','fontsize',10)
set(gcf,'color','w')
axis square

%%
figure;
load coast
axm=axesm('MapProjection','robinson');
surfm(tlat,tlon,(y3)); % this makes a color-filled plot
patchesm(lat,long,0*[1 1 1]);    % this makes grey shading over land
framem;
gridm('on')
gridm(':')
gridm('GLineWidth',1)
axis off
h = colorbar('h','location','eastoutside')
title(strcat(''));
caxis([-.5 .5])
addpath ~/cmocean_v1
cmap=cmocean('balance');
% plot lat/lon of 1d models
latlon1d1(1,:)=[330.04999999999995 -31.020715042064317 5.0 3.0];
latlon1d1(2,:)=[330.04999999999995 -47.87347880827708 5.0 3.0];
hold on
plotm(latlon1d1(:,2)+latlon1d1(:,4)./2,latlon1d1(:,1)+latlon1d1(:,3)./2,'magentao');
colormap(gca,cmap)

figure;
load coast
axm=axesm('MapProjection','robinson');
surfm(tlat,tlon,(X3)); % this makes a color-filled plot
patchesm(lat,long,0*[1 1 1]);    % this makes grey shading over land
framem;
gridm('on')
gridm(':')
gridm('GLineWidth',1)
axis off
h = colorbar('h','location','eastoutside')
title(strcat(''));
caxis([-.5 .5])
addpath ~/cmocean_v1
cmap=cmocean('balance');
hold on
plotm(latlon1d1(:,2)+latlon1d1(:,4)./2,latlon1d1(:,1)+latlon1d1(:,3)./2,'magentao');
colormap(gca,cmap)

figure;
load coast
axm=axesm('MapProjection','robinson');
%surfm(tlat,tlon,(X3-y3)./(abs(y3)+0)); % this makes a color-filled plot
surfm(tlat,tlon,(y3-X3)); % this makes a color-filled plot

patchesm(lat,long,0*[1 1 1]);    % this makes grey shading over land
framem;
gridm('on')
gridm(':')
gridm('GLineWidth',1)
axis off
h = colorbar('h','location','eastoutside')
title(strcat(''));
caxis([-0.5 0.5])
addpath ~/cmocean_v1
cmap=cmocean('balance');
colormap(gca,cmap)
%title('Normalized (\Delta u_* - \Delta MLD)/(\Delta MLD)')
title('\Delta MLD-\Delta u_*')
hold on
plotm(latlon1d1(:,2)+latlon1d1(:,4)./2,latlon1d1(:,1)+latlon1d1(:,3)./2,'magentao');

set(gcf,'color','w')

load('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholsonCarranza_public/fig_2_fricvel_fig_9_load.mat','fractotallmean');
load('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholsonCarranza_public/fig_1_stress_fig_9_load.mat','varoutoldannmean');
CVstress=fractotallmean./varoutoldclmomean;
%[r,p]=corrcoef(CVstress(:),(y3(:)-X3(:)));
figure;
load coast
axm=axesm('MapProjection','robinson');
%surfm(tlat,tlon,(X3-y3)./(abs(y3)+0)); % this makes a color-filled plot
surfm(tlat,tlon,CVstress); % this makes a color-filled plot

patchesm(lat,long,0*[1 1 1]);    % this makes grey shading over land
framem;
gridm('on')
gridm(':')
gridm('GLineWidth',1)
axis off
h = colorbar('h','location','eastoutside')
title(strcat(''));
%caxis([-0.5 0.5])
addpath ~/cmocean_v1
cmap=cmocean('balance');
colormap(gca,cmap)