% run after fig_9_ustardiffvsMLDdiff_CI.m

figure; 
odallflag =0
OUT2=squeeze(nanmean(nanmean(deltaMLD(round(stidx{6}(1)/6):(round((stidx{6}(1)+lenidx{6}(1))/6)),round(stidx{6}(2)/6):(round((stidx{6}(2)+lenidx{6}(2))/6)),:),1),2));
OUT1=squeeze(nanmean(nanmean(deltaustar(round(stidx{6}(1)/6):(round((stidx{6}(1)+lenidx{6}(1))/6)),round(stidx{6}(2)/6):(round((stidx{6}(2)+lenidx{6}(2))/6)),:),1),2));
windowSize = 3; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
if odallflag ==1
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
    plot(0:13,filtfilt(b,a,[OUT1(end) OUT1' OUT1(1)]),'r-',... % ustar
    0:13,filtfilt(b,a,[OUT2(end) OUT2' OUT2(1)]),'b-','linewidth',2); % MLD
hold on;
    %1:12,squeeze(nanmean(nanmean(deltaOUT(round(stidx{5}(1)/6):(round((stidx{5}(1)+lenidx{5}(1))/6)),round(stidx{5}(2)/6):(round((stidx{5}(2)+lenidx{5}(2))/6)),:),1),2)))
xlim([1 12])
ylim([0 1])
grid on
else
 subplot(2,1,2),...
    plot(0:13,filtfilt(b,a,[OUT1(end) OUT1' OUT1(1)]),'r-',... % ustar
    0:13,filtfilt(b,a,[OUT2(end)  OUT2'  OUT2(1)]),'b-','linewidth',2); %MLD
hold on
    %1:12,squeeze(nanmean(nanmean(deltaOUT(round(stidx{6}(1)/6):(round((stidx{6}(1)+lenidx{6}(1))/6)),round(stidx{6}(2)/6):(round((stidx{6}(2)+lenidx{6}(2))/6)),:),1),2))
xlim([1 12])
ylim([0 1])
grid on
OUT1=squeeze(nanmean(nanmean(deltaustar(round(stidx{5}(1)/6):(round((stidx{5}(1)+lenidx{5}(1))/6)),round(stidx{5}(2)/6):(round((stidx{5}(2)+lenidx{5}(2))/6)),:),1),2));
OUT2=squeeze(nanmean(nanmean(deltaMLD(round(stidx{5}(1)/6):(round((stidx{5}(1)+lenidx{5}(1))/6)),round(stidx{5}(2)/6):(round((stidx{5}(2)+lenidx{5}(2))/6)),:),1),2));
title('(B) \Delta MLD and \Delta u_*, Subpolar SATL, 3D vs 1D')
xlabel('Months')
ylabel('(CTL-LP)/CTL')

subplot(2,1,1),...
    plot(0:13,filtfilt(b,a,[OUT1(end) OUT1' OUT1(1)]),'r-',... % ustar
    0:13,filtfilt(b,a,[OUT2(end) OUT2' OUT2(1)]),'b-','linewidth',2); %MLD
hold on;
    %1:12,squeeze(nanmean(nanmean(deltaOUT(round(stidx{5}(1)/6):(round((stidx{5}(1)+lenidx{5}(1))/6)),round(stidx{5}(2)/6):(round((stidx{5}(2)+lenidx{5}(2))/6)),:),1),2)))
xlim([1 12])
ylim([0 1])
grid on 
title('(A) \Delta MLD and \Delta u_*, Subtropical SATL, 3D vs 1D')
xlabel('Months')
ylabel('(CTL-LP)/CTL')
end

%%% load ROMS MLD
tidx1=(1781-59):(2145+60)
windowSize = 90; 
b1 = (1/windowSize)*ones(1,windowSize);
a = 1;
addpath ~/ML_diagnostics/matlabfunctions/seawater/
if odallflag == 1
    idx1=[1,3]
    CASENM1={'FILE5','','FILE6'}
else
    idx1 = [1,2]
    CASENM1={'FILE5','FILE6'}

end
for itplt=idx1
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
if odallflag == 1
subplot(2,2,itplt),...
plot(15./30.15+(-59:(365+59))./(30+5./30),(filtfilt(b1,a,MXLROMSCTL)-filtfilt(b1,a,MXLROMS))./filtfilt(b1,a,MXLROMSCTL),'b--','linewidth',1);
%xlim([15 350]./(30+5./30))
xlim([1 12])
subplot(2,2,itplt),...
plot(15./30.15+(-59:(365+59))./(30+5./30),(filtfilt(b1,a,ustarromsCTL)-filtfilt(b1,a,ustarroms))./filtfilt(b1,a,ustarromsCTL),'r--','linewidth',1);
else
subplot(2,1,itplt),...
plot(15./30.15+(-59:(365+59))./(30+5./30),(filtfilt(b1,a,ustarromsCTL)-filtfilt(b1,a,ustarroms))./filtfilt(b1,a,ustarromsCTL),'r--','linewidth',1); % ustar
subplot(2,1,itplt),...
plot(15./30.15+(-59:(365+59))./(30+5./30),(filtfilt(b1,a,MXLROMSCTL)-filtfilt(b1,a,MXLROMS))./filtfilt(b1,a,MXLROMSCTL),'b--','linewidth',1); % MLD
%xlim([15 350]./(30+5./30))
xlim([1 12])
if itplt==1
legend('\Delta u_* 3D','\Delta MLD 3D','\Delta u_* 1D','\Delta MLD 1D')
end
end
%xlim([15 350]./(30+5./30))
xlim([1 12])
%wts=[31 28 31 30 31 30 31 31 30 31 30 31]';
%set(gca,'XTick',[0.5:11.5]);
end

if odallflag==1

CASENM1={'','FILE5','','FILE6'}
for itplt=(idx1+1)
   CASENM=CASENM1{itplt}
   
   clear temp salt pden b dbmax idxdbmax dbdz mask z_g hsbl MXLROMS MXLROMSCTL ustarromsCTL
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
ylim([0 1])
grid on

end
end