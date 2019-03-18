% heat map
% run fig 7 first
clear all
close all
load('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/MLD_data.mat','diffmo','varoutoldclmo','epsilon1');
clear diffmomeanyrs diffmostdyrs y3
diffmomeanyrs=squeeze(nanmean(diffmo,4));
diffmostdyrs=squeeze(nanstd(diffmo,0,4));
clear tout pvalue
tout=sqrt(5)*abs(diffmomeanyrs)./diffmostdyrs;
pvalue=1-tcdf(tout,5-1);
% set pvalue threshold for insignificance
mask=pvalue>.05;
diffmomeanyrs(mask)=nan;
y3=diffmomeanyrs./(abs(varoutoldclmo)+epsilon1);

load('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/ustar_end_workspace_nobigvar.mat','diffmo','varoutoldclmo','epsilon1','tlon','tlat');
clear diffmomeanyrs diffmostdyrs X3
diffmomeanyrs=squeeze(nanmean(diffmo,4));
diffmostdyrs=squeeze(nanstd(diffmo,0,4));
clear tout pvalue
tout=sqrt(5)*abs(diffmomeanyrs)./diffmostdyrs;
pvalue=1-tcdf(tout,5-1);
% set pvalue threshold for insignificance
mask=pvalue>.05;
diffmomeanyrs(mask)=nan;
X3=diffmomeanyrs./(abs(varoutoldclmo)+epsilon1);


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
title(['(A) Normalized monthly mean sensitivity of MLD to u_*, r^2=',num2str(r(1,2).^2,3)],'fontweight','normal','fontsize',10)



load('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/MLD_data.mat','diffmomean','varoutoldclmomean','epsilon1');
clear y3
y3=diffmomean./(abs(varoutoldclmomean)+epsilon1);
load('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/ustar_end_workspace_nobigvar.mat','diffmomean','varoutoldclmomean','epsilon1');
clear X3
X3=diffmomean./(abs(varoutoldclmomean)+epsilon1);


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
title(['(B) Normalized annual mean sensitivity of MLD to u_*, r^2=',num2str(r(1,2).^2,3)],'fontweight','normal','fontsize',10)

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
colormap(gca,cmap)

figure;
load coast
axm=axesm('MapProjection','robinson');
surfm(tlat,tlon,(X3-y3)./(abs(y3)+.05)); % this makes a color-filled plot
patchesm(lat,long,0*[1 1 1]);    % this makes grey shading over land
framem;
gridm('on')
gridm(':')
gridm('GLineWidth',1)
axis off
h = colorbar('h','location','eastoutside')
title(strcat(''));
caxis([-1 1])
addpath ~/cmocean_v1
cmap=cmocean('balance');
colormap(gca,cmap)
title('Normalized (\Delta u_* - \Delta MLD)/(\Delta MLD + \epsilon)')