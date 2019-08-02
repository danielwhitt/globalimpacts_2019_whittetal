%% plotting
if epsilon1 ~= 0
    disp('epsilon WARNING')
    pause
end
figure('Position',[1 1 800 800]);
pause(.5)
subplot(3,2,1),...
    load coast
axm=axesm('MapProjection','robinson');
surfm(tlat,tlon,log10(fractotallmean)); % this makes a color-filled plot
patchesm(lat,long,0*[1 1 1]);    % this makes grey shading over land
framem;
gridm('on')
gridm(':')
gridm('GLineWidth',1)
axis off
h = colorbar('h','location','eastoutside')
ylabel(h,cblabel1)
title(['(A) Subannual std., ',varoutname]);
colormap(gca,jet);
set(gcf,'color','w')
caxis(caxis1)
colormap(gca,jet)

subplot(3,2,2),...
    load coast
axm=axesm('MapProjection','robinson');
surfm(tlat,tlon,difftotvarmomean./(abs(fractotallmean)+epsilon1)); % this makes a color-filled plot
patchesm(lat,long,0*[1 1 1]);    % this makes grey shading over land
framem;
gridm('on')
gridm(':')
gridm('GLineWidth',1)
axis off
h = colorbar('h','location','eastoutside')
title(strcat('(B) CTL-LP/(|CTL|)'));
caxis([-1 1])
addpath ~/cmocean_v1
cmap=cmocean('balance');
colormap(gca,cmap)

subplot(3,2,3),...
    load coast
axm=axesm('MapProjection','robinson');
surfm(tlat,tlon,log10(fracintravarmean)); % this makes a color-filled plot
patchesm(lat,long,0*[1 1 1]);    % this makes grey shading over land
framem;
gridm('on')
gridm(':')
gridm('GLineWidth',1)
axis off
h = colorbar('h','location','eastoutside')
ylabel(h,cblabel1)
title(['(C) Subseasonal std., ',varoutname]);
colormap(gca,jet);
set(gcf,'color','w')
caxis(caxis1)
colormap(gca,jet)

subplot(3,2,4),...
    load coast
axm=axesm('MapProjection','robinson');
surfm(tlat,tlon,diffintravarmomean./(abs(fracintravarmean)+epsilon1)); % this makes a color-filled plot
patchesm(lat,long,0*[1 1 1]);    % this makes grey shading over land
framem;
gridm('on')
gridm(':')
gridm('GLineWidth',1)
axis off
h = colorbar('h','location','eastoutside')
title(strcat('(D) CTL-LP/(|CTL|)'));
caxis([-1 1])
addpath ~/cmocean_v1
cmap=cmocean('balance');
colormap(gca,cmap)

subplot(3,2,5),...    
loglog(abs(freq./(2.*pi)),areaweightedoldPSmean,'b-',abs(freq./(2.*pi)),areaweightedPSmean,'r-','linewidth',2)
hold on
legend('CTL','LP')
h=area(abs(freq./(2.*pi)),[(varoutoldPSiqr(1,:))' (varoutoldPSiqr(2,:)-varoutoldPSiqr(1,:))'],'linestyle','none','handlevisibility','off')
h(1).FaceColor='w';
h(2).FaceColor='b';
h(1).FaceAlpha=0;
h(2).FaceAlpha=0.2;
set(gca,'Xscale','log')
set(gca,'Yscale','log')
hold on
h=area(abs(freq./(2.*pi)),[(varoutPSiqr(1,:))' (varoutPSiqr(2,:)-varoutPSiqr(1,:))'],'linestyle','none','handlevisibility','off')
h(1).FaceColor='w';
h(2).FaceColor='r';
h(1).FaceAlpha=0;
h(2).FaceAlpha=0.2;
set(gca,'Xscale','log')
set(gca,'Yscale','log')
xlim([1./365 1])
hold on
hold on
grid on
xlabel('cyc/day')
ylabel(ylabel2)
ylim([1e-5.*((10.^caxis1(2)).^2) ((10.^caxis1(2)).^2)])
title('(E) Global area average power spectrum')

subplot(3,2,6),...
tareagnow=tareag;
tareagnow(isnan(squeeze(fracvarmean(:,:,1))))=nan;
   loglog(abs(freq((nfreq/2+2):end)./(2.*pi)),squeeze(nansum(reshape(fracvarmean.*tareagnow(:,:,(nfreq/2+2):end),[size(varoutPS,1)*size(varoutPS,2) (nfreq/2-1)]),1)./...
   nansum(reshape(tareagnow(:,:,(nfreq/2+2):end),[size(varoutPS,1)*size(varoutPS,2) (nfreq/2-1)]),1)),'linewidth',2,'color','r')
hold on
tareagnow=tareag;
tareagnow(isnan(fracvaroldmean))=nan;
   loglog(abs(freq((nfreq/2+2):end)./(2.*pi)),squeeze(nansum(reshape(fracvaroldmean.*tareagnow(:,:,(nfreq/2+2):end),[size(varoutPS,1)*size(varoutPS,2) (nfreq/2-1)]),1)./...
   nansum(reshape(tareagnow(:,:,(nfreq/2+2):end),[size(varoutPS,1)*size(varoutPS,2) (nfreq/2-1)]),1)),'linewidth',2,'color','b')

grid on
xlabel('cyc/day')
%ylabel('[N/(m^2)]^2')
xlim([1./365 1])
ylim([1e-3 1])
set(gca,'Yscale','linear');
ylim([0 1])
title('(F) Cumulative sum, fraction of variance')


export_fig(exportname,'-r300')
%%
figure('Position',[1 1 640 480]);
   load coast
axm=axesm('MapProjection','robinson');
surfm(tlat,tlon,(fracintravarmean./fractotallmean).^2); % this makes a color-filled plot
patchesm(lat,long,0*[1 1 1]);    % this makes grey shading over land
framem;
gridm('on')
gridm(':')
gridm('GLineWidth',1)
axis off
h = colorbar('h','location','eastoutside')
%ylabel(h,cblabel1)
title(['Fraction of variance that is subseasonal, ',varoutname]);
colormap(gca,jet);
set(gcf,'color','w')
caxis([0 1])
colormap(gca,flipud(CMRmap(64)))
%colormap(gca,jet)
export_fig(exportname2,'-r300')
%%%
anglediff=squeeze(angdiff(varoutoldPhS(:,:,364/2+2,:),varoutPhS(:,:,364/2+2,:)));
anglediffmean=nanmean(anglediff,3);
anglediffstd=nanstd(anglediff,0,3);
clear tout pvalue
tout=sqrt(5)*abs(anglediffmean)./anglediffstd;
pvalue=1-tcdf(tout,5-1);
% set pvalue threshold for insignificance
mask=pvalue>.05;
anglediffmean(mask)=nan;

figure('Position',[1 1 640 480]);
   load coast
axm=axesm('MapProjection','robinson');
surfm(tlat,tlon,-anglediffmean); % this makes a color-filled plot
patchesm(lat,long,0*[1 1 1]);    % this makes grey shading over land
framem;
gridm('on')
gridm(':')
gridm('GLineWidth',1)
axis off
h = colorbar('h','location','eastoutside','Ticks',[-1.047 -.5236 0 .5236 1.047],...
    'TickLabels',{'-2','-1','0','1','2'})
%ylabel(h,cblabel1)
title(['Change in phase of seasonal cycle, CTL - LP, ',varoutname]);
colormap(gca,jet);
set(gcf,'color','w')
%caxis([0 1])
caxis([-pi/2.4 pi/2.4])
colormap(gca,cmocean('balance',5))
export_fig(exportname3,'-r300');
