
figure('Position',[10 10 1200 650]);
pause(.5)
for ik = 1:4
% plot
subplot(2,2,ik),...
load coast
axm=axesm('MapProjection','robinson');
if ik == 1
surfm(tlat,tlon,varoutoldannmean); % this makes a color-filled plot
elseif ik== 2
surfm(tlat,tlon,diffmomean./(abs(varoutoldannmean)+epsilon1)); % this makes a color-filled plot
elseif ik==3
surfm(tlat,tlon,seasoldmomean); % this makes a color-filled plot
elseif ik==4
surfm(tlat,tlon,diffseasmomean./(abs(seasoldmomean)+epsilon1)); % this makes a color-filled plot
end
patchesm(lat,long,0*[1 1 1]);    % this makes grey shading over land
framem;
gridm('on')
%gridm(':')
gridm('GLineWidth',1)
axis off
h = colorbar('h','location','eastoutside')
if ik == 1 || ik == 3
ylabel(h,cblabel1)
end
if ik == 1
title(['(A) Mean ',varnameout,', CTL'],'fontsize',11);
caxis(caxis1)
elseif ik == 2
title(strcat('(B) (CTL-LP)/(|CTL|+\epsilon)'),'fontsize',11);
%caxis([-1e-8 1e-8])
caxis([-1 1])
elseif ik == 3
title(strcat('(C) Seasonal cycle amplitude,',' CTL'),'fontsize',11);
caxis(caxis2)
elseif ik == 4
title(strcat('(D) (CTL-LP)/(|CTL|+\epsilon)'),'fontsize',11);
%caxis([-5e-8 5e-8])
caxis([-1 1])
end
if ik == 1 || ik ==3
colormap(gca,jet);
else
addpath ~/cmocean_v1
cmap=cmocean('balance');
colormap(gca,cmap)
end
set(gca,'fontsize',11);
set(gcf,'color','w')
end
pause(.1)
export_fig(exportname,'-r300')