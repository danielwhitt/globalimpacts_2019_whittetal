%% plotting
figure('Position',[10 10 1000 550]);
pause(.5)
maskhotspots(isnan(tarea))=nan;
maskhotspots(maskhotspots<0)=nan;
maskhotspots(mod(maskhotspots,1)~=0)=nan;

subplot(2,3,1),...
    load coast
axm=axesm('MapProjection','robinson');
surfm(tlat,tlon,maskhotspots-0.5); % this makes a color-filled plot
patchesm(lat,long,0*[1 1 1]);    % this makes grey shading over land
framem;
caxis([0 5])
gridm('on')
gridm(':')
gridm('GLineWidth',1)
axis off
%h = colorbar('h','location','eastoutside')
%ylabel(h,cblabel1)
title(['(A) Regional Hotspots']);
cmap=colormap(jet(64));
cmap(1,:)=[1 1 1];
cmap(end,:)=[1 1 1];
colormap(gca,cmap);
set(gcf,'color','w')

for it = 1:5

clear mask maskg
mask=maskhotspots;
mask(mask>it+.1)=nan;
mask(mask<it-.1)=nan;
mask(~isnan(mask))=1;
maskg=repmat(mask,[1 1 size(tareag,3)]);
tareagnow=tareag;
tareagnow(isnan(varoutPSmean))=nan;

areaweightedPSmean=squeeze(nansum(reshape(2.*abs(varoutPSmean).*tareagnow.*maskg,...
    [size(varoutPSmean,1)*size(varoutPSmean,2) size(varoutPSmean,3)]),1)./...
    nansum(reshape(tareagnow.*maskg,[size(varoutPSmean,1)*size(varoutPSmean,2) size(varoutPSmean,3)]),1));
tareagnow=tareag;
tareagnow(isnan(varoutoldPSmean))=nan;
areaweightedoldPSmean=squeeze(nansum(reshape(2.*abs(varoutoldPSmean).*tareagnow.*maskg,...
    [size(varoutPSmean,1)*size(varoutPSmean,2) size(varoutPSmean,3)]),1)./...
    nansum(reshape(tareagnow.*maskg,[size(varoutPSmean,1)*size(varoutPSmean,2) size(varoutPSmean,3)]),1));

varoutoldPSiqr= quantile(squeeze(reshape(...
    (2.*abs(varoutoldPSmean.*maskg)),...
    [size(varoutPSmean,1)*size(varoutPSmean,2) size(varoutPSmean,3)])),[.05 .95],1);

varoutPSiqr= quantile(squeeze(reshape(...
    (2.*abs(varoutPSmean.*maskg)),...
    [size(varoutPSmean,1)*size(varoutPSmean,2) size(varoutPSmean,3)])),[.05 .95], 1);
   
 starttitlearr={'(B) Southern Ocean','(C) Kuroshio','(D) Labrador Sea','(E) Gulf Stream','(F) GIN Seas'};   
subplot(2,3,it+1),...    
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
ylim([(1e-3.*(10.^caxis1(2)).^2) (1e2.*(10.^caxis1(2)).^2)])
title([starttitlearr{it}])

end


export_fig('fig_MLD_hotspotspectra.png','-r300')


