%% bflux uncertainty
restoredefaultpath
clear all
close all;
%pause
addpath ./export_fig-master
addpath ./utility
addpath ./cmocean_v1
tic
loaddataflag =1
if loaddataflag == 1

fn0=strcat('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/',...
   'g.e20b07.2000_DATM%NYF_SLND_CICE_POP2_DROF%NYF_SGLC_SWAV.T62_t13.hybrid.016.pop.h.nday1.0012-12-01.nc');
path='/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/archive/g.e20b07.2000_DATM%NYF_SLND_CICE_POP2_DROF%NYF_SGLC_SWAV.T62_t13.hybrid.036/ocn/hist/'
fnstart='g.e20b07.2000_DATM%NYF_SLND_CICE_POP2_DROF%NYF_SGLC_SWAV.T62_t13.hybrid.036.pop.h.00'
fn{1} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/LP/SFWF-all.nc'
fn{2} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/LP/SHF-all.nc'
fn{3}='/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/LP/SST-combined.nc'
fn{5}= '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/LP/QFLUX-monmean-combined.nc'
savename='/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholsonCarranza_public/LP_bflux_uncertainty.mat'
% fn0=strcat('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/',...
%     'g.e20b07.2000_DATM%NYF_SLND_CICE_POP2_DROF%NYF_SGLC_SWAV.T62_t13.hybrid.016.pop.h.nday1.0012-12-01.nc');
% path='/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/archive/g.e20b07.2000_DATM%NYF_SLND_CICE_POP2_DROF%NYF_SGLC_SWAV.T62_t13.hybrid.016/ocn/hist/'
% fnstart='g.e20b07.2000_DATM%NYF_SLND_CICE_POP2_DROF%NYF_SGLC_SWAV.T62_t13.hybrid.016.pop.h.00'
% fn{1} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/SFWF-all.nc'
% fn{3}='/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/SST-combined.nc'
% fn{2} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/SHF-all.nc'
% fn{5}= '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/QFLUX-monmean-combined.nc'
% savename='/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholsonCarranza_public/CTL_bflux_uncertainty.mat'

yrs={'08','09','10','11','12'}
mos={'01','02','03','04','05','06','07','08','09','10','11','12'}

skct=1
epsilon1 = 0
cblabel1='m^2/s^3'
varnameout='B'
exportname='fig_4_bflux.png'
caxis1=[-1e-7 1e-7]
caxis2=[0 5e-7]

% code assumes even nt in Fourier transform/plotting:
nt=365*5
tlat = ncread(fn0,'TLAT');
tlon = ncread(fn0,'TLONG');
tlat=tlat(1:skct:end,1:skct:end);
tlon=tlon(1:skct:end,1:skct:end);
tlatg=repmat(tlat,[1 1 60]);
tarea=ncread(fn0,'TAREA')./(1e10);
tarea=tarea(1:skct:end,1:skct:end);
l_h = ncread(fn0,'latent_heat_fusion')./1e4;
g=9.81;
rhosw=ncread(fn0,'rho_sw').*1000;
rhofw=ncread(fn0,'rho_fw').*1000;
cp=ncread(fn0,'cp_sw')./1e4;
saltconst=ncread(fn0,'ocn_ref_salinity');
seaicesalt=ncread(fn0,'sea_ice_salinity');      
QFLUXold=reshape(extract_fieldfn(fn{5},'QFLUX',[0 0 0],[3600 2400 60],[skct skct 1]),[3600/skct 2400/skct 12 5]);

ik=1
var1 = 'SFWF'
var2 = 'SHF'
SFWFscaled = extract_fieldfn(fn{2.*ik-1},var1,[0 0 0],[3600 2400 nt],[skct skct 1]).*...
    saltconst.*g./rhofw;
SHFscaled = extract_fieldfn(fn{2.*ik},var2,[0 0 0],[3600 2400 nt],[skct skct 1]).*...
    g./rhosw./cp;
SFWFscaled(abs(SFWFscaled)>1e8)=nan;
SHFscaled(abs(SHFscaled)>1e8)=nan;
SSTdaily=extract_fieldfn(fn{3},'SST',[0 0 0],[3600 2400 nt],[skct skct 1]);
% fix grid
tlon(tlon>500)=nan;
tlat(tlat>500)=nan;

wts=[31 28 31 30 31 30 31 31 30 31 30 31]';
wtsg=repmat(wts,[1 size(SHFscaled,1) size(SHFscaled,2)]);
wtsg=permute(wtsg,[2 3 1]);
wtsg60=repmat(wtsg,[1 1 5]);

for yrct=1:5
    for moct = 1:12
        SFWFscmo(:,:,moct,yrct)=nansum(SFWFscaled(:,:,...
            (365*(yrct-1)+sum(wts(1:(moct)-1))+1):(365*(yrct-1)+sum(wts(1:(moct))))),3)./nansum(wts(moct));
        SHFscmo(:,:,moct,yrct)=nansum(SHFscaled(:,:,...
            (365*(yrct-1)+sum(wts(1:(moct)-1))+1):(365*(yrct-1)+sum(wts(1:(moct))))),3)./nansum(wts(moct));
    end
end
%clear SHFscaled SFWFscaled


SFWFscmo(tlatg>500)=nan;
SHFscmo(tlatg>500)=nan;

SST=zeros(size(SHFscmo));
SSS=zeros(size(SHFscmo));

for yit=1:5
    yit
    parfor moit = 1:12
    var1 = 'SALT'
    var2 = 'TEMP'
    SALT = [];
    TEMP=[];
    SALT=ncread(strcat(path,fnstart,yrs{yit},'-',mos{moit},'.nc'),var1);
    TEMP=ncread(strcat(path,fnstart,yrs{yit},'-',mos{moit},'.nc'),var2);
    SSS(:,:,moit,yit) = SALT(1:skct:end,1:skct:end,1);
    SST(:,:,moit,yit) = TEMP(1:skct:end,1:skct:end,1);
    end
end

SST=reshape(SST,[size(SST,1) size(SST,2) 60]);
SSS=reshape(SSS,[size(SST,1) size(SST,2) 60]);
wts=[31 28 31 30 31 30 31 31 30 31 30 31];
wts5yr=[wts wts wts wts wts];
wts5yr = wts5yr';
days5yr=cumsum(wts5yr)-wts5yr./2;
daystgt = 0.5:1:(365*5);
SSSdaily=permute(reshape(interp1(days5yr,reshape(permute(SSS,[3 1 2]),[60 400*600]),daystgt,'linear','extrap'),[1825 600 400]),[2 3 1]);

SST(tlatg>500)=nan;
SSS(tlatg>500)=nan;
[PDEN,DRHODT,DRHODS] = mjwfstate(0,SST,SSS);
[PDENdaily,DRHODTdaily,DRHODSdaily] = mjwfstate(0,SSTdaily,SSSdaily);

ALPHA=-DRHODT./rhosw;
ALPHAdaily=-DRHODTdaily./rhosw;
BETA=DRHODS./rhosw;
BETAdaily=DRHODSdaily./rhosw;
clear PDENdaily DRHODTdaily DRHODSdaily SSTdaily SSSdaily

%PDENmean=nansum(wtsg60.*PDEN,3)./...
%    nansum(wtsg60,3);
ALPHAmean=repmat(nansum(wtsg60.*ALPHA,3)./...
    nansum(wtsg60,3),[1 1 60]);
BETAmean=repmat(nansum(wtsg60.*BETA,3)./...
    nansum(wtsg60,3),[1 1 60]);

SHFQFLUXsc= QFLUXold.*g./rhosw./cp;
SFWFQFLUXsc=QFLUXold./l_h.*(saltconst-seaicesalt).*g./rhofw;

save(savename,'-v7.3')
SFWFscmo=reshape(SFWFscmo,[size(SFWFscmo,1) size(SFWFscmo,2) size(SFWFscmo,3)*size(SFWFscmo,4)]);
SHFscmo=reshape(SHFscmo,[size(SHFscmo,1) size(SHFscmo,2) size(SHFscmo,3)*size(SHFscmo,4)]);
SFWFQFLUXsc=reshape(SFWFQFLUXsc,[size(SFWFQFLUXsc,1) size(SFWFQFLUXsc,2) size(SFWFQFLUXsc,3)*size(SFWFQFLUXsc,4)]);
SHFQFLUXsc=reshape(SHFQFLUXsc,[size(SHFQFLUXsc,1) size(SHFQFLUXsc,2) size(SHFQFLUXsc,3)*size(SHFQFLUXsc,4)]);

else
load(savename);
SFWFscmo=reshape(SFWFscmo,[size(SFWFscmo,1) size(SFWFscmo,2) size(SFWFscmo,3)*size(SFWFscmo,4)]);
SHFscmo=reshape(SHFscmo,[size(SHFscmo,1) size(SHFscmo,2) size(SHFscmo,3)*size(SHFscmo,4)]);
SFWFQFLUXsc=reshape(SFWFQFLUXsc,[size(SFWFQFLUXsc,1) size(SFWFQFLUXsc,2) size(SFWFQFLUXsc,3)*size(SFWFQFLUXsc,4)]);
SHFQFLUXsc=reshape(SHFQFLUXsc,[size(SHFQFLUXsc,1) size(SHFQFLUXsc,2) size(SHFQFLUXsc,3)*size(SHFQFLUXsc,4)]);

end
toc;
SFWFdailynoQF=(BETAdaily.*SFWFscaled);
SHFdailynoQF=(ALPHAdaily.*SHFscaled);
fullSFWF=(BETA.*(SFWFscmo+SFWFQFLUXsc));
fullSHF=(ALPHA.*(SHFscmo+SHFQFLUXsc));
meanbetaSFWF=(BETAmean.*(SFWFscmo+SFWFQFLUXsc));
meanalphaSHF=(ALPHAmean.*(SHFscmo+SHFQFLUXsc));
meanbetaSFWFnoQF=(BETAmean.*(SFWFscmo));
meanalphaSHFnoQF=(ALPHAmean.*(SHFscmo));

%%
xidx=[318 197 281 382 33];
yidx=[100 560 548 463 117];
figure;
for ik = 1:3
    subplot(3,4,1+4*(ik-1)),...
        plot(squeeze(ALPHA(yidx(ik),xidx(ik),:)));
    xlabel('Months')
    ylabel('\alpha')
    if ik == 1
        title('Lab. Sea')
    elseif ik == 2
        title('Eastern Equatorial Pacific');
    elseif ik == 3
        title('Eastern subpolar N. Pacific');
    elseif ik == 4
        title('Arctic/Beaufort Gyre');
    elseif ik == 5
        title('Ross Sea');
    end
    subplot(3,4,2+4*(ik-1)),...
        plot(squeeze(BETA(yidx(ik),xidx(ik),:)));
    xlabel('Months')
    ylabel('\beta')
    subplot(3,4,3+4*(ik-1)),...
        plot(squeeze(SST(yidx(ik),xidx(ik),:)));
    xlabel('Months')
    ylabel('SST')
    subplot(3,4,4+4*(ik-1)),...
        plot(squeeze(SSS(yidx(ik),xidx(ik),:)));
    xlabel('Months')
    ylabel('SSS')
end

figure
for ik = 1:3
    subplot(3,1,ik),...
        plot(squeeze(fullSFWF(yidx(ik),xidx(ik),:) + fullSHF(yidx(ik),xidx(ik),:)),'k');
    hold on
    plot(squeeze(fullSFWF(yidx(ik),xidx(ik),:)),'b--','linewidth',2);
    plot(squeeze(fullSHF(yidx(ik),xidx(ik),:)),'r--','linewidth',2);
    plot(squeeze(meanbetaSFWF(yidx(ik),xidx(ik),:)),'b:','linewidth',2);
    plot(squeeze(meanalphaSHF(yidx(ik),xidx(ik),:)),'r:','linewidth',2);
    plot(squeeze(meanbetaSFWFnoQF(yidx(ik),xidx(ik),:)),'cyan--','linewidth',2);
    plot(squeeze(meanalphaSHFnoQF(yidx(ik),xidx(ik),:)),'magenta--','linewidth',2);
    if ik == 1
        title('Lab. Sea')
    elseif ik == 2
        title('Eastern Equatorial Pacific');
    elseif ik == 3
        title('Eastern subpolar N. Pacific');
    elseif ik == 4
        title('Arctic/Beaufort Gyre');
    elseif ik == 5
        title('Ross Sea');
    end
end
%%


         
