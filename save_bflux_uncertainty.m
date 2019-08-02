%% bflux uncertainty
restoredefaultpath
clear all
close all;
%pause
addpath ./export_fig-master
addpath ./utility
addpath ./cmocean_v1
tic
loaddataflag =0
if loaddataflag == 1
fn0=strcat('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/',...
    'g.e20b07.2000_DATM%NYF_SLND_CICE_POP2_DROF%NYF_SGLC_SWAV.T62_t13.hybrid.016.pop.h.nday1.0012-12-01.nc');

path='/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/archive/g.e20b07.2000_DATM%NYF_SLND_CICE_POP2_DROF%NYF_SGLC_SWAV.T62_t13.hybrid.016/ocn/hist/'
fnstart='g.e20b07.2000_DATM%NYF_SLND_CICE_POP2_DROF%NYF_SGLC_SWAV.T62_t13.hybrid.016.pop.h.00'
yrs={'08','09','10','11','12'}
mos={'01','02','03','04','05','06','07','08','09','10','11','12'}

fn{1} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/SFWF-all.nc'
fn{2} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/SHF-all.nc'
%fn{3} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/LP/SFWF-all.nc'
%fn{4} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/LP/SHF-all.nc'
fn{5}= '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/QFLUX-monmean-combined.nc'
%fn{6}= '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/LP/QFLUX-monmean-combined.nc'

load('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/alphabetamean.mat','alphaCTL','betaCTL','alphaLP','betaLP');

skct=6
epsilon1=1e-9
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
SHFscaled = -extract_fieldfn(fn{2.*ik},var2,[0 0 0],[3600 2400 nt],[skct skct 1]).*...
    g./rhosw./cp;
SFWFscaled(abs(SFWFscaled)>1e8)=nan;
SHFscaled(abs(SHFscaled)>1e8)=nan;

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
SST(tlatg>500)=nan;
SSS(tlatg>500)=nan;
[PDEN,DRHODT,DRHODS] = mjwfstate(0,SST,SSS);
ALPHA=-DRHODT./rhosw;
BETA=DRHODS./rhosw;

%PDENmean=nansum(wtsg60.*PDEN,3)./...
%    nansum(wtsg60,3);
ALPHAmean=repmat(nansum(wtsg60.*ALPHA,3)./...
    nansum(wtsg60,3),[1 1 60]);
BETAmean=repmat(nansum(wtsg60.*BETA,3)./...
    nansum(wtsg60,3),[1 1 60]);

SHFQFLUXsc=-QFLUXold.*g./rhosw./cp;
SFWFQFLUXsc=QFLUXold./l_h.*(saltconst-seaicesalt).*g./rhofw;

save('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholsonCarranza_public/bflux_uncertainty.mat','-v7.3')
else
load('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholsonCarranza_public/bflux_uncertainty.mat','-v7.3')
end
   toc           
