%% bflux uncertainty
restoredefaultpath
clear all
close all;
%pause
addpath ./export_fig-master
addpath ./utility
addpath ./cmocean_v1

fn0=strcat('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/',...
   'g.e20b07.2000_DATM%NYF_SLND_CICE_POP2_DROF%NYF_SGLC_SWAV.T62_t13.hybrid.016.pop.h.nday1.0012-12-01.nc');
path='/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/archive/g.e20b07.2000_DATM%NYF_SLND_CICE_POP2_DROF%NYF_SGLC_SWAV.T62_t13.hybrid.036/ocn/hist/'
fnstart='g.e20b07.2000_DATM%NYF_SLND_CICE_POP2_DROF%NYF_SGLC_SWAV.T62_t13.hybrid.036.pop.h.00'
fn{1} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/LP/SFWF-all.nc'
fn{2} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/LP/SHF-all.nc'
fn{3}='/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/LP/SST-combined.nc'
fn{5}= '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/LP/QFLUX-monmean-combined.nc'
savename='/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/LP/BFLUX-all.nc'
% fn0=strcat('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/',...
%     'g.e20b07.2000_DATM%NYF_SLND_CICE_POP2_DROF%NYF_SGLC_SWAV.T62_t13.hybrid.016.pop.h.nday1.0012-12-01.nc');
% path='/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/archive/g.e20b07.2000_DATM%NYF_SLND_CICE_POP2_DROF%NYF_SGLC_SWAV.T62_t13.hybrid.016/ocn/hist/'
% fnstart='g.e20b07.2000_DATM%NYF_SLND_CICE_POP2_DROF%NYF_SGLC_SWAV.T62_t13.hybrid.016.pop.h.00'
% fn{1} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/SFWF-all.nc'
% fn{3}='/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/SST-combined.nc'
% fn{2} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/SHF-all.nc'
% fn{5}= '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/QFLUX-monmean-combined.nc'
%savename='/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/BFLUX-all.nc'

yrs={'08','09','10','11','12'}
mos={'01','02','03','04','05','06','07','08','09','10','11','12'}

skct=1
epsilon1 = 0
cblabel1='m^2/s^3'
varnameout='B'

% code assumes even nt in Fourier transform/plotting:
tlat = ncread(fn0,'TLAT');
tlon = ncread(fn0,'TLONG');
tlat=tlat(1:skct:end,1:skct:end);
tlon=tlon(1:skct:end,1:skct:end);
ncwrite(savename,'TLAT',tlat);
ncwrite(savename,'TLONG',tlon);
l_h = ncread(fn0,'latent_heat_fusion')./1e4;
g=9.81;
rhosw=ncread(fn0,'rho_sw').*1000;
rhofw=ncread(fn0,'rho_fw').*1000;
cp=ncread(fn0,'cp_sw')./1e4;
saltconst=ncread(fn0,'ocn_ref_salinity');
seaicesalt=ncread(fn0,'sea_ice_salinity');      
nt=365
for sttx = 0:365:1820
sttx
SSTdaily=extract_fieldfn(fn{3},'SST',[0 0 sttx],[3600 2400 nt],[skct skct 1]);
SSTdaily=single(SSTdaily);
% fix grid
SSS=zeros([3600/skct 2400/skct 12]);
    yit=round(sttx/365)+1
    for moit = 1:12
    var1 = 'SALT'
    SALT = [];
    SALT=ncread(strcat(path,fnstart,yrs{yit},'-',mos{moit},'.nc'),var1);
    SSS(:,:,moit) = SALT(1:skct:end,1:skct:end,1);
    end
clear SALT

wts=[31 28 31 30 31 30 31 31 30 31 30 31]';
days=cumsum(wts)-wts./2;
daystgt = 0.5:1:(365);
disp('calculating SSS daily...')
SSSdaily=permute(reshape(interp1(days,reshape(permute(SSS,[3 1 2]),[12 3600*2400/(skct*skct)]),daystgt','linear','extrap'),[nt 3600/skct 2400/skct]),[2 3 1]);
clear SSS
SSSdaily=single(SSSdaily);
[PDENdaily,DRHODTdaily,DRHODSdaily] = mjwfstate(0,SSTdaily,SSSdaily);
clear PDENdaily SSTdaily SSSdaily
ALPHAdaily=-DRHODTdaily./rhosw;
clear DRHODTdaily
BETAdaily=DRHODSdaily./rhosw;
clear DRHODSdaily
disp('writing A/B...')
ncwrite(savename,'ALPHA',ALPHAdaily,[1 1 (1+sttx)]);
ncwrite(savename,'BETA',BETAdaily,[1 1 (1+sttx)]);
ik=1
var1 = 'SFWF'
var2 = 'SHF'
SFWFscaled = BETAdaily.*extract_fieldfn(fn{2.*ik-1},var1,[0 0 sttx],[3600 2400 nt],[skct skct 1]).*...
    saltconst.*g./rhofw;
SHFscaled = ALPHAdaily.*extract_fieldfn(fn{2.*ik},var2,[0 0 sttx],[3600 2400 nt],[skct skct 1]).*...
    g./rhosw./cp;
SFWFscaled(abs(SFWFscaled)>1e4)=nan;
SHFscaled(abs(SHFscaled)>1e4)=nan;
disp('writing FLUX scaled...')
ncwrite(savename,'SHFscaled',SHFscaled,[1 1 (1+sttx)]);
clear SHFscaled
ncwrite(savename,'SFWFscaled',SFWFscaled,[1 1 (1+sttx)]);
clear SFWFscaled
end
