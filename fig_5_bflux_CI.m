restoredefaultpath
clear all
%close all;
addpath ./export_fig-master
addpath ./utility

fn0=strcat('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/',...
    'g.e20b07.2000_DATM%NYF_SLND_CICE_POP2_DROF%NYF_SGLC_SWAV.T62_t13.hybrid.016.pop.h.nday1.0012-12-01.nc');

fn{1} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/SFWF-all.nc'
fn{2} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/SHF-all.nc'
fn{3} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/LP/SFWF-all.nc'
fn{4} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/LP/SHF-all.nc'
%fn{5}='/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/QFLUX-monmean-combined.nc'
%fn{6}='/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/LP/QFLUX-monmean-combined.nc'
cblabel1='log_{10} m^2/s^3'
varoutname='B'
caxis1=log10([1e-2.*1.76E-7 (1.76E-7)])
epsilon1=1e-9
epsilon1=0
exportname='fig_5_bflux.png'
exportname2='fig_fracvarB.png'
exportname3='fig_phasevarB.png'
ylabel2='m^4/(s^6 cyc/day)'

% this skct has to be 6 because of the saved alpha/beta parameters
skct=6
nt=1825
nfreq=364
tlat=extract_fieldfn(fn0,'TLAT',[0 0],[3600 2400],[skct skct]);
tlon=extract_fieldfn(fn0,'TLONG',[0 0],[3600 2400],[skct skct]);
tlon(tlon>500)=nan;
tlat(tlat>500)=nan;
tarea=extract_fieldfn(fn0,'TAREA',[0 0],[3600 2400],[skct skct])./1e10;
tareag=repmat(tarea,[1 1 nfreq]);
tlatg=repmat(tlat,[1 1 nt]);

l_h = ncread(fn0,'latent_heat_fusion')./1e4;
g=9.81;
rhosw=ncread(fn0,'rho_sw').*1000;
rhofw=ncread(fn0,'rho_fw').*1000;
cp=ncread(fn0,'cp_sw')./1e4;
saltconst=ncread(fn0,'ocn_ref_salinity');
seaicesalt=ncread(fn0,'sea_ice_salinity');      

%% extract ALPHA/BETA
wts=[31 28 31 30 31 30 31 31 30 31 30 31];
wts5yr=[wts wts wts wts wts];
wts5yr = wts5yr';
days5yr=cumsum(wts5yr)-wts5yr./2;
daystgt = 0.5:1:(365*5);
dailyalpha=0
if dailyalpha==1
    load('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholsonCarranza_public/LP_bflux_uncertainty.mat','ALPHAdaily','BETAdaily');
    ALPHALPdaily=-ALPHAdaily;
    BETALPdaily=BETAdaily;
    clear ALPHAdaily BETAdaily
    load('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholsonCarranza_public/CTL_bflux_uncertainty.mat','ALPHAdaily','BETAdaily');
    ALPHACTLdaily=-ALPHAdaily;
    BETACTLdaily=BETAdaily;
    clear ALPHAdaily BETAdaily
else
    load('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholsonCarranza_public/LP_bflux_uncertainty.mat','ALPHA','BETA');
    ALPHALP=-ALPHA;
    ALPHALPdaily=permute(reshape(interp1(days5yr,reshape(permute(ALPHALP,[3 1 2]),[60 400*600]),daystgt,'linear','extrap'),[1825 600 400]),[2 3 1]);
    BETALP=BETA;
    BETALPdaily=permute(reshape(interp1(days5yr,reshape(permute(BETALP,[3 1 2]),[60 400*600]),daystgt,'linear','extrap'),[1825 600 400]),[2 3 1]);
    clear ALPHA BETA
    load('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholsonCarranza_public/CTL_bflux_uncertainty.mat','ALPHA','BETA');
    ALPHACTL=-ALPHA;
    ALPHACTLdaily=permute(reshape(interp1(days5yr,reshape(permute(ALPHACTL,[3 1 2]),[60 400*600]),daystgt,'linear','extrap'),[1825 600 400]),[2 3 1]);
    BETACTL=BETA;
    BETACTLdaily=permute(reshape(interp1(days5yr,reshape(permute(BETACTL,[3 1 2]),[60 400*600]),daystgt,'linear','extrap'),[1825 600 400]),[2 3 1]);
end
alphaboth = zeros(size(ALPHALPdaily,1),size(ALPHALPdaily,2),size(ALPHALPdaily,3),2);
alphaboth(:,:,:,2)=ALPHALPdaily;
alphaboth(:,:,:,1)=ALPHACTLdaily;
betaboth = zeros(size(BETALPdaily,1),size(BETALPdaily,2),size(BETALPdaily,3),2);
betaboth(:,:,:,2)=BETALPdaily;
betaboth(:,:,:,1)=BETACTLdaily;
clear ALPHACTLdaily BETACTLdaily ALPHALPdaily BETALPdaily
%%


for ik = 1:2
    
    var1 = 'SFWF'
    var2 = 'SHF'
    varout = extract_fieldfn(fn{2.*ik-1},var1,[0 0 0],[3600 2400 nt],[skct skct 1]).*...
        saltconst.*g.*(squeeze(betaboth(:,:,:,ik)))./rhofw-...
        extract_fieldfn(fn{2.*ik},var2,[0 0 0],[3600 2400 nt],[skct skct 1]).*...
        g.*(squeeze(alphaboth(:,:,:,ik)))./rhosw./cp;

    time1=0.5:(nt-0.5);
    varout(abs(varout)>1)=nan;
    % fix grid
    varout(tlatg>500)=nan;

    dt=time1(2)-time1(1);
    Lt = time1(end)-time1(1)+1;
    Ltfreq=time1(nfreq)-time1(1)+1;
    varoutPS=zeros(size(varout,1),size(varout,2),nfreq);
    varoutPhS=zeros(size(varout,1),size(varout,2),nfreq);
    numitspec=round(nt./nfreq)
    for it=1:numitspec
        % then compute the ffts
        clear varhat
        [varhat,freq] = fftfuncore(varout(:,:,(it-1)*365+(1:nfreq)),1./dt);
        freq = freq.*2.*pi;
        dfreq=freq(2)-freq(1);
        varoutPS(:,:,:,it)=(varhat.*conj(varhat))./(2.*pi.*Ltfreq);
        varoutPhS(:,:,:,it)=angle(varhat);
    end
    if ik == 1
        varoutold=varout;
        varoutoldPS=varoutPS;
        varoutoldPhS=varoutPhS;
        clear varhat varoutPS varout varoutPhS
    end
end
tlatg=repmat(tlat,[1 1 nfreq]);
clear alphaboth betaboth
save('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholsonCarranza_public/fig_5_bflux_CI.mat','-v7.3')

%% calculate some stats
create_stats_2;

%% plotting
plot_fig_2;

%% more stats
create_stats_2b;


