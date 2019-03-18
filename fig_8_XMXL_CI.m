restoredefaultpath
clear all
%close all;
addpath ./export_fig-master
addpath ./utility
addpath ./cmocean_v1/

fn0=strcat('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/',...
    'g.e20b07.2000_DATM%NYF_SLND_CICE_POP2_DROF%NYF_SGLC_SWAV.T62_t13.hybrid.016.pop.h.nday1.0012-12-01.nc');

fn{1} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/XMXL-all.nc'
fn{2} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/LP/XMXL-all.nc'
cblabel1='log_{10} m'
varoutname='MLD'
caxis1=log10([1e-3.*600./2.8 (600./2.8)])
epsilon1=1
exportname='fig_8_MLD.png'
exportname2='fig_fracvarMLD.png'
exportname3='fig_phasevarMLD.png'
skct=6
ylabel2='m^2/(cyc/day)'
nt=1825
nfreq=364 % only does an even number of frequencies
tlat=extract_fieldfn(fn0,'TLAT',[0 0],[3600 2400],[skct skct]);
tlon=extract_fieldfn(fn0,'TLONG',[0 0],[3600 2400],[skct skct]);
tlon(tlon>500)=nan;
tlat(tlat>500)=nan;
tlatg=repmat(tlat,[1 1 nfreq]);
tarea=extract_fieldfn(fn0,'TAREA',[0 0],[3600 2400],[skct skct])./1e10;
tareag=repmat(tarea,[1 1 nfreq]);
maskhotspots=extract_fieldfn('regional_mask_all1.nc','mask_hotspots',[0 0],[3600 2400],[skct skct]);


for ik = 1:2
    var = 'XMXL'
    if ik == 1
        varout = 1e-2.*extract_fieldfn(fn{ik},var,[0 0 0],[3600 2400 nt],[skct skct 1]);
        time1=0.5:(nt-0.5);
    else
        varout = 1e-2.*extract_fieldfn(fn{ik},var,[0 0 0],[3600 2400 nt],[skct skct 1]);
        time1=0.5:(nt-0.5);
    end
    varout(varout>1e10)=nan;
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
save('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholsonCarranza_public/fig_8_XMXL_CI.mat','-v7.3')

create_stats_2;

%% plotting
plot_fig_2;
%% more stats
create_stats_2b;
plot_fig_hotspotspectra;
 
