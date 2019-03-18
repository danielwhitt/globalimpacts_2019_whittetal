restoredefaultpath
clear all
%close all;
figure('Position',[20 20 1000 600]);
%pause
addpath ./export_fig-master
addpath ./utility
gladeflag =1
if gladeflag == 1

fn0=strcat('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/',...
    'g.e20b07.2000_DATM%NYF_SLND_CICE_POP2_DROF%NYF_SGLC_SWAV.T62_t13.hybrid.016.pop.h.nday1.0012-12-01.nc');

fn{1} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/TAUmag-all.nc'
fn{2} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/LP/TAUmag-all.nc'

cblabel1='log_{10} N/m^2'
varoutname='|\tau|'
caxis1=log10([1e-2.*.3 .3])
epsilon1=.01
exportname='fig_2_stress.png'
exportname2='fig_fracvartau.png'
exportname3='fig_phasevartau.png'
ylabel2='N^2/(m^4 cyc/day)'
skct=6
% code assumes even nt in Fourier transform/plotting:
nt=1825
nfreq=364
tlat=extract_fieldfn(fn0,'ULAT',[0 0],[3600 2400],[skct skct]);
tlon=extract_fieldfn(fn0,'ULONG',[0 0],[3600 2400],[skct skct]);
tlon(tlon>500)=nan;
tlat(tlat>500)=nan;
tlatg=repmat(tlat,[1 1 nfreq]);
tarea=extract_fieldfn(fn0,'UAREA',[0 0],[3600 2400],[skct skct])./1e10;
tareag=repmat(tarea,[1 1 nfreq]);

for ik = 1:2
    var = 'TAUMAG'
    if ik == 1
        varout = 1e-1.*extract_fieldfn(fn{ik},var,[0 0 0],[3600 2400 nt],[skct skct 1]);
        time1=0.5:(nt-0.5);
    else
        varout = 1e-1.*extract_fieldfn(fn{ik},var,[0 0 0],[3600 2400 nt],[skct skct 1]);
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
save('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholsonCarranza_public/fig_2_stress_CI.mat','-v7.3')
else
load('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholsonCarranza_public/fig_2_stress_CI.mat')
end
%%
create_stats_2;
%% plotting
plot_fig_2;

%% more stats
create_stats_2b;

