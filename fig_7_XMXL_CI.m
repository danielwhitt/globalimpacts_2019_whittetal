restoredefaultpath
clear all
%close all
addpath ~/ML_diagnostics/matlabfunctions/m_map
fn0='/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/g.e20b07.2000_DATM%NYF_SLND_CICE_POP2_DROF%NYF_SGLC_SWAV.T62_t13.hybrid.016.pop.h.nday1.0012-12-01.nc'
fn{1} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/XMXL-all.nc'
fn{2} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/LP/XMXL-all.nc'
%pause
addpath ./export_fig-master
addpath ./utility
epsilon1 = 0
cblabel1='m'
varnameout='MLD'
exportname='fig_7_MLD.png'
caxis1=[0 200]
caxis2=[0 600]
skct=6

tlat = ncread(fn0,'TLAT');
tlon = ncread(fn0,'TLONG');
tlat=tlat(1:skct:end,1:skct:end);
tlon=tlon(1:skct:end,1:skct:end);
tlon(tlon>500)=nan;
tlat(tlat>500)=nan;
tlatg=repmat(tlat,[1 1 5*365]);
tarea=ncread(fn0,'TAREA')./(1e10);
tarea=tarea(1:skct:end,1:skct:end);
tareag=repmat(tarea,[1 1 12]);


%% extract variables
for ik = 1:2
    var = 'XMXL'
    if ik == 1
        varout = 1e-2.*extract_fieldfn(fn{ik},var,[0 0 0],[3600 2400 5*365],[skct skct 1]);
        time1=0.5:1:(365*5-0.5);
    else
        varout = 1e-2.*extract_fieldfn(fn{ik},var,[0 0 0],[3600 2400 5*365],[skct skct 1]);
        time1=0.5:1:(365*5-0.5);
    end
    varout(varout>1e10)=nan;
    % fix grid
    tlon(tlon>500)=nan;
    tlat(tlat>500)=nan;
    varout(tlatg>500)=nan;
    if ik == 1
        varoutold=varout;
        clear varout varnow
    end
end
tlatg=repmat(tlat,[1 1 12]);
save('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholsonCarranza_public/fig_7_XMXL_CI.mat','-v7.3')

%%
create_monthly_means;
%%
create_stats_1;
%% plot
plot_fig_1;
%% Some area stats
create_stats_1b;
save('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholsonCarranza_public/fig_7_XMXL_fig_9_load.mat')

