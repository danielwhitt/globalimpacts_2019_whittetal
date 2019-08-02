restoredefaultpath
clear all
%close all
fricvel = 1
addpath ./utility
addpath ./cmocean_v1
addpath ./export_fig-master/

fn0= strcat('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/',...
    'g.e20b07.2000_DATM%NYF_SLND_CICE_POP2_DROF%NYF_SGLC_SWAV.T62_t13.hybrid.016.pop.h.nday1.0012-12-01.nc')
%CTL
fn{1} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/TAUmag-all.nc'
% LP:
fn{2} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/LP/TAUmag-all.nc'

skct=6
if fricvel == 1
    % friction velocity
    epsilon1=1e-4
    epsilon1=0
    cblabel1='m/s'
    varnameout='u_*'
    exportname='fig_1_fricvel.png'
    caxis1=[0 .018]
    caxis2=[0 .018]
else
    epsilon1=1e-2
    epsilon1 = 0
    cblabel1='N/m^2'
    varnameout='|\tau|'
    exportname='fig_1_stress.png'
    caxis1=[0 .3]
    caxis2=[0 .3]
end



% grid variables
tlat = ncread(fn0,'ULAT');
tlon = ncread(fn0,'ULONG');
tlat=tlat(1:skct:end,1:skct:end);
tlon=tlon(1:skct:end,1:skct:end);
tlatg=repmat(tlat,[1 1 5*365]);
tarea=ncread(fn0,'UAREA')./(1e10);
tarea=tarea(1:skct:end,1:skct:end);

%% extract variables
for ik = 1:2
    var = 'TAUMAG'
    if ik == 1
        varout = 1e-1.*extract_fieldfn(fn{ik},var,[0 0 0],[3600 2400 5*365],[skct skct 1]);
    else
        varout = 1e-1.*extract_fieldfn(fn{ik},var,[0 0 0],[3600 2400 5*365],[skct skct 1]);
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
%%
if fricvel==1
    % friction velocity
    varoutold=sqrt(varoutold./1026);
    varout=sqrt(varout./1026);
end
tlatg=repmat(tlat,[1 1 12]);
save('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholsonCarranza_public/fig_1_fricvel_CI.mat','-v7.3')
%%
create_monthly_means;
%%
create_stats_1;
%% plot
plot_fig_1;
%% Some area stats
create_stats_1b;
save('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholsonCarranza_public/fig_1_fricvel_fig_9_load.mat')



