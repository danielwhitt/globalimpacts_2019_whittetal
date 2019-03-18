restoredefaultpath
clear all
close all;
%pause
addpath ./export_fig-master
addpath ./utility
addpath ./cmocean_v1


fn0=strcat('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/',...
    'g.e20b07.2000_DATM%NYF_SLND_CICE_POP2_DROF%NYF_SGLC_SWAV.T62_t13.hybrid.016.pop.h.nday1.0012-12-01.nc');

fn{1} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/SFWF-all.nc'
fn{2} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/SHF-all.nc'
fn{3} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/LP/SFWF-all.nc'
fn{4} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/LP/SHF-all.nc'
fn{5}='/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/QFLUX-monmean-combined.nc'
fn{6}='/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/LP/QFLUX-monmean-combined.nc'

skct=6
epsilon1=1e-9
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
tlatg=repmat(tlat,[1 1 12]);
tarea=ncread(fn0,'TAREA')./(1e10);
tarea=tarea(1:skct:end,1:skct:end);
load('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/alphabetamean.mat','alphaCTL','betaCTL','alphaLP','betaLP');
l_h = ncread(fn0,'latent_heat_fusion')./1e4;
g=9.81;
rhosw=ncread(fn0,'rho_sw').*1000;
rhofw=ncread(fn0,'rho_fw').*1000;
cp=ncread(fn0,'cp_sw')./1e4;
saltconst=ncread(fn0,'ocn_ref_salinity');
seaicesalt=ncread(fn0,'sea_ice_salinity');      
alphaboth(:,:,2)=alphaLP;
betaboth(:,:,2)=betaLP;
alphaboth(:,:,1)=alphaCTL;
betaboth(:,:,1)=betaCTL;
QFLUXold=reshape(extract_fieldfn(fn{5},'QFLUX',[0 0 0],[3600 2400 60],[skct skct 1]),[600 400 12 5]);
QFLUX=reshape(extract_fieldfn(fn{6},'QFLUX',[0 0 0],[3600 2400 60],[skct skct 1]),[600 400 12 5]);

for ik = 1:2
    var1 = 'SFWF'
    var2 = 'SHF'
    varout = extract_fieldfn(fn{2.*ik-1},var1,[0 0 0],[3600 2400 nt],[skct skct 1]).*...
        repmat(saltconst.*g.*(squeeze(betaboth(:,:,ik)))./rhofw,[1 1 nt])-...
        extract_fieldfn(fn{2.*ik},var2,[0 0 0],[3600 2400 nt],[skct skct 1]).*...
        repmat(g.*(squeeze(alphaboth(:,:,ik)))./rhosw./cp,[1 1 nt]);
    varout(abs(varout)>1)=nan;
    % fix grid
    tlon(tlon>500)=nan;
    tlat(tlat>500)=nan;
    varout(tlatg>500)=nan;
    if ik == 1
        varoutold=varout;
    end
end



%%
create_monthly_means;
% add qflux at monthly mean stage
varoutmo=varoutmo-QFLUX.*repmat(g.*(squeeze(alphaboth(:,:,2)))./rhosw./cp,[1 1 12 5])+...
                  QFLUX./l_h.*repmat((saltconst-seaicesalt).*g.*(squeeze(betaboth(:,:,2)))./rhofw,[1 1 12 5]);
varoutoldmo=varoutoldmo-QFLUXold.*repmat(g.*(squeeze(alphaboth(:,:,1)))./rhosw./cp,[1 1 12 5])+...
                  QFLUXold./l_h.*repmat((saltconst-seaicesalt).*g.*(squeeze(betaboth(:,:,1)))./rhofw,[1 1 12 5]);
save('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholsonCarranza_public/fig_4_bflux_CI.mat','-v7.3')

%%
create_stats_1;
%% plot
plot_fig_1;
%% Some area stats
create_stats_1b;

