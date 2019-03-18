restoredefaultpath
clear all
%close all;
figure('Position',[20 20 1000 600]);
%pause
addpath ./export_fig-master
addpath ./utility

fn0=strcat('/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/',...
    'g.e20b07.2000_DATM%NYF_SLND_CICE_POP2_DROF%NYF_SGLC_SWAV.T62_t13.hybrid.016.pop.h.nday1.0012-12-01.nc');

fn{1} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/TAUX-all.nc'
fn{2} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/LP/TAUX-all.nc'
fn{3} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/TAUY-all.nc'
fn{4} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/LP/TAUY-all.nc'
fn{5} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/CTL/TAUmag-all.nc'
fn{6} = '/glade/p/cgd/oce/people/dwhitt/nsfsubmeso/WhittNicholson_data/LP/TAUmag-all.nc'


skct=12
% code assumes even nt in Fourier transform/plotting:
nt=365
tlat=extract_fieldfn(fn0,'ULAT',[0 0],[3600 2400],[skct skct]);
tlon=extract_fieldfn(fn0,'ULONG',[0 0],[3600 2400],[skct skct]);
ANGLE=extract_fieldfn(fn0,'ANGLE',[0 0],[3600 2400],[skct skct]);
tlon(tlon>500)=nan;
tlat(tlat>500)=nan;
tlatg=repmat(tlat,[1 1 nt]);
tarea=extract_fieldfn(fn0,'UAREA',[0 0],[3600 2400],[skct skct])./1e10;
tareag=repmat(tarea,[1 1 nt]);

for ik = 1:2
    var1 = 'TAUX'
    var2 = 'TAUY'
    var3 = 'TAUMAG'
    varout1 = 1e-1.*extract_fieldfn(fn{ik},var1,[0 0 0],[3600 2400 nt],[skct skct 1]);
    varout2 = 1e-1.*extract_fieldfn(fn{2+ik},var2,[0 0 0],[3600 2400 nt],[skct skct 1]);
    varout3 = 1e-1.*extract_fieldfn(fn{4+ik},var3,[0 0 0],[3600 2400 nt],[skct skct 1]);
    time1=0.5:(nt-0.5);

    varout1(varout1>1e10)=nan;
    varout1(tlatg>500)=nan;
    varout2(varout2>1e10)=nan;
    varout2(tlatg>500)=nan;    
    varout3(varout3>1e10)=nan;
    varout3(tlatg>500)=nan;
    dt=time1(2)-time1(1)
    if ik == 1
        varoutold1=varout1;
        varoutold2=varout2;
        varoutold3=varout3;
        clear varout1 varout2 varout3
    end
end
%%
figure('Position',[20 20 1000 600])
for il = 1:3
    clear idx1 idx2 titlestr titlestr1 titlestr2
    if il == 1
        idx1=52
        idx2=158
        titlestr=strcat('Lab. Sea, lat=',num2str(tlat(idx1,idx2),3),', lon=',num2str(tlon(idx1,idx2),3))
        titlestr1=['(A) \tau_x, ',titlestr];
        titlestr2=['(B) |\tau|, ',titlestr];
    elseif il == 2
        idx1=73
        idx2=101
        titlestr=strcat('Equatorial Atl., lat=',num2str(tlat(idx1,idx2),3),', lon=',num2str(tlon(idx1,idx2),3)) 
        titlestr1=['(C) \tau_x, ',titlestr];
        titlestr2=['(D) |\tau|, ',titlestr];
    elseif il == 3
        idx1=90
        idx2=50
        titlestr=strcat('S. Atl., lat=',num2str(tlat(idx1,idx2),3),', lon=',num2str(tlon(idx1,idx2),3))  
        titlestr1=['(E) \tau_x, ',titlestr];
        titlestr2=['(F) |\tau|, ',titlestr];
    end
subplot(3,2,2*il-1),...
plot(time1,squeeze(varoutold1(idx1,idx2,:)).*cos(-ANGLE(idx1,idx2))-...
squeeze(varoutold2(idx1,idx2,:)).*sin(-ANGLE(idx1,idx2)),'b-','linewidth',1);
hold on
plot(time1,squeeze(varout1(idx1,idx2,:)).*cos(-ANGLE(idx1,idx2))-...
squeeze(varout2(idx1,idx2,:)).*sin(-ANGLE(idx1,idx2)),'r-','linewidth',1);
xlabel('day of year')
%if il == 1
%legend('\tau_x, CTL','\tau_x, LP');
%end
ylabel('N/m^2')
title(titlestr1)
grid on
xlim([0 365])
subplot(3,2,2*il),...
plot(time1,squeeze(varoutold3(idx1,idx2,:)),'b-','linewidth',1);
hold on
plot(time1,squeeze(varout3(idx1,idx2,:)),'r-','linewidth',1);
xlabel('day of year')
if il == 1
legend('CTL','LP')
end
ylabel('N/m^2')
title(titlestr2)
grid on
xlim([0 365])
end
set(gcf,'color','w')

export_fig fig_stress_timeseries_supplement.pdf -painters
