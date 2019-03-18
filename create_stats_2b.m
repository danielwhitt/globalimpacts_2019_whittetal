tarea(isnan(squeeze(varoutold(:,:,1))))=nan;
clear varoutclmo varoutoldclmo
varoutclmo=squeeze(nanmean(sqrt(sum(abs(varoutPS).*dfreq,3)-abs(varoutPS(:,:,nfreq/2+1,:)).*dfreq),4));
varoutoldclmo=squeeze(nanmean(sqrt(sum(abs(varoutoldPS).*dfreq,3)-abs(varoutoldPS(:,:,nfreq/2+1,:)).*dfreq),4));

display('........')
display('Difference, global area avg. std.: ')
CTLglobal=nansum(nansum(tarea.*varoutoldclmo,2),1)./nansum(nansum(tarea,2),1)
LPglobal=nansum(nansum(tarea.*varoutclmo,2),1)./nansum(nansum(tarea,2),1)
ratioout=1-nansum(nansum(tarea.*varoutclmo./varoutoldclmo,2),1)./nansum(nansum(tarea,2),1)
1-LPglobal/CTLglobal

display('<-66:')
tlatgmask=ones(size(tlat));
tlatgmask(tlat>-66)=nan;
CTLglobal=nansum(nansum(tarea.*varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
LPglobal=nansum(nansum(tarea.*varoutclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
ratioglobal=1-nansum(nansum(tarea.*varoutclmo./varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
1-LPglobal/CTLglobal

display('>66:')
tlatgmask=ones(size(tlat));
tlatgmask(tlat<66)=nan;
CTLglobal=nansum(nansum(tarea.*varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
LPglobal=nansum(nansum(tarea.*varoutclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
ratioglobal=1-nansum(nansum(tarea.*varoutclmo./varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
1-LPglobal/CTLglobal


display('abs(lat) < 30:')
tlatgmask=ones(size(tlat));
tlatgmask(abs(tlat)>30)=nan;
CTLglobal=nansum(nansum(tarea.*varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
LPglobal=nansum(nansum(tarea.*varoutclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
ratioglobal=1-nansum(nansum(tarea.*varoutclmo./varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
1-LPglobal/CTLglobal

display('> 30:')
tlatgmask=ones(size(tlat));
tlatgmask(tlat<30)=nan;
CTLglobal=nansum(nansum(tarea.*varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
LPglobal=nansum(nansum(tarea.*varoutclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
ratioglobal=1-nansum(nansum(tarea.*varoutclmo./varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
1-LPglobal/CTLglobal


display('< -30:')
tlatgmask=ones(size(tlat));
tlatgmask(tlat>-30)=nan;
CTLglobal=nansum(nansum(tarea.*varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
LPglobal=nansum(nansum(tarea.*varoutclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
ratioglobal=1-nansum(nansum(tarea.*varoutclmo./varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
1-LPglobal/CTLglobal








clear varoutclmo varoutoldclmo
varoutclmo=squeeze(nanmean(sqrt(sum(2.*abs(varoutPS(:,:,(nfreq/2+8):end,:)).*dfreq,3)),4));
varoutoldclmo=squeeze(nanmean(sqrt(sum(2.*abs(varoutoldPS(:,:,(nfreq/2+8):end,:)).*dfreq,3)),4));

display('........')
display('Difference, global area avg. intraseasonal std.: ')
CTLglobal=nansum(nansum(tarea.*varoutoldclmo,2),1)./nansum(nansum(tarea,2),1)
LPglobal=nansum(nansum(tarea.*varoutclmo,2),1)./nansum(nansum(tarea,2),1)
ratioout=1-nansum(nansum(tarea.*varoutclmo./varoutoldclmo,2),1)./nansum(nansum(tarea,2),1)
1-LPglobal/CTLglobal

display('<-66:')
tlatgmask=ones(size(tlat));
tlatgmask(tlat>-66)=nan;
CTLglobal=nansum(nansum(tarea.*varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
LPglobal=nansum(nansum(tarea.*varoutclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
ratioglobal=1-nansum(nansum(tarea.*varoutclmo./varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
1-LPglobal/CTLglobal

display('>66:')
tlatgmask=ones(size(tlat));
tlatgmask(tlat<66)=nan;
CTLglobal=nansum(nansum(tarea.*varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
LPglobal=nansum(nansum(tarea.*varoutclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
ratioglobal=1-nansum(nansum(tarea.*varoutclmo./varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
1-LPglobal/CTLglobal


display('abs(lat) < 30:')
tlatgmask=ones(size(tlat));
tlatgmask(abs(tlat)>30)=nan;
CTLglobal=nansum(nansum(tarea.*varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
LPglobal=nansum(nansum(tarea.*varoutclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
ratioglobal=1-nansum(nansum(tarea.*varoutclmo./varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
1-LPglobal/CTLglobal

display('> 30:')
tlatgmask=ones(size(tlat));
tlatgmask(tlat<30)=nan;
CTLglobal=nansum(nansum(tarea.*varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
LPglobal=nansum(nansum(tarea.*varoutclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
ratioglobal=1-nansum(nansum(tarea.*varoutclmo./varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
1-LPglobal/CTLglobal


display('< -30:')
tlatgmask=ones(size(tlat));
tlatgmask(tlat>-30)=nan;
CTLglobal=nansum(nansum(tarea.*varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
LPglobal=nansum(nansum(tarea.*varoutclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
ratioglobal=1-nansum(nansum(tarea.*varoutclmo./varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
1-LPglobal/CTLglobal





clear varoutclmo varoutoldclmo
varoutclmo=squeeze(nanmean((sum(abs(varoutPS).*dfreq,3)-abs(varoutPS(:,:,nfreq/2+1,:)).*dfreq),4));
varoutoldclmo=squeeze(nanmean((sum(abs(varoutoldPS).*dfreq,3)-abs(varoutoldPS(:,:,nfreq/2+1,:)).*dfreq),4));

display('........')
display('Difference, global area avg. variance.: ')
CTLglobal=nansum(nansum(tarea.*varoutoldclmo,2),1)./nansum(nansum(tarea,2),1)
LPglobal=nansum(nansum(tarea.*varoutclmo,2),1)./nansum(nansum(tarea,2),1)
ratioout=1-nansum(nansum(tarea.*varoutclmo./varoutoldclmo,2),1)./nansum(nansum(tarea,2),1)
1-LPglobal/CTLglobal

display('<-66:')
tlatgmask=ones(size(tlat));
tlatgmask(tlat>-66)=nan;
CTLglobal=nansum(nansum(tarea.*varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
LPglobal=nansum(nansum(tarea.*varoutclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
ratioglobal=1-nansum(nansum(tarea.*varoutclmo./varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
1-LPglobal/CTLglobal

display('>66:')
tlatgmask=ones(size(tlat));
tlatgmask(tlat<66)=nan;
CTLglobal=nansum(nansum(tarea.*varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
LPglobal=nansum(nansum(tarea.*varoutclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
ratioglobal=1-nansum(nansum(tarea.*varoutclmo./varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
1-LPglobal/CTLglobal


display('abs(lat) < 30:')
tlatgmask=ones(size(tlat));
tlatgmask(abs(tlat)>30)=nan;
CTLglobal=nansum(nansum(tarea.*varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
LPglobal=nansum(nansum(tarea.*varoutclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
ratioglobal=1-nansum(nansum(tarea.*varoutclmo./varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
1-LPglobal/CTLglobal

display('> 30:')
tlatgmask=ones(size(tlat));
tlatgmask(tlat<30)=nan;
CTLglobal=nansum(nansum(tarea.*varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
LPglobal=nansum(nansum(tarea.*varoutclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
ratioglobal=1-nansum(nansum(tarea.*varoutclmo./varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
1-LPglobal/CTLglobal


display('< -30:')
tlatgmask=ones(size(tlat));
tlatgmask(tlat>-30)=nan;
CTLglobal=nansum(nansum(tarea.*varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
LPglobal=nansum(nansum(tarea.*varoutclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
ratioglobal=1-nansum(nansum(tarea.*varoutclmo./varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
1-LPglobal/CTLglobal



clear varoutclmo varoutoldclmo
varoutclmo=squeeze(nanmean((sum(2.*abs(varoutPS(:,:,(nfreq/2+8):end,:)).*dfreq,3)),4));
varoutoldclmo=squeeze(nanmean((sum(2.*abs(varoutoldPS(:,:,(nfreq/2+8):end,:)).*dfreq,3)),4));

display('........')
display('Difference, global area avg. intraseasonal variance: ')
CTLglobal=nansum(nansum(tarea.*varoutoldclmo,2),1)./nansum(nansum(tarea,2),1)
LPglobal=nansum(nansum(tarea.*varoutclmo,2),1)./nansum(nansum(tarea,2),1)
ratioout=1-nansum(nansum(tarea.*varoutclmo./varoutoldclmo,2),1)./nansum(nansum(tarea,2),1)
1-LPglobal/CTLglobal

display('<-66:')
tlatgmask=ones(size(tlat));
tlatgmask(tlat>-66)=nan;
CTLglobal=nansum(nansum(tarea.*varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
LPglobal=nansum(nansum(tarea.*varoutclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
ratioglobal=1-nansum(nansum(tarea.*varoutclmo./varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
1-LPglobal/CTLglobal

display('>66:')
tlatgmask=ones(size(tlat));
tlatgmask(tlat<66)=nan;
CTLglobal=nansum(nansum(tarea.*varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
LPglobal=nansum(nansum(tarea.*varoutclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
ratioglobal=1-nansum(nansum(tarea.*varoutclmo./varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
1-LPglobal/CTLglobal


display('abs(lat) < 30:')
tlatgmask=ones(size(tlat));
tlatgmask(abs(tlat)>30)=nan;
CTLglobal=nansum(nansum(tarea.*varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
LPglobal=nansum(nansum(tarea.*varoutclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
ratioglobal=1-nansum(nansum(tarea.*varoutclmo./varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
1-LPglobal/CTLglobal

display('> 30:')
tlatgmask=ones(size(tlat));
tlatgmask(tlat<30)=nan;
CTLglobal=nansum(nansum(tarea.*varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
LPglobal=nansum(nansum(tarea.*varoutclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
ratioglobal=1-nansum(nansum(tarea.*varoutclmo./varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
1-LPglobal/CTLglobal


display('< -30:')
tlatgmask=ones(size(tlat));
tlatgmask(tlat>-30)=nan;
CTLglobal=nansum(nansum(tarea.*varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
LPglobal=nansum(nansum(tarea.*varoutclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
ratioglobal=1-nansum(nansum(tarea.*varoutclmo./varoutoldclmo.*tlatgmask,2),1)./nansum(nansum(tarea.*tlatgmask,2),1)
1-LPglobal/CTLglobal


