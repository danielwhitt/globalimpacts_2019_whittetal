% global area and regional average statistics
tarea(isnan(squeeze(varoutold(:,:,1))))=nan;
tareag=repmat(tarea,[1 1 12]);
clear varoutclmo varoutoldclmo
varoutclmo=squeeze(nanmean(varoutmo,4));
varoutoldclmo=squeeze(nanmean(varoutoldmo,4));

display('Difference, global area avg. annual mean: ')

CTLglobal=nansum(nansum(nansum(wtsg.*tareag.*varoutoldclmo,3),2),1)./nansum(nansum(nansum(wtsg.*tareag,3),2),1)
LPglobal=nansum(nansum(nansum(wtsg.*tareag.*varoutclmo,3),2),1)./nansum(nansum(nansum(wtsg.*tareag,3),2),1)
ratioglobal=1-nansum(nansum(nansum(wtsg.*tareag.*varoutclmo./varoutoldclmo,3),2),1)./nansum(nansum(nansum(wtsg.*tareag,3),2),1)
1-LPglobal/CTLglobal

display('<-66:')
tlatgmask=ones(size(tlatg));
tlatgmask(tlatg>-66)=nan;
CTLglobal=nansum(nansum(nansum(wtsg.*tareag.*varoutoldclmo.*tlatgmask,3),2),1)./nansum(nansum(nansum(wtsg.*tareag.*tlatgmask,3),2),1)
LPglobal=nansum(nansum(nansum(wtsg.*tareag.*varoutclmo.*tlatgmask,3),2),1)./nansum(nansum(nansum(wtsg.*tareag.*tlatgmask,3),2),1)
ratioglobal=1-nansum(nansum(nansum(wtsg.*tareag.*varoutclmo./varoutoldclmo.*tlatgmask,3),2),1)./nansum(nansum(nansum(wtsg.*tareag.*tlatgmask,3),2),1)
1-LPglobal/CTLglobal

display('>66:')
tlatgmask=ones(size(tlatg));
tlatgmask(tlatg<66)=nan;
CTLglobal=nansum(nansum(nansum(wtsg.*tareag.*varoutoldclmo.*tlatgmask,3),2),1)./nansum(nansum(nansum(wtsg.*tareag.*tlatgmask,3),2),1)
LPglobal=nansum(nansum(nansum(wtsg.*tareag.*varoutclmo.*tlatgmask,3),2),1)./nansum(nansum(nansum(wtsg.*tareag.*tlatgmask,3),2),1)
ratioglobal=1-nansum(nansum(nansum(wtsg.*tareag.*varoutclmo./varoutoldclmo.*tlatgmask,3),2),1)./nansum(nansum(nansum(wtsg.*tareag.*tlatgmask,3),2),1)
1-LPglobal/CTLglobal


display('abs(lat) < 30:')
tlatgmask=ones(size(tlatg));
tlatgmask(abs(tlatg)>30)=nan;
CTLglobal=nansum(nansum(nansum(wtsg.*tareag.*varoutoldclmo.*tlatgmask,3),2),1)./nansum(nansum(nansum(wtsg.*tareag.*tlatgmask,3),2),1)
LPglobal=nansum(nansum(nansum(wtsg.*tareag.*varoutclmo.*tlatgmask,3),2),1)./nansum(nansum(nansum(wtsg.*tareag.*tlatgmask,3),2),1)
ratioglobal=1-nansum(nansum(nansum(wtsg.*tareag.*varoutclmo./varoutoldclmo.*tlatgmask,3),2),1)./nansum(nansum(nansum(wtsg.*tareag.*tlatgmask,3),2),1)
1-LPglobal/CTLglobal

display('> 30:')
tlatgmask=ones(size(tlatg));
tlatgmask(tlatg<30)=nan;
CTLglobal=nansum(nansum(nansum(wtsg.*tareag.*varoutoldclmo.*tlatgmask,3),2),1)./nansum(nansum(nansum(wtsg.*tareag.*tlatgmask,3),2),1)
LPglobal=nansum(nansum(nansum(wtsg.*tareag.*varoutclmo.*tlatgmask,3),2),1)./nansum(nansum(nansum(wtsg.*tareag.*tlatgmask,3),2),1)
ratioglobal=1-nansum(nansum(nansum(wtsg.*tareag.*varoutclmo./varoutoldclmo.*tlatgmask,3),2),1)./nansum(nansum(nansum(wtsg.*tareag.*tlatgmask,3),2),1)
1-LPglobal/CTLglobal


display('< -30:')
tlatgmask=ones(size(tlatg));
tlatgmask(tlatg>-30)=nan;
CTLglobal=nansum(nansum(nansum(wtsg.*tareag.*varoutoldclmo.*tlatgmask,3),2),1)./nansum(nansum(nansum(wtsg.*tareag.*tlatgmask,3),2),1)
LPglobal=nansum(nansum(nansum(wtsg.*tareag.*varoutclmo.*tlatgmask,3),2),1)./nansum(nansum(nansum(wtsg.*tareag.*tlatgmask,3),2),1)
ratioglobal=1-nansum(nansum(nansum(wtsg.*tareag.*varoutclmo./varoutoldclmo.*tlatgmask,3),2),1)./nansum(nansum(nansum(wtsg.*tareag.*tlatgmask,3),2),1)
1-LPglobal/CTLglobal




clear varoutclmo varoutoldclmo
varoutclmo=squeeze(nanmean(squeeze(max(varoutmo,[],3)-min(varoutmo,[],3)),3));
varoutoldclmo=squeeze(nanmean(squeeze(max(varoutoldmo,[],3)-min(varoutoldmo,[],3)),3));



display('Difference, global area avg. seas cycl ampl.: ')
CTLglobal=nansum(nansum(tarea.*varoutoldclmo,2),1)./nansum(nansum(tarea,2),1)
LPglobal=nansum(nansum(tarea.*varoutclmo,2),1)./nansum(nansum(tarea,2),1)
ratioglobal=1-nansum(nansum(tarea.*varoutclmo./varoutoldclmo,2),1)./nansum(nansum(tarea,2),1)
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
