
tareagnow=tareag;
tareagnow(isnan(varoutoldclmo))=nan;
display('Difference, global area mean: ')
CTLglobal=nansum(nansum(nansum(wtsg.*tareag.*varoutoldclmo,3),2),1)./nansum(nansum(nansum(wtsg.*tareagnow,3),2),1)
tareagnow=tareag;
tareagnow(isnan(varoutclmo))=nan;
LPglobal=nansum(nansum(nansum(wtsg.*tareag.*varoutclmo,3),2),1)./nansum(nansum(nansum(wtsg.*tareagnow,3),2),1)
1-LPglobal./CTLglobal

display('<-66:')
tlatgmask=ones(size(tlatg));
tlatgmask(tlatg>-66)=nan;
tareagnow=tareag;
tareagnow(isnan(varoutoldclmo))=nan;
CTLglobal=nansum(nansum(nansum(wtsg.*tareag.*varoutoldclmo.*tlatgmask,3),2),1)./nansum(nansum(nansum(wtsg.*tareagnow.*tlatgmask,3),2),1)
tareagnow=tareag;
tareagnow(isnan(varoutclmo))=nan;
LPglobal=nansum(nansum(nansum(wtsg.*tareag.*varoutclmo.*tlatgmask,3),2),1)./nansum(nansum(nansum(wtsg.*tareagnow.*tlatgmask,3),2),1)
1-LPglobal./CTLglobal

display('>66:')
tlatgmask=ones(size(tlatg));
tlatgmask(tlatg<66)=nan;
tareagnow=tareag;
tareagnow(isnan(varoutoldclmo))=nan;
CTLglobal=nansum(nansum(nansum(wtsg.*tareag.*varoutoldclmo.*tlatgmask,3),2),1)./nansum(nansum(nansum(wtsg.*tareagnow.*tlatgmask,3),2),1)
tareagnow=tareag;
tareagnow(isnan(varoutclmo))=nan;
LPglobal=nansum(nansum(nansum(wtsg.*tareag.*varoutclmo.*tlatgmask,3),2),1)./nansum(nansum(nansum(wtsg.*tareagnow.*tlatgmask,3),2),1)
1-LPglobal./CTLglobal

display('abs(lat) < 30:')
tlatgmask=ones(size(tlatg));
tlatgmask(abs(tlatg)>30)=nan;
tareagnow=tareag;
tareagnow(isnan(varoutoldclmo))=nan;
CTLglobal=nansum(nansum(nansum(wtsg.*tareag.*varoutoldclmo.*tlatgmask,3),2),1)./nansum(nansum(nansum(wtsg.*tareagnow.*tlatgmask,3),2),1)
tareagnow=tareag;
tareagnow(isnan(varoutclmo))=nan;
LPglobal=nansum(nansum(nansum(wtsg.*tareag.*varoutclmo.*tlatgmask,3),2),1)./nansum(nansum(nansum(wtsg.*tareagnow.*tlatgmask,3),2),1)
1-LPglobal./CTLglobal

display('> 30:')
tlatgmask=ones(size(tlatg));
tlatgmask(tlatg<30)=nan;
tareagnow=tareag;
tareagnow(isnan(varoutoldclmo))=nan;
CTLglobal=nansum(nansum(nansum(wtsg.*tareag.*varoutoldclmo.*tlatgmask,3),2),1)./nansum(nansum(nansum(wtsg.*tareagnow.*tlatgmask,3),2),1)
tareagnow=tareag;
tareagnow(isnan(varoutclmo))=nan;
LPglobal=nansum(nansum(nansum(wtsg.*tareag.*varoutclmo.*tlatgmask,3),2),1)./nansum(nansum(nansum(wtsg.*tareagnow.*tlatgmask,3),2),1)
1-LPglobal./CTLglobal


display('< -30:')
tlatgmask=ones(size(tlatg));
tlatgmask(tlatg>-30)=nan;
tareagnow=tareag;
tareagnow(isnan(varoutoldclmo))=nan;
CTLglobal=nansum(nansum(nansum(wtsg.*tareag.*varoutoldclmo.*tlatgmask,3),2),1)./nansum(nansum(nansum(wtsg.*tareagnow.*tlatgmask,3),2),1)
tareagnow=tareag;
tareagnow(isnan(varoutclmo))=nan;
LPglobal=nansum(nansum(nansum(wtsg.*tareag.*varoutclmo.*tlatgmask,3),2),1)./nansum(nansum(nansum(wtsg.*tareagnow.*tlatgmask,3),2),1)
1-LPglobal./CTLglobal