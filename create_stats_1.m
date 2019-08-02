% stats for plotting

varoutoldann=squeeze(nansum(varoutoldmo.*wtsg,3)./nansum(wtsg,3));
varoutoldannmean=squeeze(nanmean(varoutoldann,3));
varoutoldannstd=squeeze(nanstd(varoutoldann,0,3));


seasoldmo=squeeze(max(varoutoldmo,[],3)-min(varoutoldmo,[],3));
seasoldmomean=nanmean(seasoldmo,3);
seasoldmostd=nanstd(seasoldmo,0,3);




diffmo = varoutoldmo-varoutmo;
diffmomean1=nansum(wtsg60.*reshape(diffmo,[size(diffmo,1) size(diffmo,2) size(diffmo,3).*size(diffmo,4)]),3)./...
    nansum(wtsg60,3);
diffmostd=nanstd(reshape(diffmo,[size(diffmo,1) size(diffmo,2) size(diffmo,3).*size(diffmo,4)]),0,3);
clear tempvar tempmeanvar
tempvar=diffmo;
tempvar=reshape(tempvar,[size(tempvar,1)*size(tempvar,2) size(tempvar,3)*size(tempvar,4)]);
tempmeanvar=remove_insignificant_points(tempvar,1);
diffmomean=reshape(tempmeanvar,[size(diffmo,1) size(diffmo,2)]);

diffseasmo=squeeze(max(varoutoldmo,[],3)-min(varoutoldmo,[],3)-max(varoutmo,[],3)+min(varoutmo,[],3));
diffseasmomean1=nanmean(diffseasmo,3);
diffseasmostd=nanstd(diffseasmo,0,3);
clear tempvar tempmeanvar
tempvar=diffseasmo;
tempvar=reshape(tempvar,[size(tempvar,1)*size(tempvar,2) size(tempvar,3)*size(tempvar,4)]);
tempmeanvar=remove_insignificant_points(tempvar,1);
diffseasmomean=reshape(tempmeanvar,[size(diffseasmo,1) size(diffseasmo,2)]);
