% stats for plotting

varoutoldann=squeeze(nansum(varoutoldmo.*wtsg,3)./nansum(wtsg,3));
varoutoldannstd=squeeze(nanstd(varoutoldann,0,3));
% clear tout pvalue
% tout=sqrt(5*12)*abs(varoutoldannmean)./varoutoldannstd;
% pvalue=1-tcdf(tout,12*5-1);
% % set pvalue threshold for insignificance
% mask=pvalue>.05;
% varoutoldannmean(mask)=nan;
clear tempvar tempmeanvar
tempvar=varoutoldann;
tempvar=reshape(tempvar,[size(tempvar,1)*size(tempvar,2) size(tempvar,3)]);
tempmeanvar=remove_insignificant_points(tempvar);
varoutoldannmean=reshape(tempmeanvar,[size(varoutoldann,1) size(varoutoldann,2)]);


seasoldmo=squeeze(max(varoutoldmo,[],3)-min(varoutoldmo,[],3));
seasoldmomean1=nanmean(seasoldmo,3);
seasoldmostd=nanstd(seasoldmo,0,3);
% clear tout pvalue
% tout=sqrt(5)*abs(seasoldmomean)./seasoldmostd;
% pvalue=1-tcdf(tout,5-1);
% % set pvalue threshold for insignificance
% mask=pvalue>.05;
% seasoldmomean(mask)=nan;
clear tempvar tempmeanvar
tempvar=seasoldmo;
tempvar=reshape(tempvar,[size(tempvar,1)*size(tempvar,2) size(tempvar,3)]);
tempmeanvar=remove_insignificant_points(tempvar);
seasoldmomean=reshape(tempmeanvar,[size(seasoldmo,1) size(seasoldmo,2)]);



diffmo = varoutoldmo-varoutmo;
diffmomean1=nansum(wtsg60.*reshape(diffmo,[size(diffmo,1) size(diffmo,2) size(diffmo,3).*size(diffmo,4)]),3)./...
    nansum(wtsg60,3);
diffmostd=nanstd(reshape(diffmo,[size(diffmo,1) size(diffmo,2) size(diffmo,3).*size(diffmo,4)]),0,3);
% clear tout pvalue
% tout=sqrt(5*12)*abs(diffmomean)./diffmostd;
% pvalue=1-tcdf(tout,5*12-1);
% % set pvalue threshold for insignificance
% mask=pvalue>.05;
% diffmomean(mask)=nan;
clear tempvar tempmeanvar
tempvar=diffmo;
tempvar=reshape(tempvar,[size(tempvar,1)*size(tempvar,2) size(tempvar,3)*size(tempvar,4)]);
tempmeanvar=remove_insignificant_points(tempvar);
diffmomean=reshape(tempmeanvar,[size(diffmo,1) size(diffmo,2)]);

diffseasmo=squeeze(max(varoutoldmo,[],3)-min(varoutoldmo,[],3)-max(varoutmo,[],3)+min(varoutmo,[],3));
diffseasmomean1=nanmean(diffseasmo,3);
diffseasmostd=nanstd(diffseasmo,0,3);
% clear tout pvalue mask
% tout=sqrt(5)*abs(diffseasmomean)./diffseasmostd;
% pvalue=1-tcdf(tout,5-1);
% % set pvalue threshold for insignificance
% mask=pvalue>.05;
% diffseasmomean(mask)=nan;
clear tempvar tempmeanvar
tempvar=diffseasmo;
tempvar=reshape(tempvar,[size(tempvar,1)*size(tempvar,2) size(tempvar,3)*size(tempvar,4)]);
tempmeanvar=remove_insignificant_points(tempvar);
diffseasmomean=reshape(tempmeanvar,[size(diffseasmo,1) size(diffseasmo,2)]);
