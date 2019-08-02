% stats for plotting
% 3,2,1
fractotall=sqrt(sum(abs(varoutoldPS).*dfreq,3)-abs(varoutoldPS(:,:,nfreq/2+1,:)).*dfreq);
fractotallmean=squeeze(nanmean(fractotall,4));
fractotallstd=squeeze(nanstd(fractotall,0,4));
clear tout pvalue
tout=sqrt(5)*abs(fractotallmean)./fractotallstd;
pvalue=1-tcdf(tout,5-1);
% set pvalue threshold for insignificance
mask=pvalue>.05;
fractotallmean(mask)=nan;


% fractotallLP=sqrt(sum(abs(varoutPS).*dfreq,3)-abs(varoutPS(:,:,nfreq/2+1,:)).*dfreq);
% fractotallLPmean=squeeze(nanmean(fractotallLP,4));
% fractotallLPstd=squeeze(nanstd(fractotallLP,0,4));
% clear tout pvalue
% tout=sqrt(5)*abs(fractotallLPmean)./fractotallLPstd;
% pvalue=1-tcdf(tout,5-1);
% % set pvalue threshold for insignificance
% mask=pvalue>.05;
% fractotallLPmean(mask)=nan;

% 3,2,2
difftotvarmo= fractotall -sqrt((sum(abs(varoutPS).*dfreq,3)-abs(varoutPS(:,:,nfreq/2+1,:)).*dfreq));
difftotvarmomean=squeeze(nanmean(difftotvarmo,4));
difftotvarmostd=squeeze(nanstd(difftotvarmo,0,4));
clear tout pvalue
tout=sqrt(5)*abs(difftotvarmomean)./difftotvarmostd;
pvalue=1-tcdf(tout,5-1);
% set pvalue threshold for insignificance
mask=pvalue>.05;
difftotvarmomean(mask)=nan;

% 3,2,3
fracintravar=sqrt(sum(2.*abs(varoutoldPS(:,:,(nfreq/2+8):end,:)).*dfreq,3));
fracintravarmean=squeeze(nanmean(fracintravar,4));
fracintravarstd=squeeze(nanstd(fracintravar,0,4));
clear tout pvalue
tout=sqrt(5)*abs(fracintravarmean)./fracintravarstd;
pvalue=1-tcdf(tout,5-1);
% set pvalue threshold for insignificance
mask=pvalue>.05;
fracintravarmean(mask)=nan;

% fracintravarLP=sqrt(sum(2.*abs(varoutPS(:,:,(nfreq/2+8):end,:)).*dfreq,3));
% fracintravarLPmean=squeeze(nanmean(fracintravarLP,4));
% fracintravarLPstd=squeeze(nanstd(fracintravarLP,0,4));
% clear tout pvalue
% tout=sqrt(5)*abs(fracintravarLPmean)./fracintravarLPstd;
% pvalue=1-tcdf(tout,5-1);
% % set pvalue threshold for insignificance
% mask=pvalue>.05;
% fracintravarLPmean(mask)=nan;

% 3,2,4
diffintravarmo=fracintravar-sqrt(sum(2.*abs(varoutPS(:,:,(nfreq/2+8):end,:)).*dfreq,3));
diffintravarmomean=squeeze(nanmean(diffintravarmo,4));
diffintravarmostd=squeeze(nanstd(diffintravarmo,0,4));
clear tout pvalue
tout=sqrt(5)*abs(diffintravarmomean)./diffintravarmostd;
pvalue=1-tcdf(tout,5-1);
% set pvalue threshold for insignificance
mask=pvalue>.05;
diffintravarmomean(mask)=nan;


% 3,2,5
varoutPSmean=squeeze(nanmean(varoutPS,4));
varoutoldPSmean=squeeze(nanmean(varoutoldPS,4));

tareag(isnan(varoutoldPSmean))=nan;

tareagnow=tareag;
tareagnow(isnan(varoutPSmean))=nan;

areaweightedPSmean=squeeze(nansum(reshape(2.*abs(varoutPSmean).*tareagnow,...
    [size(varoutPSmean,1)*size(varoutPSmean,2) size(varoutPSmean,3)]),1)./...
    nansum(reshape(tareagnow,[size(varoutPSmean,1)*size(varoutPSmean,2) size(varoutPSmean,3)]),1));

tareagnow=tareag;
tareagnow(isnan(varoutoldPSmean))=nan;
areaweightedoldPSmean=squeeze(nansum(reshape(2.*abs(varoutoldPSmean).*tareagnow,...
    [size(varoutPSmean,1)*size(varoutPSmean,2) size(varoutPSmean,3)]),1)./...
    nansum(reshape(tareagnow,[size(varoutPSmean,1)*size(varoutPSmean,2) size(varoutPSmean,3)]),1));

varoutoldPSiqr= quantile(squeeze(reshape(...
    (2.*abs(varoutoldPSmean)),...
    [size(varoutPSmean,1)*size(varoutPSmean,2) size(varoutPSmean,3)])),[.05 .95],1);

varoutPSiqr= quantile(squeeze(reshape(...
    (2.*abs(varoutPSmean)),...
    [size(varoutPSmean,1)*size(varoutPSmean,2) size(varoutPSmean,3)])),[.05 .95], 1);

% 3,2,6
fracvar=cumsum(2.*abs(varoutPS(:,:,(nfreq/2+2):nfreq,:)).*dfreq,3);
fracvar=fracvar./repmat(fracvar(:,:,end,:),[1 1 (nfreq/2-1) 1]);
fracvarold=cumsum(2.*abs(varoutoldPS(:,:,(nfreq/2+2):nfreq,:)).*dfreq,3);
fracvarold=fracvarold./repmat(fracvarold(:,:,end,:),[1 1 (nfreq/2-1) 1]);
fracvarmean=squeeze(nanmean(fracvar,4));
fracvaroldmean=squeeze(nanmean(fracvarold,4));

%
display('some statistics:')
display('global area avg. ratio of variance in the intraseasonal band:')
tareagnow=tareag;
tareagnow(isnan(varoutoldPSmean))=nan;
squeeze(nansum(reshape(...
    (nansum(abs(varoutoldPSmean(:,:,(nfreq/2+8):nfreq)),3)./...
    nansum(abs(varoutoldPSmean(:,:,(nfreq/2+2):nfreq)),3))...
    .*tareagnow(:,:,1),...
    [size(varoutPSmean,1)*size(varoutPSmean,2) 1]),1)./...
    nansum(reshape(tareagnow(:,:,1),[size(varoutPSmean,1)*size(varoutPSmean,2) 1]),1))
tempvar=(nansum(abs(varoutoldPSmean(:,:,(nfreq/2+8):nfreq)),3)./...
    nansum(abs(varoutoldPSmean(:,:,(nfreq/2+2):nfreq)),3));


tareagnow=tareag;
tareagnow(isnan(varoutoldPSmean))=nan;
tareagnow=squeeze(tareagnow(:,:,1));
tareagnow2=tareagnow;
display('frac of ocean area where intraseasonal band > half the variance:')
tareagnow2(tempvar<0.5)=nan;
squeeze(nansum(reshape(...
    tareagnow2,...
    [size(varoutPSmean,1)*size(varoutPSmean,2) 1]),1)./...
    nansum(reshape(tareagnow,[size(varoutPSmean,1)*size(varoutPSmean,2) 1]),1))

