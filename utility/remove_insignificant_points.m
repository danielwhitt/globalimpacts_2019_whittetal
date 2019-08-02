function [Amean] = remove_insignificant_points(Ain,flag)
%%
% flag=0 = t test, 1 = bootstrap
Amean=squeeze(nanmean(Ain,2));
Ain(isnan(Ain))=0;
if flag == 0
    Astd=squeeze(nanstd(Ain,0,2));
    clear tout pvalue
    tout=sqrt(size(Ain,2))*abs(Amean)./Astd;
    pvalue=1-tcdf(tout,size(Ain,2)-1);
    % set pvalue threshold for insignificance
    mask=pvalue>.05;
    Amean(mask)=nan;
elseif flag == 1
    options1=statset('UseParallel',logical(1));
    ci = bootci(1000,{@nanmean,Ain'},'Options',options1);
    clear mask
    mask = max(abs(ci),[],1)==0;
    ci(:,mask)=nan;
    Amean(mask,:)=nan;
    clear mask
    mask = sign(min(ci,[],1))~=sign(max(ci,[],1));
    Amean(mask,:)=nan;
    ci(:,mask)=nan;
    ci = ci';
end
end

