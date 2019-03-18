%% create monthly means

wts=[31 28 31 30 31 30 31 31 30 31 30 31]';
wtsg=repmat(wts,[1 size(varoutold,1) size(varoutold,2)]);
wtsg=permute(wtsg,[2 3 1]);
wtsg60=repmat(wtsg,[1 1 5]);

for yrct=1:5
    for moct = 1:12
        varoutmo(:,:,moct,yrct)=nansum(varout(:,:,...
            (365*(yrct-1)+sum(wts(1:(moct)-1))+1):(365*(yrct-1)+sum(wts(1:(moct))))),3)./nansum(wts(moct));
        varoutoldmo(:,:,moct,yrct)=nansum(varoutold(:,:,...
            (365*(yrct-1)+sum(wts(1:(moct)-1))+1):(365*(yrct-1)+sum(wts(1:(moct))))),3)./nansum(wts(moct));
    end
end
