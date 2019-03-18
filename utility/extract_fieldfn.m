function var = extract_fieldfn(fn,varname,stidxs,lenidxs,strideidxs)
if ~exist('strideidxs','var') || isempty(strideidxs)
  strideidxs=ones(1,length(lenidxs));
end
ncid = netcdf.open(fn,'NC_NOWRITE');
varID = netcdf.inqVarID(ncid,varname);
var = netcdf.getVar(ncid,varID,stidxs,lenidxs./strideidxs,strideidxs,'double');
netcdf.close(ncid);
end
