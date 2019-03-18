% compute spectral components
function [y,freq] = fftfuncore(x,Fs)
x=permute(x,[3 1 2]);
N = size(x,1);
%Qf = .94;
Qf=1;
%wind = tukeywin(N,.2);
wind = ones(N,1);
wind = repmat(wind, [1 size(x,2) size(x,3)]); 
xdft = 1./Fs.*fft(wind.*x,[],1); 
xdft = fftshift(xdft(:,:,:),1);
y = (1/(Qf)).*xdft;
dF = Fs/size(x,1);
freq = (-Fs/2):dF:(Fs/2-dF);
y=permute(y,[2 3 1]);
end

