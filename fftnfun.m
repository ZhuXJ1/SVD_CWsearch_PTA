function [TRwr0 w]=fftnfun(dat,t,oS,sigmaT)
% normalized fft for evenly-sampled PTA data (assuming equal error bars and only white noise)

% dat - data Np*Nobs, t - time Nobs*1 (assumed to be the same for pulsars)
% oS - frequency oversampling factor, sigmaT - noise rms Np*1
% TRwr0 - normalized fft of data, w - frequency vector
Np=length(sigmaT);
N=length(t);
Nmid=floor(N/2);
T=t(N)-t(1);

Df=1/(oS*T);
w=Df.*(1:(Nmid*oS));
TRwr=zeros(Np,N*oS);
DmTRt=zeros(Np,oS*N);
DmTRt(:,(oS*Nmid+1):(oS*Nmid+N))=dat;
for i=1:Np
    TRwr(i,:)=(fft(DmTRt(i,:)))/(sqrt(N/2)*sigmaT(i));
end
TRwr0=TRwr(:,2:(Nmid*oS+1));
return