%% source parameters, a circuar SMBHB at the Virgo cluster
raj=3.2594; % Virgo cluster
dec=0.2219;
cosi=0.88;
psi=0.5;
phi0=2.89;
evlInd=0; % non-evolving, 1 for evolving
% h0=1e-13/7.44; % scaled for desired S/N, for 12p
h0=1e-13/7.97; % scaled for desired S/N, for 12p
f=1e-8;
%% observations
fid = fopen('IPTA30p.txt');
psrDat=textscan(fid,'%s %f64 %f64 %f64 %f64');
fclose(fid);
Jname=psrDat{1};
alphap=psrDat{2};
deltap=psrDat{3};
dp=psrDat{4};
sigmaT=psrDat{5};
Np=length(alphap);
% observing deltapan and time vector
datlen=15*365.25;
t=86400.*(0:14:datlen);% 15 yr, two-week cadence
N=length(t);
% some constants
C = 299792458;
kpc=1e3 * 3.08568025 * 10^(16) ;
solar_mass = 1.989e30; 
G = 6.67384e-11;
Mpc=1e6 * 3.08568025 * 10^(16) ;
Dl=16.5*Mpc; % work out chirp mass for given h0
Mc=(h0/(2*((G*solar_mass/C^2).^(5/3))*(1/Dl)*((pi*f).^(2/3))*(C^(-2/3))))^(3/5);
%% calculate the statistic, assuming all source parameters are known
oS=4; % frequency over-Sampling factor
mTRt=zeros(Np,N);
A=zeros(Np,2); % Earth-term search
A1=zeros(Np,2);
nt0=randn(Np,N);
nt=zeros(Np,N);
sigmaW=sqrt(N/2)*sigmaT;
SNRsq=zeros(1,Np);
DeltaPhi=zeros(1,Np);
for ipsr=1:Np
    [Fp Fc cst] = Fpcfun(dec,raj,deltap(ipsr),alphap(ipsr)); 
    [ret rpt TRt]=TimResfun(dec,raj,deltap(ipsr),alphap(ipsr),dp(ipsr),h0,f,cosi,psi,phi0,t,evlInd,Mc);
    mTRt(ipsr,:)=TRt;
    A(ipsr,:)=[Fp,Fc]/sigmaW(ipsr);
    nt(ipsr,:) = sigmaT(ipsr).*nt0(ipsr,:);
    
    SNRsq1=(1/sigmaT(ipsr)).*TRt;
    SNRsq(ipsr)=sum(SNRsq1.^2);
    delta1=mod((dp(ipsr)*kpc*(1-cst)/C*(2*pi*f)),2*pi); %phase difference between ret & rpt
    DeltaPhi(ipsr)=delta1;
    A1(ipsr,1)=Fp*sin(delta1/2)*(sin(delta1/2)+1i*cos(delta1/2))/sigmaW(ipsr);
    A1(ipsr,2)=Fc*sin(delta1/2)*(sin(delta1/2)+1i*cos(delta1/2))/sigmaW(ipsr);
end
%
indSNR=sqrt(SNRsq);
netSNR=sqrt(sum(SNRsq));
dat=mTRt;%+nt
[TRwr0 w]=fftnfun(dat,t,oS,sigmaT);
% % reduce the matrix size, for a better-look plot
% indi=find(w>f*0.6 & w<f*1.4);
% w=w(indi);
% TRwr0=TRwr0(:,indi);
%%%%%%%%%%%%%
[U D V] = svd(A1); % replace A1 with A for Earth-term search
UTRwr=U'*TRwr0;
Udw=(abs(UTRwr)).^2;
svdSNR=Udw(1,:)+Udw(2,:);
realDS=max(svdSNR); % detection statistic
figure
semilogx(w,svdSNR,'b')
%% perform a search using PSO, consider both ET and full-signal search
xmaxmin=zeros(Np+2,2);  % x_max, x_min for each parameter x
xmaxmin(2,1)=2*pi;  % alpha
xmaxmin(2,2)=0.0;
xmaxmin(1,1)=pi/2;  % delta 
xmaxmin(1,2)=-pi/2;
xmaxmin(3:Np+2,1)=2*pi;  % DeltaPhi 
xmaxmin(3:Np+2,2)=0;
% frequency has been internally maximized in the fitness funtion,
% so the dimension is 2 and Np+2 for E-T search and full search respectively
nDim = length(xmaxmin(:,1));
nDim1=2;
[SNR4,freq4]=svdDSetfun(dec,raj,deltap,alphap,w,TRwr0,t,sigmaT);
[SNR2,freq2]=svdDSfun(dec,raj,deltap,alphap,DeltaPhi,w,TRwr0,t,sigmaT);
inParams = struct('Np',Np,'N',N,'deltap',deltap,'alphap',alphap,...
                    'sigmaT',sigmaT,'TRwr0',TRwr0,'xmaxmin',xmaxmin);%,'dp',dp
fHandle = @(x) svd2pso(x,inParams);
fHandle1 = @(x) svd2psoET(x,inParams);
tic
psoResults=ptapso(fHandle,nDim);
outParam=psoResults.bestLocation;
maxSNR = psoResults.bestSNR;
psoResults1=ptapso1(fHandle1,nDim1);
outParam1=psoResults1.bestLocation;
maxSNR1 = psoResults1.bestSNR;
pos1=outParam(1)*(xmaxmin(1,1)-xmaxmin(1,2))+xmaxmin(1,2);
pos2=outParam(2)*(xmaxmin(2,1)-xmaxmin(2,2))+xmaxmin(2,2);
pha=outParam(3:Np+2).*(xmaxmin(3:Np+2,1)-xmaxmin(3:Np+2,2))'+(xmaxmin(3:Np+2,2))'; % pulsar distance
SNR1=-maxSNR;
pos3=outParam1(1)*(xmaxmin(1,1)-xmaxmin(1,2))+xmaxmin(1,2);
pos4=outParam1(2)*(xmaxmin(2,1)-xmaxmin(2,2))+xmaxmin(2,2);
SNR3=-maxSNR1;
% just double check, such that
% SNR5==SNR1; SNR6==SNR3
[SNR5,freq1]=svdDSfun(pos1,pos2,deltap,alphap,pha,w,TRwr0,t,sigmaT);
[SNR6,freq3]=svdDSetfun(pos3,pos4,deltap,alphap,w,TRwr0,t,sigmaT);
toc