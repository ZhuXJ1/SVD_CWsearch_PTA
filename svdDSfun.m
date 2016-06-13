function [DS,varargout]=svdDSfun(delta,alpha,deltap,alphap,DeltaPhi,w,TRwr0,t,sigmaT)
% calculate the maximum Detection Statistic (DS) for a set of PTA data
% INPUT: (delta,alpha) - (DEC,RA) of the GW source
% (DeltaPhi,alphap) - vector of (DEC,RA) for pulsars
% DeltaPhi (0,2pi) - 2*pi*f*dp*(1-cos\theta)/C, assuming non-evolving source
% data - Np*N, sigmaT - noise rms
% f - GW freq., only used to reduced the size of matrix to aid calculation
% OUTPUT: the maximum DS (over freq.)

% Ref. - Zhu et al MNRAS (2016)

Np=length(deltap);
A=zeros(Np,2);
N=length(t);
sigmaW=sqrt(N/2).*sigmaT;

for ipsr=1:Np
    [Fp,Fc] = Fpcfun(delta,alpha,deltap(ipsr),alphap(ipsr)); 
    delta1=DeltaPhi(ipsr); % phase difference between ret & rpt
    A(ipsr,1)=Fp*sin(delta1/2)*(sin(delta1/2)+1i*cos(delta1/2))/sigmaW(ipsr);
    A(ipsr,2)=Fc*sin(delta1/2)*(sin(delta1/2)+1i*cos(delta1/2))/sigmaW(ipsr);
end
%%%%%%%%%%%%%
[U D V] = svd(A);
UTRwr=U'*TRwr0;
Udw=(abs(UTRwr)).^2;
svdSNR=Udw(1,:)+Udw(2,:);
DS=max(svdSNR);
if nargout > 1
    indi= svdSNR==DS;
    varargout{1}=w(indi);
end
end
