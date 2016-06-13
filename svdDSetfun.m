function [DS,varargout]=svdDSetfun(delta,alpha,deltap,alphap,w,TRwr0,t,sigmaT)
% ONLY search for earth-term signals
% calculate the maximum Detection Statistic (DS) for a set of PTA data
% INPUT: (delta,alpha) - (DEC,RA) of the GW source
% (deltaPhi,alphap) - vector of (DEC,RA) for pulsars
% data - Np*N, sigmaT - noise rms
% f - GW freq., only used to reduced the size of matrix to aid calculation
% OUTPUT: the maximum DS (over freq.)

% Ref. - Zhu et al MNRAS (2015)

Np=length(deltap);
A=zeros(Np,2);
N=length(t);
sigmaW=sqrt(N/2).*sigmaT;

for ipsr=1:Np
    [Fp,Fc] = Fpcfun(delta,alpha,deltap(ipsr),alphap(ipsr)); 
    A(ipsr,:)=[Fp,Fc]/sigmaW(ipsr);
end
%
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