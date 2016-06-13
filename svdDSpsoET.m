% calculate the maximum Detection Statistic (DS) for a set of PTA data, ET search
% INPUT: 1. x - standard coordinates (0,1) for (delta,alpha)
% 2. inParams - pulsar information and freq-domain data
% OUTPUT: the maximum DS (over freq.)

% Ref. - Zhu et al MNRAS (2016)

function [DS,varargout]=svdDSpsoET(x,inParams)
%% transfer parameters from structure inParams
Np = inParams.Np;
N = inParams.N;
deltap = inParams.deltap;
alphap = inParams.alphap;
sigmaT = inParams.sigmaT;
TRwr0 = inParams.TRwr0; % Np*length(w)
xmaxmin = inParams.xmaxmin;

sita0=x(1)*(xmaxmin(1,1)-xmaxmin(1,2))+xmaxmin(1,2);  % [-pi/2, pi/2]
phi0=x(2)*(xmaxmin(2,1)-xmaxmin(2,2))+xmaxmin(2,2);  % [0, 2*pi]

A=zeros(Np,2);
sigmaW=sqrt(N/2).*sigmaT;
for ipsr=1:Np
    [Fp,Fc] = Fpcfun(sita0,phi0,deltap(ipsr),alphap(ipsr)); 
    A(ipsr,:)=[Fp,Fc]/sigmaW(ipsr);
end
%
%%%%%%%%%%%%%
[U D V] = svd(A);
UTRwr=U'*TRwr0;
Udw=(abs(UTRwr)).^2;
svdSNR=Udw(1,:)+Udw(2,:);
DS=-max(svdSNR);

if nargout > 1
    varargout{1}=[sita0,phi0];  % in real coord.
end
end