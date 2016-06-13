% calculate the maximum Detection Statistic (DS) for a set of PTA data
% INPUT: 1. x - standard coordinates (0,1) for (delta,alpha,DeltaPhi)
% DeltaPhi (0,2pi) = 2*pi*f*dp*(1-cos\theta)/C, assuming non-evolving source
% 2. inParams - pulsar information and freq-domain data
% OUTPUT: the maximum DS (over freq.)

% Ref. - Zhu et al MNRAS (2016)
function [DS,varargout]=svdDSpso1(x,inParams)
%% transfer parameters from structure inParams
Np = inParams.Np;
N = inParams.N;
deltap = inParams.deltap;
alphap = inParams.alphap;
sigmaT = inParams.sigmaT;
TRwr0 = inParams.TRwr0; % Np*length(w)
xmaxmin = inParams.xmaxmin;

delta=x(1)*(xmaxmin(1,1)-xmaxmin(1,2))+xmaxmin(1,2);  % [-pi/2, pi/2]
alpha=x(2)*(xmaxmin(2,1)-xmaxmin(2,2))+xmaxmin(2,2);  % [0, 2*pi]
DeltaPhi=(x(3:Np+2))'.*(xmaxmin(3:Np+2,1)-xmaxmin(3:Np+2,2))+xmaxmin(3:Np+2,2);

A=zeros(Np,2);
sigmaW=sqrt(N/2).*sigmaT;
for ipsr=1:Np
    [Fp,Fc] = Fpcfun(delta,alpha,deltap(ipsr),alphap(ipsr)); 
    %delta=mod((dp(ipsr)*kpc*(1-cst)/C*(2*pi*omega)),2*pi); % phase difference between ret & rpt
    delta1=DeltaPhi(ipsr);
    A(ipsr,1)=Fp*sin(delta1/2)*(sin(delta1/2)+1i*cos(delta1/2))/sigmaW(ipsr);
    A(ipsr,2)=Fc*sin(delta1/2)*(sin(delta1/2)+1i*cos(delta1/2))/sigmaW(ipsr);
end
%
%%%%%%%%%%%%%
[U D V] = svd(A);
UTRwr=U'*TRwr0;
Udw=(abs(UTRwr)).^2;
svdSNR=Udw(1,:)+Udw(2,:);
DS=-max(svdSNR);

if nargout > 1
    varargout{1}=[delta,alpha,DeltaPhi'];  % in real coord.
end
end