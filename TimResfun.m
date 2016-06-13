%----------return GW-induced pulsar timing redisuals from circular binaries-------
% ret(earth term) rpt(pulsar term) TRt(full term) 

function [ret rpt TRt]=TimResfun(delta,alpha,deltap,alphap,dp,h0,f,cosi,psi,phi0,t,evlInd,Mc)
% pulsar timing residuals induced by a circular SMBH binary at (RA,DEC) - (alpha,delta)
% evlInd - indicator for evolution (0 for non-evolving, 1 for evolving) over pulsar-Earth baseline
% (alphap,deltap) - (RA,DEC) for pulsar, dp - pulsar distance in kpc
% h0 - GW strain amplitude, f - GW frequency, cosi - cosin of inclination angle
% psi - polarization angle, phi0 - phase constant
% t - time vector in seconds, Mc - binary chirp mass in solar mass

% Refs: Zhu et al MNRAS (2015,2016), Contact - Xingjiang Zhu (zhuxingjiang@gmail.com)
% Babak&Sesana PRD 2012, Ellis et al ApJ2012,  Wang et al ApJ (2014,2015)

% some constants
C = 299792458;
kpc=1e3 * 3.08568025 * 10^(16) ;
solar_mass = 1.989*10^(30); 
G = 6.67384.*10.^(-11);
M=Mc*solar_mass;

[Fp Fc cst]=Fpcfun(delta,alpha,deltap,alphap);

Ap=(h0/(2*pi*f)).*((1+(cosi^2)).*cos(2*psi).*sin(2*pi*f.*t+phi0)+2*cosi*sin(2*psi).*cos(2*pi*f.*t+phi0));
Ac=(h0/(2*pi*f)).*((1+(cosi^2)).*sin(2*psi).*sin(2*pi*f.*t+phi0)-2*cosi*cos(2*psi).*cos(2*pi*f.*t+phi0));

ret=Fp.*Ap+Fc.*Ac;      % earth term signal
tp=t-dp*kpc*(1-cst)/C;  % time vector for pulsar term

if (evlInd==0)
    Aalphap=(h0/(2*pi*f)).*((1+(cosi^2)).*cos(2*psi).*sin(2*pi*f.*tp+phi0)+2*cosi*sin(2*psi).*cos(2*pi*f.*tp+phi0));
    Acp=(h0/(2*pi*f)).*((1+(cosi^2)).*sin(2*psi).*sin(2*pi*f.*tp+phi0)-2*cosi*cos(2*psi).*cos(2*pi*f.*tp+phi0));
    rpt=Fp.*Aalphap+Fc.*Acp;
else
%     df=1.5e-8*((f/5e-8).^(11/3))*dp*(1-cst);
%     fp=f-df; % frequency of pulsar term signal, aalphaproximation hold only
%     when df is small
    fp=((dp*kpc*(1-cst)/C)*(256/5)*(pi^(8/3))*(C^(-5))*((G*M)^(5/3))+f^(-8/3)).^(-3/8);
    hp=h0*((fp/f)^(2/3));
    Aalphap=(hp/(2*pi*fp)).*((1+(cosi^2)).*cos(2*psi).*sin(2*pi*fp.*tp+phi0)+2*cosi*sin(2*psi).*cos(2*pi*fp.*tp+phi0));
    Acp=(hp/(2*pi*fp)).*((1+(cosi^2)).*sin(2*psi).*sin(2*pi*fp.*tp+phi0)-2*cosi*cos(2*psi).*cos(2*pi*fp.*tp+phi0));
    rpt=Fp.*Aalphap+Fc.*Acp;
end
TRt=ret-rpt;    % full term signal
return

