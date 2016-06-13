function [Fp,Fc,varargout]=Fpcfun(delta,alpha,deltap,alphap)
% Antenna pattern function for pulsar timing, F_plus, F_cross
% defined in Eqs(2-3) in Zhu et al.MNRAS(2015,2016), see also Lee et al.MNRAS(2011)

% delta alpha - GW source Dec and RA (in radians)
% deltap alphap - pulsar Dec and RA (in radians)
% delta & deltap [-pi/2,pi/2], alpha & alphap & psi [0,2*pi)

cst=cos(delta)*cos(deltap)*cos(alpha-alphap)+sin(delta)*sin(deltap);
% cosin of delta (angle between GW source and Pulsar)
if (cst==1)
    Fp=0;
    Fc=0;
else
Fp0=(1+(sin(delta)).^2)*((cos(deltap)).^2)*cos(2*(alpha-alphap))-sin(2*delta)*sin(2*deltap)*cos(alpha-alphap)+(2-3*((cos(deltap)).^2))*((cos(delta)).^2);
Fc0=2*cos(delta)*sin(2*deltap)*sin(alpha-alphap)-2*sin(delta)*((cos(deltap)).^2)*sin(2*(alpha-alphap));
Fp=Fp0./(4*(1-cst));
Fc=Fc0./(4*(1-cst));
end

if nargout > 1
    varargout{1}=cst;
end
return