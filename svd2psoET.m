% Calculate the statistic of the SVD method using particle swarm optimazition (PSO)
% PSO code was developed by Prof. Soumya Mohanty (soumya.mohanty@utrgv.edu)
% used in Zhu et al. MNRAS (2016), also in Wang et al. (2014,2015)
% connector between SVD and PSO, for ET search

function [fitVal,varargout]=svd2psoET(xVec,inParams)
%function fitVal=LLR_PSO(xVec,inParams)

%rows: points
%columns: coordinates of a point
[nrows,ncols]=size(xVec);

%storage for fitness values
fitVal = zeros(nrows,1);

validPts = chkstdsrchrng(xVec);
% %Set fitness for invalid points to infty
fitVal(~validPts)=inf;
% if isempty(inParams)
%     %use default range of coordinates
%     %(only the numerical values below should be changed for different
%     %fitness functions)
%     xVec(validPts,:) = s2rscalar(xVec(validPts,:),-5,5);
% else
%     xVec(validPts,:) = s2rvector(xVec(validPts,:),inParams);
% end

%x=zeros(1,ncols);
% ===============================
%realCoord = zeros(size(xVec));
realCoord = zeros(1,2);
for lpc = 1:nrows
    if validPts(lpc)
    % Only the body of this block should be replaced for different fitness
    % functions
        x = xVec(lpc,:);
        % x(1) alpha; x(2) delta; x(3) omega; x(4) phi0
        % x(5) phi1; x(6) phi2; x(7) phi3; x(8) phi4 ...
        % 0<=x<=1 is convert to physical unit in 'LogLikelihoodRatio'
            %
            [ft,dummy]=svdDSpsoET(x,inParams);  % search sky, earth term only
            fitVal(lpc)= ft;
            realCoord(lpc,:)=dummy;
            %a=dummy1;
    end
end
% ===============================


%Return real coordinates if requested
if nargout > 1
    varargout{1}=realCoord;
end
% END of function