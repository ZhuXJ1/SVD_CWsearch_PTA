function returnData=ptapso(fitfuncHandle,nDim)
%% Simple PSO minimizer for PTA data analysis
% S=PTAPSO(Fhandle,N)
% Runs PSO on the fitness function with handle Fhandle. If Fname is the name
% of the function, Fhandle = @(x) <Fname>(x, FP). The optional
% parameters to supply to the fitness function should be given in the structure FP.
% N is the dimensionality of the  fitness function. If N is empty, then N=2. 
% The output is returned in the structure S. The field of S are:
%  bestLocation : Best location found (in standardized coordinates)
%  bestSNR: Best fitness value found
%  totalFuncEvals: Total number of fitness function evaluations
%% Authors
% Soumya D. Mohanty, Sep 2014
% Lightweight version of an existing code called PSO (author: S. D.
% Mohanty) created for a project on PTA data analysis involving SDM and Y.
% Wang. The use of this code outside of this project requires the permission
% of the author (SDM).
% Used in Zhu et al. MNRAS (2016) more recently.
%% References
% Eberhart R. C., Kennedy J., 1995, in Proceedings of the sixth international symposium on micro machine and human science, Vol. 1, New York, NY, pp. 39-43
% Wang Y., Mohanty S. D., 2010, Phys. Rev. D, 81, 063002
% Wang Y., Mohanty S. D., Jenet F. A., 2014, ApJ, 795, 96
% Wang Y., Mohanty S. D., Jenet F. A., 2015, ApJ, 815, 125
% Mohanty S. D. 2012 Astronomical Review, Volume 7, p. 29-35
% Mohanty S. D. 2012 Astronomical Review, Volume 7, Issue 4, Page 4-25
% Zhu et al. 2016, MNRAS
popsize=100;
maxSteps= 2000;
c1=2;
c2=2;
max_initial_velocity = 0.5;
max_velocity = 0.2;
dcLaw_a = 0.9;
dcLaw_b = 0.4;
dcLaw_c = maxSteps;
dcLaw_d = 0.2;
nmOptions = optimset('fminsearch');
nmOptions.MaxFunEvals=200;

partCoordCols = 1:nDim;
partVelCols = (nDim+1):2*nDim;
partPbestCols = (2*nDim+1):3*nDim;
partSnrPbestCols = 3*nDim+1;
partSnrCurrCols = partSnrPbestCols+1;
partSnrLbestCols = partSnrCurrCols+1;
partInertiaCols = partSnrLbestCols+1;
partLocalBestCols = (partInertiaCols+1):(partInertiaCols+nDim);
partFlagFitEvalCols = partLocalBestCols(end)+1;
partFitEvalsCols = partFlagFitEvalCols+1;

nColsPop = length([partCoordCols,partVelCols,partPbestCols,partSnrPbestCols,...
                   partSnrCurrCols,partSnrLbestCols,partInertiaCols,partLocalBestCols,...
                   partFlagFitEvalCols,partFitEvalsCols]);
pop=zeros(popsize,nColsPop);

best_in_history = inf;       
gbest_loc_in_history = 2*ones(1,nDim);
best_fitness = inf;
pop(:,partCoordCols)=rand(popsize,nDim);
pop(:,partVelCols)= -pop(:,partCoordCols)+rand(popsize,nDim);
pop(:,partPbestCols)=pop(:,partCoordCols);
pop(:,partSnrPbestCols)= inf;
pop(:,partSnrCurrCols)=0;
pop(:,partSnrLbestCols)= inf;
pop(:,partLocalBestCols) = 0;
pop(:,partFlagFitEvalCols)=1;
pop(:,partInertiaCols)=0;
pop(:,partFitEvalsCols)=0;
bestLocHistPart = 1;

fitnessValues = zeros(popsize,1);
for lpc_steps=2:maxSteps    
    fitnessValues = fitfuncHandle(pop(:,partCoordCols));
    for k=1:popsize
        startCoord = pop(k,partCoordCols);
        endFitnessVal = pop(k,partSnrCurrCols);
        computeOK = pop(k,partFlagFitEvalCols);
        funcCount = 0;
        if computeOK
            endFitnessVal = fitnessValues(k);
            endCoord = startCoord;
            funcCount = 1;
        else
            endCoord = startCoord;
            funcCount = 0;
        end
        pop(k,partCoordCols)=endCoord;
        pop(k,partSnrCurrCols)=endFitnessVal;
        pop(k,partFitEvalsCols)=pop(k,partFitEvalsCols)+funcCount;
        if pop(k,partSnrPbestCols) > pop(k,partSnrCurrCols)
            pop(k,partSnrPbestCols) = pop(k,partSnrCurrCols);
            pop(k,partPbestCols) = pop(k,partCoordCols);
        end
    end
    [best_fitness,bestfitParticle]=min(pop(:,partSnrCurrCols)); 
    if best_in_history > best_fitness
        best_in_history     =   best_fitness;
        gbest_loc_in_history  =  pop(bestfitParticle,partCoordCols);
        startCoord = gbest_loc_in_history;
        [endCoord,endFitnessVal,~,nmInfo] = fminsearch(fitfuncHandle,...
            startCoord,...
            nmOptions);
        funcCount = nmInfo.funcCount;
        pop(bestfitParticle,partCoordCols)=endCoord;
        pop(bestfitParticle,partSnrCurrCols)=endFitnessVal;
        pop(bestfitParticle,partFitEvalsCols)=pop(bestfitParticle,partFitEvalsCols)+funcCount;
        gbest_loc_in_history = endCoord;
        best_in_history = endFitnessVal;
    else
    end
    for k = 1:popsize
           switch k
               case 1
                   ringNbrs = [popsize,k,k+1];
               case popsize
                   ringNbrs = [k-1,k,1];
               otherwise
                   ringNbrs = [k-1,k,k+1];
           end
           [~,lbestPart] = min(pop(ringNbrs,partSnrCurrCols));
           lbestTruIndx = ringNbrs(lbestPart);
           lbestCurrSnr = pop(lbestTruIndx,partSnrCurrCols);
           if lbestCurrSnr < pop(k,partSnrLbestCols)
               pop(k,partSnrLbestCols) = lbestCurrSnr;
               pop(k,partLocalBestCols) = pop(lbestTruIndx,partCoordCols);
           end
    end
     inertiaWt = max(dcLaw_a-(dcLaw_b/dcLaw_c)*(lpc_steps-1),dcLaw_d);
    for k=1:popsize
        pop(k,partInertiaCols)=inertiaWt;
        partInertia = pop(k,partInertiaCols);
        chi1 = diag(rand(1,nDim));
        chi2 = diag(rand(1,nDim));
        pop(k,partVelCols)=partInertia*pop(k,partVelCols)+...
            c1*(pop(k,partPbestCols)-pop(k,partCoordCols))*chi1+...
            c2*(pop(k,partLocalBestCols)-pop(k,partCoordCols))*chi2;
        maxvBustCompPos = find(pop(k,partVelCols) > max_velocity);
        maxvBustCompNeg = find(pop(k,partVelCols) < -max_velocity);
        if ~isempty(maxvBustCompPos)
            pop(k,partVelCols(maxvBustCompPos))= max_velocity;
        end
        if ~isempty(maxvBustCompNeg)
            pop(k,partVelCols(maxvBustCompNeg))= -max_velocity(1);
        end
        pop(k,partCoordCols)=pop(k,partCoordCols)+pop(k,partVelCols);        
        if any(pop(k,partCoordCols)> 1 | ...
                pop(k,partCoordCols)< 0)
                    pop(k,partSnrCurrCols)= inf;
                    pop(k,partFlagFitEvalCols)= 0;
        else
            pop(k,partFlagFitEvalCols)=1;
        end
    end
end

stepsToEnd=lpc_steps-1;
actualEvaluations = sum(pop(:,partFitEvalsCols));
returnData = struct('totalSteps',stepsToEnd,...
    'totalFuncEvals',actualEvaluations,...
    'bestLocation',gbest_loc_in_history,...
    'bestSNR',best_in_history);





