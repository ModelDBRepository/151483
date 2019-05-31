function [f,Conductances] =  EulerMaruyama(x,dt,ModelSwitch,NoiseSwitches,Params,NumChannelTypes,ActivationVarsPerChannelType,NumActivationVars,StatesPerChannelType,InputCurrent,NumChannelsEachType)

%get activation and inactivation rates (Alphas and Betas) for current membrane potential
%also find the A matrices (returned in a cell array, where each cell is a different channel type type)
%also find the state probabilities (returns in a cell array, p, where each cell is a different channel type)
[A,p,Alphas,Betas] = UpdateDriftAndOccupancies(x,NumActivationVars,StatesPerChannelType,ModelSwitch);

%extract the conductances from the last element of the p vectors
Conductances = zeros(1,NumChannelTypes);
for i  = 1:NumChannelTypes
    Conductances(1,i) = p{i}(end);
end

%update V and activation variables, using their previous values and fluctuations at previous time step
%for these variables, we just use the Euler method
f = x + dt*UpdateEquations(x,p,Params,NumChannelTypes,ActivationVarsPerChannelType,NumActivationVars,StatesPerChannelType,InputCurrent,Alphas,Betas);

%update the associated SDEs using the Euler-Maruyama method
j = NumActivationVars + 1;
for i  = 1:NumChannelTypes
    k = StatesPerChannelType(i);
    if NoiseSwitches(i) == 1
        %DiffusionFunc finds the necessary matrix square roots of the diffusion matrices formed from the A matrices
        f(j+1:j+k) = f(j+1:j+k) + dt * A{i}*x(j+1:j+k) + sqrt(dt) * UpdateDiffusionMatrixSquareRoot(A{i},p{i},NumChannelsEachType(i)) * randn(k,1);
    end
    j = j + k;
end


