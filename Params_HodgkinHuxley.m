
%Hodgkin and Huxley model
ICs = [0,0.06,0.6,0.35];

C = 1; % muF /cm^2
gNa = 120; % mS/cm^2
ENa = 115; % mV
gK = 36; % mS/cm^2
EK = -12; % mV
gL = 0.3; % mS / cm^2
EL = 10.0; % mV

gs = [gNa,gK,gL];
Es = [ENa,EK,EL];

NumChannelTypes = 2;
ActivationVarsPerChannel = [2,1];
NumActivationVars = sum(ActivationVarsPerChannel);
NumGatesPerActivationVariable{1} = [3,1];
NumGatesPerActivationVariable{2} = [4];

SpikeThreshold = 50;