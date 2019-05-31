
%Rothman and Manis Type II model
ICs = [-63.6269,0.5122,0.6618,0.0077,0.0011,0.0251,0.4430,0.1458];

C = 12;

gNa=2000;
gLT=400; 
gN=255; 
gP = 45;
gR = 40; 
glk = 4; 

ENa=55;
ELT=-70; 
EN=-70; 
EP = -70;
ER = -43; 
Elk = -65;

gs = [gNa,gLT,gN,gP,gR,glk];
Es = [ENa,ELT,EN,EP,ER,Elk];

NumChannelTypes = 5;
ActivationVarsPerChannel = [2,2,1,1,1];
NumActivationVars = sum(ActivationVarsPerChannel);
NumGatesPerActivationVariable{1} = [3,1];
NumGatesPerActivationVariable{2} = [4,1];
NumGatesPerActivationVariable{3} = [2];
NumGatesPerActivationVariable{4} = [1];
NumGatesPerActivationVariable{5} = [1];

SpikeThreshold = -35;