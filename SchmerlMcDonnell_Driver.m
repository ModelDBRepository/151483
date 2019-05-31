clear all

%selection of deterministic or stochastic
IsDeterministic = 0;%set this to 1 for the SSE model, 0 for the deterministic model

%selection of model: only one of the following should be non-zero
IsHodgkin = 1;
IsRothmanII = 0;
IsRothmanI_II = 0;
IsRothmanI_C = 0;

%get parameters for selected model
if IsHodgkin
    %load parameters for Hodgkin-Huxley model
    Switch = 1;
    Params_HodgkinHuxley
    InputCurrent = 7.2; %can be modified as desired
elseif IsRothmanII
    %load parameters for Rothman-Manis Type II model
    Switch = 2;
    Params_RothmanManisTypeII
    InputCurrent = 230; %can be modified as desired
elseif IsRothmanI_II
    %load parameters for Rothman-Manis Type I-II model
    Switch = 2;
    Params_RothmanManisTypeI_II
    InputCurrent = 231; %can be modified as desired
elseif IsRothmanI_C
    %load parameters for Rothman-Manis Type I-C model
    Switch = 2;
    Params_RothmanManisTypeI_C
    InputCurrent = 25.75; %can be modified as desired
end

%the user can decide to make some channel types stochastic and some deterministic if desired
NoiseSwitches = zeros(NumChannelTypes,1);
if IsDeterministic
    %deterministic model
    NoiseSwitches(:) = 0;
else
    %SSE model
    %for any channels that are to be deterministic, set the corresponding index to 0
    NoiseSwitches(:) = 1;
end

%simulation times, initial and input conditions
NumPoints = 1; %this is how many different channel numbers we want to run the sim for
NumRepeats = 1;
%NumChannelsToTry = round(logspace(3,8,NumPoints));
NumChannelsToTry = [10000];
dt = 0.01;
Min_t = 0;
Max_t = 400;
ts = [Min_t:dt:Max_t];
Num_t = length(ts);
OnsetTime = 0; %ms
OffsetTime = 0; %ms

%count the number of variables
StatesPerChannelType = zeros(NumChannelTypes,1);
TotalStates = 0;
for j = 1:NumChannelTypes
    NumStatesPerActivationVariable{j} = NumGatesPerActivationVariable{j} +1;
    StatesPerChannelType(j) = prod(NumStatesPerActivationVariable{j});
    TotalStates = TotalStates + StatesPerChannelType(j);
end
Num_Vars = 1 + NumActivationVars + TotalStates;
Params{1} = gs;
Params{2} = Es;
Params{3} = C;

%set up the input constant current, with onset and offset times
InputCurrents = zeros(Num_t,1);
for i = 2:Num_t
    if ts(i) < OnsetTime | ts(i) > Max_t-OffsetTime
        InputCurrents(i) = 0;
    else
        InputCurrents(i) = InputCurrent;
    end
end
    
%do the simulation for the desired number of repeats and number of channels
SpikeCount = zeros(NumPoints,NumRepeats);
InRefrac = 0;
TotalISICount = 0;
ISIs = zeros(30*NumRepeats,1);
for k = 1:NumRepeats
    for j = 1:NumPoints
        
        Solution = zeros(Num_Vars,Num_t);
        Conductances = zeros(NumChannelTypes,Num_t);
        Solution(1:1+NumActivationVars,1) = ICs;
        
        %set the number of channels of each type
        %can set arbitrarily different numbers for each type
        %could be better to do this in the params file
        NumChannelsEachType(1:NumChannelTypes) = NumChannelsToTry(j);            
        
        %do the simulation
        LastSpikeTime = -10^5;
        for i = 2:Num_t
                        
            %solve the equations
            [Solution(:,i),Conductances(:,i)] = EulerMaruyama(Solution(:,i-1),dt,Switch,NoiseSwitches,Params,NumChannelTypes,ActivationVarsPerChannel,NumActivationVars,StatesPerChannelType,InputCurrents(i),NumChannelsEachType);
            
            %count spikes, using a rudimentary spike detector
            if InRefrac == 0 & Solution(1,i) > SpikeThreshold & Solution(1,i-1) < SpikeThreshold
                SpikeCount(j,k) = SpikeCount(j,k) + 1;
                InRefrac = 1;
                if LastSpikeTime > 0
                    TotalISICount = TotalISICount+1;
                    ISIs(TotalISICount) = ts(i)-LastSpikeTime;
                end
                LastSpikeTime = ts(i);
            end
            if InRefrac == 1 & Solution(1,i) < SpikeThreshold-5
                InRefrac = 0;
            end
        end
    end
end

%plot the membrane potential as a function of time
f0= figure;hold on
plot(ts,Solution(1,:))
box on
xlabel('Time (ms)','fontsize',14)
ylabel('Membrane potential (mV)','fontsize',14)
