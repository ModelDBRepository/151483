function [y]  = UpdateEquations(x,p,Params,NumChannelTypes,ActivationVarsPerChannelType,NumActivationVars,StatesPerChannelType,InputCurrent,Alphas,Betas)

gs=Params{1};
Es=Params{2};
C=Params{3};

V_old = x(1);
Activations_old = x(2:NumActivationVars+1);

V_rhs = 0;
Activation_rhs = zeros(NumActivationVars,1);
ActivationIndex = 1;
Ind = 1+NumActivationVars;
for i = 1:NumChannelTypes
    Ind = Ind+StatesPerChannelType(i);
    %get mean and fluctuations of fractions of open channel of this type
    Fracs_old = p{i}(StatesPerChannelType(i));
    Flucts_old = x(Ind);
    
    %update the membrane potential to include this ionic current
    %ensure the total conductance is bounded to the interval [0,gs(i)]
    V_rhs = V_rhs - gs(i)*(max(0,min(1,Fracs_old + Flucts_old))) *(V_old-Es(i));
    
    %update the activation variables for this ionic current
    for j = 1:ActivationVarsPerChannelType(i)
        Activation_rhs(ActivationIndex) = Alphas(ActivationIndex)*(1-Activations_old(ActivationIndex))-Betas(ActivationIndex)*Activations_old(ActivationIndex);
        ActivationIndex = ActivationIndex + 1;
    end 
end
%update with input current, leak current and capacitance
V_rhs = (V_rhs + InputCurrent - gs(end)*(V_old-Es(end)))/C;

%update membrane potential and activation variables
y = zeros(size(x));
y(1) = V_rhs;
y(2:NumActivationVars+1) = Activation_rhs;


