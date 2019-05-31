function [A,X,Alphas,Betas] = UpdateDriftAndOccupancies(x,NumActivationVars,StatesPerChannelType,Switch)

Alphas = zeros(NumActivationVars,1);
Betas = zeros(NumActivationVars,1);

NumChannelTypes = length(StatesPerChannelType);
A = cell(1,NumChannelTypes);
X = cell(1,NumChannelTypes);
for i = 1:NumChannelTypes
    A{i} = zeros(StatesPerChannelType(i),StatesPerChannelType(i));
    X{i} = zeros(StatesPerChannelType(i),1);
end


if Switch == 1
    %Hodgkin-Huxley model:
    
    V = x(1);
    
    Alphas(1) = (25-V)/(10*(exp((25-V)/10)-1));
    Betas(1) = 4*exp(-V/18);
    
    %h
    Alphas(2) = 7*exp(-(V/20))/100;
    Betas(2) = 1/(exp((30-V)/10)+1);
    
    %n
    Alphas(3) = (10-V)/(100*(exp((10-V)/10)-1));
    Betas(3) = exp(-(V/80))/8;
    
    m0 = x(2);
    h0 = x(3);
    n0 = x(4);
    NaBar = [(1-m0)^3*(1-h0) ; 3*(1-m0)^2*m0*(1-h0) ; 3*(1-m0)*m0^2*(1-h0) ; m0^3*(1-h0) ; (1-m0)^3*h0 ; 3*(1-m0)^2*m0*h0 ; 3*(1-m0)*m0^2*h0 ; m0^3*h0];
    KBar  = [(1-n0)^4 ; 4*n0*(1-n0)^3 ; 6*n0^2*(1-n0)^2  ; 4*n0^3*(1-n0)  ; n0^4];
            
    % Drift Na
    ANa =  [ -3*Alphas(1)-Alphas(2)       , Betas(1)                , 0              , 0                      ,  Betas(2)                 , 0                       , 0          , 0 ;
        3*Alphas(1)               ,-2*Alphas(1)-Betas(1)-Alphas(2), 2*Betas(1)        , 0                      ,  0                     , Betas(2)                   , 0          , 0 ;
        0                      , 2*Alphas(1)             , -Alphas(1)-2*Betas(1)-Alphas(2),  3*Betas(1)        ,  0                     ,  0                      , Betas(2)      , 0 ;
        0                      , 0                    , Alphas(1)         , -3*Betas(1)-Alphas(2)        ,  0                     ,  0                      , 0          , Betas(2)    ;
        Alphas(2)                 , 0                    , 0              , 0                      ,  -3*Alphas(1) - Betas(2)     , Betas(1)                   , 0          , 0 ;
        0                      , Alphas(2)               , 0              , 0                      ,  3*Alphas(1)              ,  -2*Alphas(1)-Betas(1)-Betas(2)  ,   2*Betas(1)  , 0 ;
        0                      , 0                    , Alphas(2)         , 0                      ,  0                     ,  2*Alphas(1)               ,   -Alphas(1)-2*Betas(1)-Betas(2) , 3*Betas(1)  ;
        0                      , 0                    , 0              , Alphas(2)                 ,           0            ,  0                      ,  Alphas(1)    , -3*Betas(1)-Betas(2)];
    
    % Drift K
    AK = [-4*Alphas(3), Betas(3)             , 0                , 0                  , 0;
        4*Alphas(3), -3*Alphas(3)-Betas(3),  2*Betas(3)               , 0,                   0;
        0,        3*Alphas(3),        -2*Alphas(3)-2*Betas(3), 3*Betas(3),          0;
        0,        0,               2*Alphas(3),          -Alphas(3)-3*Betas(3), 4*Betas(3);
        0,        0,               0,                 Alphas(3),          -4*Betas(3)];
    
    A{1}=ANa;
    A{2} = AK;
    X{1} = NaBar;
    X{2} = KBar;
else
    %Rothman and Manis model
    
    V = x(1);
    
    winf = (1+exp(-(V+48)/6))^(-1/4);
    tw = ((100*(6*exp((V+60)/6)+16*exp(-(V+60)/45))^(-1)+1.5)/3);
    
    zinf = 0.5*(1+exp((V+71)/10))^(-1)+0.5;
    tz = ((1000*(exp((V+60)/20)+exp(-(V+60)/8))^(-1)+50)/3);
    
    ninf = (1+exp(-(V+15)/5))^(-1/2);
    tn = ((100*(11*exp((V+60)/24)+21*exp(-(V+60)/23))^(-1)+0.7)/3);
    
    pinf = (1+exp(-(V+23)/6))^(-1);
    tp =(( 100*(4*exp((V+60)/32)+5*exp(-(V+60)/22))^(-1)+5)/3);
    
    minf = (1+exp(-(V+38)/7))^(-1);
    tm =(( 10*(5*exp((V+60)/18)+36*exp(-(V+60)/25))^(-1)+0.04)/3);
    
    hinf = (1+exp((V+65)/6))^(-1);
    th =(( 100*(7*exp((V+60)/11)+10*exp(-(V+60)/25))^(-1)+0.6)/3);
    
    rinf = (1+exp((V+76)/7))^(-1);
    tr =(( (10^5)*(237*exp((V+60)/12)+21*exp(-(V+60)/23))^(-1)+0.7)/3);
    
    
    Alphas(1) = winf/tw;
    Alphas(2) = zinf/tz;
    Alphas(3) = ninf/tn;
    Alphas(4) = pinf/tp;
    Alphas(5) = minf/tm;
    Alphas(6) = hinf/th;
    Alphas(7) = rinf/tr;
    
    Betas(1) = (1-winf)/tw;
    Betas(2) = (1-zinf)/tz;
    Betas(3) = (1-ninf)/tn;
    Betas(4) = (1-pinf)/tp;
    Betas(5) = (1-minf)/tm;
    Betas(6) = (1-hinf)/th;
    Betas(7) = (1-rinf)/tr;
    
    w = x(2);
    z = x(3);
    n = x(4);
    p = x(5);
    m = x(6);
    h = x(7);
    r = x(8);
    
    %State occupancy probabilities
    NaBar = [(1-m)^3*(1-h) ; 3*(1-m)^2*m*(1-h) ; 3*(1-m)*m^2*(1-h) ; m^3*(1-h) ; (1-m)^3*h ; 3*(1-m)^2*m*h ; 3*(1-m)*m^2*h ; m^3*h];
    LTBar=[((1-w)^4)*(1-z); 4*((1-w)^3)*w*(1-z); 6*((1-w)^2 )*(w^2)*(1-z); 4*((1-w)*(w^3)*(1-z)); (1-z)*w^4;((1-w)^4)*z; 4*((1-w)^3)*w*z; 6*((1-w)^2 *(w^2)*z); 4*((1-w)*(w^3)*z); (z*w^4)];
    NBar=[(1-n)^2; 2*(1-n)*n;n^2];
    PBar=[(1-p);p];
    RBar=[(1-r);r];
    
    % Drift Na
    ANa = [ -3*Alphas(5)-Alphas(6)       , Betas(5)                , 0              , 0                      ,  Betas(6)                 , 0                       , 0          , 0 ;
        3*Alphas(5)               ,-2*Alphas(5)-Betas(5)-Alphas(6), 2*Betas(5)        , 0                      ,  0                     , Betas(6)                   , 0          , 0 ;
        0                      , 2*Alphas(5)             , -Alphas(5)-2*Betas(5)-Alphas(6),  3*Betas(5)        ,  0                     ,  0                      , Betas(6)      , 0 ;
        0                      , 0                    , Alphas(5)         , -3*Betas(5)-Alphas(6)        ,  0                     ,  0                      , 0          , Betas(6)    ;
        Alphas(6)                 , 0                    , 0              , 0                      ,  -3*Alphas(5) - Betas(6)     , Betas(5)                   , 0          , 0 ;
        0                      , Alphas(6)               , 0              , 0                      ,  3*Alphas(5)              ,  -2*Alphas(5)-Betas(5)-Betas(6)  ,   2*Betas(5)  , 0 ;
        0                      , 0                    , Alphas(6)         , 0                      ,  0                     ,  2*Alphas(5)               ,   -Alphas(5)-2*Betas(5)-Betas(6) , 3*Betas(5)  ;
        0                      , 0                    , 0              , Alphas(6)                 ,           0            ,  0                      ,  Alphas(5)    , -3*Betas(5)-Betas(6)];
    
    %Drift N
    AN=[-2*Alphas(3), Betas(3),0; 2*Alphas(3), -(Alphas(3)+Betas(3)),2*Betas(3); 0, Alphas(3), -2*Betas(3)];
    
    %Drift P
    AP=[-Alphas(4),Betas(4); Alphas(4),-Betas(4)];
    
    %Drift R
    AR=[-Alphas(7),Betas(7); Alphas(7),-Betas(7)];
    
    %Drift LT
    ALT = zeros(StatesPerChannelType(2),StatesPerChannelType(2));
    ALT(1,1)=  -4*Alphas(1)-Alphas(2);
    ALT(1,2)=Betas(1);
    ALT(2,1)=4*Alphas(1);
    ALT(1,6)=Betas(2);
    ALT(6,1)=Alphas(2);
    
    ALT(2,2)=-3*Alphas(1)-Betas(1)-Alphas(2);
    ALT(2,3)=2*Betas(1);
    ALT(3,2)=3*Alphas(1);
    ALT(2,7)=Betas(2);
    ALT(7,2)=Alphas(2);
    
    ALT(3,3)=-2*Alphas(1)-2*Betas(1)-Alphas(2);
    ALT(3,4)=3*Betas(1);
    ALT(4,3)=2*Alphas(1);
    ALT(3,8)=Betas(2);
    ALT(8,3)=Alphas(2);
    
    ALT(4,4)=-Alphas(1)-3*Betas(1)-Alphas(2);
    ALT(4,5)=4*Betas(1);
    ALT(5,4)=Alphas(1);
    ALT(4,9)=Betas(2);
    ALT(9,4)=Alphas(2);
    
    ALT(5,5)=-4*Betas(1) - Alphas(2);
    ALT(5,10)=Betas(2);
    ALT(10,5)=Alphas(2);
    
    ALT(6,6)=-4*Alphas(1) - Betas(2);
    ALT(6,7)=Betas(2);
    ALT(7,6)=4*Alphas(1);
    
    ALT(7,7)=-3*Alphas(1)-Betas(1)-Betas(2);
    ALT(7,8)=2*Betas(2);
    ALT(8,7)=3*Alphas(1);
    
    ALT(8,8)=-2*Betas(1)-Betas(2)-2*Alphas(1);
    ALT(8,9)=3*Betas(2);
    ALT(9,8)=2*Alphas(1);
    
    ALT(9,9)=-3*Betas(1)-Betas(2)-Alphas(1);
    ALT(9,10)=4*Betas(2);
    ALT(10,9)=Alphas(1);
    
    ALT(10,10)=-4*Betas(1)-Betas(2);
    
    A{1} = ANa;
    A{2} = ALT;
    A{3} = AN;
    A{4} = AP;
    A{5} = AR;
    
    X{1} = NaBar;
    X{2} = LTBar;
    X{3} = NBar;
    X{4} = PBar;
    X{5} = RBar;
    
end