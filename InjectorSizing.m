%% Main Injector Sizing
addpath('/Users/JBR132/Documents/_Stanford/AA284B Propulsion System Design Lab/MATLAB-PnID/Gas Feed System Model')
% addpath('D:\Documents\GitHub\MATLAB-PnID\Gas Feed System Model\')
FluidDatabase


PIF = 1.6e6;% Pa, Methane injector feed pressure
TIF = 244;% K, Methane injector feed temperature
rhoIF = PREoS(Methane,"rho",PIF,TIF);% In manifold
IdPF = (PIF-10e5)/PIF;
rhoIF2 = PREoS(Methane,"rho",PIF*(1-IdPF),TIF);% After orifice

PIO = 20*1e5;% Pa, Oxygen injector feed pressure
TIO = 140;% K, Oxygen injector feed temperature
rhoIO = PREoS(Oxygen,"rho",PIO,TIO);% In manifold
rhoIO = rhoIO(2);
IdPO = (PIO-10e5)/PIO;
rhoIO2 = PREoS(Oxygen,"rho",PIO*(1-IdPO),TIO);% After orifice
rhoIO2 = rhoIO2(2);


n_elements = 10;
O_per_element = 2;
mdotF = 0.07945;% kg/s, Methane mass flow rate
mdotO = 0.1986;% kg/s, Oxygen mass flow rate

CdIF = 0.6;% Cd of thin plate orifice
CdIO = 0.6;

AIF = mdotF/(CdIF*sqrt(2*rhoIF*(PIF*IdPF)));% From Cd equation
AIO = mdotF/(CdIO*sqrt(2*rhoIO*(PIO*IdPO)));

% Injector orifice diameters
DIO = sqrt((AIO/n_elements/O_per_element)/pi)*2;% n circular orifices
t = max([PIO*(1-IdPO)*DIO/2*1/convpres(125e3,"psi","Pa"),0.381e-3]);
DIFi = DIO+2*t;
% AIF = pi*(DIFo/2)^2 - pi*(DIFi/2)^2
DIFo = sqrt((AIF/n_elements + pi*(DIFi/2)^2)/pi)*2;

% Injection velocities
VIF = mdotF/(rhoIF2*AIF);
VIO = mdotO/(rhoIO2*AIO);

% Momentum ratio
dp = VIO*mdotO/(VIF*mdotF);

fprintf(['Ox Orifice Diameter: %0.2f mm \n',...
    'Fuel Annulus Inner Diameter: %0.2f mm \n',...
    'Fuel Annulus Outer Diameter: %0.2f mm \n',...
    'Ox Injection Velocity: %0.2f m/s \n',...
    'Fuel Injection Velocity: %0.2f m/s \n',...
    'O/F Momentum Flux Ratio: %0.4f\n'],...
    DIO*1000,DIFi*1000,DIFo*1000,VIO,VIF,dp)

ts = deg2rad(0:1:360);
rs = zeros(length(ts),3);
rs(:,1) = DIO/2;
rs(:,2) = DIFi/2;
rs(:,3) = DIFo/2;
polarplot(ts,rs(:,1))
hold on
polarplot(ts,rs(:,2))
polarplot(ts,rs(:,3))
hold off