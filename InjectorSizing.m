%% Main Injector Sizing
addpath('/Users/JBR132/Documents/_Stanford/AA284B Propulsion System Design Lab/MATLAB-PnID/')
FluidDatabase

IdP = 0.2;

PIF = 10/(1-IdP)*1e5;% Pa, Methane injector feed pressure
TIF = 230;% K, Methane injector feed temperature
rhoIF = PREoS(Methane,"rho",PIF,TIF);% In manifold
rhoIF2 = PREoS(Methane,"rho",PIF*(1-IdP),TIF);% After orifice

PIO = 10/(1-IdP)*1e5;% Pa, Oxygen injector feed pressure
TIO = 120;% K, Oxygen injector feed temperature
rhoIO = PREoS(Oxygen,"rho",PIO,TIO);% In manifold
rhoIO = rhoIO(2);
rhoIO2 = PREoS(Oxygen,"rho",PIO*(1-IdP),TIO);% After orifice
rhoIO2 = rhoIO2(2);

n_elements = 6;
mdotF = 0.108140110052637;% kg/s, Methane mass flow rate
mdotO = 0.270350275131593;% kg/s, Oxygen mass flow rate

CdIF = 0.6;% Cd of thin plate orifice
CdIO = 0.6;

AIF = mdotF/(CdIF*sqrt(2*rhoIF*(PIF-0.8*PIF)));% From Cd equation
AIO = mdotF/(CdIO*sqrt(2*rhoIO*(PIO-0.8*PIO)));

% Injector orifice diameters
DIO = sqrt((AIO/n_elements)/pi)*2;% n circular orifices
t = max([PIO*(1-IdP)*DIO/2*1/convpres(125e3,"psi","Pa"),0.4e-3]);
DIFi = DIO+2*t;
% AIF = pi*(DIFo/2)^2 - pi*(DIFi/2)^2
DIFo = sqrt((AIF/n_elements + pi*(DIFi/2)^2)/pi)*2;

% Injection velocities
VIF = mdotF/(rhoIF2*AIF);
VIO = mdotO/(rhoIO2*AIO);

% Momentum ratio
dp = VIF*mdotF/(VIO*mdotO);

fprintf(['Ox Orifice Diameter: %0.2f mm \n',...
    'Fuel Annulus Inner Diameter: %0.2f mm \n',...
    'Fuel Annulus Outer Diameter: %0.2f mm \n',...
    'Ox Injection Velocity: %0.2f m/s \n',...
    'Fuel Injection Velocity: %0.2f m/s \n',...
    'F/O Momentum Ratio: %0.4f\n'],...
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