% Defines LH2 model parameters
% here : (ST) or 1 is trailer (horizontal cylinder)
%        (ET) or 2 is station storage (vertical cylinder)

name = evalin ('base', 'LH2Model.name');
clear LH2Model;
HydrogenTransfer = evalin ('base','HydrogenTransfer');
LH2Model.name = name;

% Constant parameters
psiToPa = 6894.75729;        % conversion factor, psi to Pascals
barToPa = 100000;              % conversion factor, bar to Pascals
galTom3 = 0.00378541;        % conversion factor, from gallons to cubic meters
litersTom3=0.001;            % conversion factor, from liters to cubic meters
inTom = 0.0254;              % conversion factor, from inches to meters
dayTosec = 24*60*60;         % conversion factor, from day to sec
LH2Model.p_atm = 1.01325e5;	 % [Pa] pressure of atmosphere
LH2Model.g = 9.81;			 % [m/s^2] acceleration due to gravity

% Tank geometry parameters
LH2Model.VTotal1 = 18000*litersTom3;          % [m^3]total volume of (ST)    
LH2Model.VTotal2 = 18000*litersTom3;           % [m^3] total volume of (ET)   
LH2Model.R1 = 1; 	         				  % [m] radius of (ST)
LH2Model.R2 = 1;  				          % [m] radius of (ET)
LH2Model.A1 = pi*(LH2Model.R1)^2;	          % [m^2] cross section area of (ST)
LH2Model.A2 = pi*(LH2Model.R2)^2;	          % [m^2] cross section area of (ET)
LH2Model.Lcyl = LH2Model.VTotal1/LH2Model.A1; % [m] length of (ST)
LH2Model.H = LH2Model.VTotal2 /LH2Model.A2;   % [m] height of (ET)

% LH2 parameters
LH2Model.T_c = 32.938;				% [K] critical temperature of hydrogen, updated using REFPROP 9.1
LH2Model.p_c = 186.49 * psiToPa;	% [Pa] critical pressure of hydrogen, updated using REFPROP 9.1
LH2Model.lambda = 5;				% dimensionless exponent for saturated H2 vapor, used to determine film temperature assuming saturated hydrogen vapor
LH2Model.rho_L = 70.9;				% [kg/m^3] liquid density, updated , updated using REFPROP 9.1 @ 1 bar
LH2Model.c_L = 9702.5;				% [J/kg/K] specific heat of liquid, updated , updated using REFPROP 9.1 @ 1 bar
LH2Model.kappa_L = 0.10061;			% [W/mK] thermal conductivity of liquid, updated using REFPROP 9.1 @ 1 bar
LH2Model.mu_L = 13.54e-6;			% [Pa*s] dynamic viscosity (liquid), updated using REFPROP 9.1 @ 1 bar

% GH2 (vapor) parameters
LH2Model.R_v = 4124;				% [J/kg/K] vapor constant
LH2Model.c_v = 6490;				% [J/kg/K] specific heat of hydrogen at V=const; rotational degrees are frozen
LH2Model.c_p = LH2Model.c_v + LH2Model.R_v;
LH2Model.gamma_ = 5/3;				% dimensionless parameter
LH2Model.Gamma_ = ((LH2Model.gamma_+1)/2)^((LH2Model.gamma_+1)/2/(LH2Model.gamma_-1)); % dimensionless parameter
LH2Model.mu_v = 0.98e-6;			% [Pa*s] dynamic viscosity (vapor), updated using REFPROP 9.1 @ 1 bar
LH2Model.kappa_v = 0.0166;			% [W/mK] thermal conductivity of saturated vapor, updated using REFPROP 9.1 @ 1 bar

% Grid parameters
LH2Model.nL1 = 3;					% [] ST grid size (liquid)
LH2Model.tminL1 = 0.1;				% [s] ST grid time constant (liquid)
LH2Model.nL2 = 3;					% [] ET grid size (liquid)
LH2Model.tminL2 = 0.1;				% [s] ET grid time constant (liquid)
LH2Model.nV1 = 3;					% [] ST grid size (vapor)
LH2Model.tminV1 = 0.1;				% [s] ST grid time constant (vapor)
LH2Model.nV2 = 4;					% [] ET grid size (vapor)
LH2Model.tminV2 = 0.1;				% [s] ET grid time constant (vapor)

% Initial conditons (ST)
LH2Model.p10 = 3*barToPa;        % [Pa] initial pressure !! VALUE SHOULD BE BELOW PRESSURE RATING OR CRITICAL PRESSURE (186 psia), WHICHEVER IS LOWEST !!
LH2Model.TL10 = 21;                 % [K] initial liquid temperature
LH2Model.initialfillfraction = 0.9; % [0-1] initial filling fraction
LH2Model.totalmass10 = LH2Model.initialfillfraction*LH2Model.VTotal1*LH2Model.rho_L;        % [kg] inital total mass (liquid + vapor) !! using density that considers only liquid, but this is only an estimate value for total mass.

LH2Model.Tv10 = 0.1+(-1.603941638811E-11*(LH2Model.p10/psiToPa)^6 + 7.830478134841E-09*(LH2Model.p10/psiToPa)^5 - 1.549372675881E-06*(LH2Model.p10/psiToPa)^4 + 1.614567978153E-04*(LH2Model.p10/psiToPa)^3 - 9.861776990784E-03*(LH2Model.p10/psiToPa)^2 + 4.314905904166E-01*(LH2Model.p10/psiToPa)^1 + 1.559843335080E+01); % saturation temperature of vapor, from Refprop
LH2Model.Ts10 = LH2Model.T_c*(LH2Model.p10/LH2Model.p_c)^(1/LH2Model.lambda);      % [K] initial film temperature. From Osipov 2008, see ref in Readme file
LH2Model.rhov10 = refpropm('D','T',LH2Model.Tv10,'P',LH2Model.p10/1000,'PARAHYD'); % [g/L] initial density of the vapor phase in ST
LH2Model.rhoL10 = refpropm('D','T',LH2Model.TL10,'Q',0,'PARAHYD');                 % [g/L] inital density of liquid phase in ST
LH2Model.Vullage10 = (LH2Model.totalmass10- LH2Model.rhoL10*LH2Model.VTotal1)/(LH2Model.rhov10-LH2Model.rhoL10); % [m^3] initial ullage volume in ST
LH2Model.mL10 = LH2Model.rhoL10 * (LH2Model.VTotal1 - LH2Model.Vullage10);         % [kg] initial liquid mass in ST
LH2Model.mv10 = LH2Model.totalmass10 - LH2Model.mL10;                              % [kg] initial vapor mass in ST

% Initial conditions (ET)
LH2Model.p20 = 3*barToPa;             % [Pa] initial pressure !! VALUE SHOULD BE BELOW PRESSURE RATING OR CRITICAL PRESSURE (186 psia), WHICHEVER IS LOWEST !!
LH2Model.Tv20 = 0.1+(-1.603941638811E-11*(LH2Model.p20/psiToPa)^6 + 7.830478134841E-09*(LH2Model.p20/psiToPa)^5 - 1.549372675881E-06*(LH2Model.p20/psiToPa)^4 + 1.614567978153E-04*(LH2Model.p20/psiToPa)^3 - 9.861776990784E-03*(LH2Model.p20/psiToPa)^2 + 4.314905904166E-01*(LH2Model.p20/psiToPa)^1 + 1.559843335080E+01); % saturation temperature of vapor, from Refprop
LH2Model.TL20 = 21;                  % [K] initial liquid temperature
LH2Model.Tw20 = 21;                    % [K] initial wall temperature     
LH2Model.pct_hL20 = 0.01;              % initial level of liquid, measured in inH2O but reported as a fraction (i.e. 0.5 is 5 out of 10 inH2O). Value should be between >0 and 1 (=0 may trigger errors)

LH2Model.Ts20 = LH2Model.T_c*(LH2Model.p20/LH2Model.p_c)^(1/LH2Model.lambda);      % [K] initial film temperature.  From Osipov 2008, see reference in Readme file
LH2Model.hL20 = LH2Model.pct_hL20 * LH2Model.H;                                    % [m] initial level of liquid in ET
LH2Model.VL20 = LH2Model.hL20 * LH2Model.A2;                                       % [m^3] inital volume of liquid in ET
LH2Model.rhoL20 = refpropm('D','T',LH2Model.TL20,'Q',0,'PARAHYD');                 % [g/L] inital density of liquid phase in ET
LH2Model.mL20 = LH2Model.rhoL20 *LH2Model.VL20;                                    % [kg] initial mass of liquid in ET
LH2Model.rhov20 = refpropm('D','T',LH2Model.Tv20,'P',LH2Model.p20/1000,'PARAHYD'); % [g/L] initial density of the vapor phase in ET
LH2Model.mv20 = LH2Model.rhov20 * (LH2Model.VTotal2 - LH2Model.VL20);              % [kg] initial mass of vapor in ET

% Initial flows
LH2Model.Jboil0 = 0;                   % [kg/s] initial boiling flow of vaporizer
LH2Model.Jtr0 = 0;                     % [kg/s] initial transmission line flow

% Heat transfer coefficients (ST)
LH2Model.QdotEL1 = 200*(LH2Model.VTotal1/64.5);                % [W] heat transfer from environment to liquid. Using initial value as a reference, scaled to actual volume.
LH2Model.QdotEV1 = 40*(LH2Model.VTotal1/64.5);                 % [W] heat transfer from environment to vapor. Using initial value as a reference, scaled to actual volume.
LH2Model.STBoiloffrate = 0.015/dayTosec;                       % [%] flowrate as a %/s of the total liquid mass.

% Heat transfer coefficients (ET)
LH2Model.QdotEW2 = -7.462776654302E-02*LH2Model.VTotal2^2 + 4.445867251697E+00*LH2Model.VTotal2 + 3.108170556297E+01;  % [W] heat transfer from environment to vessel's wall depending on the volume. Correlation from LLNL.
LH2Model.mw2 = LH2Model.VTotal2*5000/12.5;       % [kg] mass of ET inner vessel. Using as a reference: tank of 12500L weights 5000kg
LH2Model.ETBoiloffrate = 0.015/dayTosec;                       % [%] flowrate as a %/s of the total liquid mass.

% Vaporizer parameters, in (ST)
LH2Model.mVap0 = 0;                    % [kg] initial liquid mass in vaporizer tubes
LH2Model.Tboil = 22;                   % [K] temperature of vapor bubbles
LH2Model.tau_vap = 0.3;                % [s] vaporizer time constant !! MAY NEED TO INCREASE THIS VALUE IF ODE SOLVER CRASHES (TOO STIFF) !!
LH2Model.c_vap = 4e-4;                 % [] vaporizer valve flow coefficient
LH2Model.VapValveState = 0;            % initial vaporizer valve state

% Transmission line parameters, between (ST) and (ET)
LH2Model.DPipe = 2*inTom;              % [m] diameter of pipe
LH2Model.LPipe = 10;                   % [m] transmission line length
LH2Model.drPipe = 1e-6;                % roughness of pipe
%LH2Model.ReCrit = 3e3;                 % critical Re value
LH2Model.f = 1.3/log(LH2Model.DPipe/2/LH2Model.drPipe)^2; % friction parameter
LH2Model.tau_tr = 5;                  % [s] transmission line delay constant !! MAY NEED TO INCREASE THIS VALUE IF ODE SOLVER CRASHES (TOO STIFF) !!
LH2Model.dE = 0.4*inTom;		       % [m] transfer line valve diameter (10 in.)
LH2Model.kE = 4;				       % [] transfer line valve coefficient

% Top fill (ignore if deactivated)
if Topfill
    LH2Model.ETinletdiameter = 0.02;         % [m] Diameter of inlet nozzle
    LH2Model.ETnozzleamout=1;
    LH2Model.ConvCoeffTopfill = 000;           % Number of nozzles to simulate spray
    LH2Model.initial_ratio_top_bottom=0.05*(LH2Model.ConvCoeffTopfill/20000);        % ratio between top and bottom fill to (ET) at the beginning of the process. 0.5 = 50% of liquid goes to top. Values between 0.01 and 0.03 are best...
else
    LH2Model.initial_ratio_top_bottom=0;
    LH2Model.ConvCoeffTopfill = 0;           % Number of nozzles to simulate spray
    LH2Model.bulkevap_ratio_top_bottom=0;
end

% Vent valves
LH2Model.S_valve1 = 0.001;              % [m^2] orifice area of ST vent valve !! MAY NEED TO REDUCE THIS VALUE IF ODE SOLVER CRASHES (TOO STIFF) !!
LH2Model.STVentState = 0;               % initial ST vent state (starts closed)
LH2Model.S_valve2 = 3.1416*(2.5*0.0035)^2; % [m^2] orifice area of ET vent valve !! MAY NEED TO REDUCE THIS VALUE IF ODE SOLVER CRASHES (TOO STIFF) !!
LH2Model.ETVentState = 0;               % initial ET vent state (starts closed) 

% Pressure settings
LH2Model.p_ST_slow = 12*barToPa;       % [Pa] threshold pressure for slow fill 
LH2Model.p_ST_fast = 12*barToPa;       % [Pa] threshold pressure for fast fill
LH2Model.p_ST_final = 12*barToPa;      % [Pa] final venting pressure for (ST) (= pressure in the trailer before leaving the station)

LH2Model.p_ET_low = 9*barToPa;        % [Pa] ET vent valve lower pressure threshold (= PRD hysteresis pressure)
LH2Model.p_ET_high = 10*barToPa;       % [Pa] ET vent valve upper pressure threshold (= PRD set pressure) 
LH2Model.p_ET_final = 10*barToPa;      % [Pa] final venting pressure for (ET) (= pressure in the truck's tank before leaving the station)


LH2Model.TopET = 0.9;                  % [] maximum fraction full for ET, when fill stops. Value should be between 0 and 1 (0.9 = 90%)

% Solver options
if HydrogenTransfer==0      % If hydrogen transfer is deactivated, duration is set
    LH2Model.tFinal = 1*60*60;
    LH2Model.namedetail=name+"_"+string(LH2Model.VTotal1)+"-"+string(LH2Model.VTotal2)+"m3";
else  % If hydrogen transfer is activation, simulation duration is approximate to the transfer process duration
    LH2Model.tFinal = 0.7*0.75*0.75*60*60*(LH2Model.VTotal2/12.5)*(LH2Model.TopET-LH2Model.pct_hL20)/1.25;  % [s] Simulation final time (adaptative with ET volume and fill percentage)
    LH2Model.namedetail=name+"_"+string(LH2Model.VTotal1)+"-"+string(LH2Model.VTotal2)+"m3_"+string(LH2Model.p_ST_fast/barToPa)+"-"+string(LH2Model.p_ET_high/barToPa)+"bar";
end

LH2Model.relTol = 5e-5;                % relative tolerance for ode solver. Doesnt work well with 1e-3
