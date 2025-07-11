function K = GEOSCHEM_K_v1452(Met)
% This function calculates the rate constants for the F0AM-compliant
% GEOS-Chem Mechanism version 14.5.2-Hg. This file was generated
% by the GEOS-Chem Emulator on 2025-07-10 from the following KPP file:
%    /uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/GEOSChem_Emulator/Input_Files/GC_Mech_Files/v14.5.2-Hg/Hg.eqn 
% This file contains (only) the definitions of not-commented-out, named
% reactions rates that are referenced within:
%   /uufs/chpc.utah.edu/common/home/u6044586/MATLAB/F0AM_v4.4.2/Chem/GEOSChem/v14.5.2-Hg-reg/GEOSCHEM_GasRxns_v1452.m
% 
% INPUTS:
%   (1) Met: A structure containing environment variables (indexed in time)
%            containing variables used in rate functions.
%
% OUTPUTS:
%    (1) K: structure of rate constants. Each is size length(Met.T) x 1

% All of the following functions are referenced within this file when defining each
% named reaction's rate. Each of these functions must be defined in the
% "import_GC_rates.m" file and imported here for your mechanism to work.
% This import should be automatically done for you, using import_GC_rates()
% below, so this list is just for your reference.
%
%   required_functions={GCARR_ab(),GCARR_abc(),GCJPLPR_abab(),GCARR_ac()};


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%    Import Functions Used in Rate Constants    %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import the MATALB dictionary containing a cell array with the handles to
% individual rate constant functions. In this dictionary, the keys are the
% function names and the  values are the function handles (which are used
% to set rate constants).
rate_dict = import_GC_rates_v1452({'GCARR_ab','GCARR_abc','GCJPLPR_abab','GCARR_ac'});

% Unpack function handles from the dictionary we just returned to use herein.
[GCARR_ab,GCARR_abc,GCJPLPR_abab,GCARR_ac]= rate_dict{:};

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%    Initialize Rate Vars/Functions    %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the total number of rate constants
nk = 26;

% Create Empty matrix of Nans:  len = # of dif conditions (aka time)
%                             width = # of rate constants
krx = nan(length(Met.T),nk);  

% Create an empy cell of "Knames" that will correspond to the column of this rate in krx
%( A hacky way to define a dictionary before MATLAB implemented them in v2022b... )
Knames = cell(nk,1); 

% Extract all vars that might be used directly in rates in this file: 
[TEMP, PRESS, NUMDEN, H2O,TEMP_OVER_K300, K300_OVER_TEMP,INV_TEMP,SR_TEMP,RELHUM] = extract_MetVars(Met);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%   BEGIN NON-HET RATE DEFINITIONS    %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=0; % Set initial counter var for rates...

i=i+1;
Knames{i} = 'Hg0_Br_to_HgBr';
krx(:,i) = GCARR_ab(Met, 1.46e-32, 1.86) * NUMDEN;

i=i+1;
Knames{i} = 'HgBr_to_Hg0';
krx(:,i) = GCARR_abc(Met, 1.6e-9, 1.86, -7801.0) * NUMDEN;

i=i+1;
Knames{i} = 'HgBr_NO2_to_HgBrNO2';
krx(:,i) = GCJPLPR_abab(Met, 4.3e-30, 5.9, 1.2e-10, 1.90, 0.6);

i=i+1;
Knames{i} = 'HgBr_HO2_to_HgBrHO2';
krx(:,i) = GCJPLPR_abab(Met, 4.3e-30, 5.9, 6.9e-11, 2.40, 0.6);

i=i+1;
Knames{i} = 'HgBr_ClO_to_HgBrClO';
krx(:,i) = GCJPLPR_abab(Met, 4.3e-30, 5.9, 6.9e-11, 2.40, 0.6);

i=i+1;
Knames{i} = 'HgBr_BrO_to_HgBrBrO';
krx(:,i) = GCJPLPR_abab(Met, 4.3e-30, 5.9, 6.9e-11, 2.40, 0.6);

i=i+1;
Knames{i} = 'HgBr_OH_to_HgBrOH';
krx(:,i) = GCJPLPR_abab(Met, 4.3e-30, 5.9, 6.9e-11, 2.40, 0.6);

i=i+1;
Knames{i} = 'HgBrO_CH4_to_HgBrOH';
krx(:,i) = GCARR_ac(Met, 4.1e-12, -856.0);

i=i+1;
Knames{i} = 'HgBrO_CO_to_HgBr';
krx(:,i) = GCARR_ac(Met, 6.0e-11, -550.0);

i=i+1;
Knames{i} = 'Hg0_Cl_to_HgCl';
krx(:,i) = GCARR_ac(Met, 2.25e-33, 680.0) * NUMDEN;

i=i+1;
Knames{i} = 'HgCl_NO2_to_HgClNO2';
krx(:,i) = GCJPLPR_abab(Met, 4.3e-30, 5.9, 1.2e-10, 1.90, 0.6);

i=i+1;
Knames{i} = 'HgCl_HO2_to_HgClHO2';
krx(:,i) = GCJPLPR_abab(Met, 4.3e-30, 5.9, 6.9e-11, 2.40, 0.6);

i=i+1;
Knames{i} = 'HgCl_ClO_to_HgClClO';
krx(:,i) = GCJPLPR_abab(Met, 4.3e-30, 5.9, 6.9e-11, 2.40, 0.6);

i=i+1;
Knames{i} = 'HgCl_BrO_to_HgClBrO';
krx(:,i) = GCJPLPR_abab(Met, 4.3e-30, 5.9, 6.9e-11, 2.40, 0.6);

i=i+1;
Knames{i} = 'HgClO_CH4_to_HgClOH';
krx(:,i) = GCARR_ac(Met, 4.1e-12, -856.0);

i=i+1;
Knames{i} = 'HgClO_CO_to_HgCl';
krx(:,i) = GCARR_ac(Met, 6.0e-11, -550.0);

i=i+1;
Knames{i} = 'Hg0_OH_to_HgOH';
krx(:,i) = GCARR_ac(Met, 3.34e-33, 43.0) * NUMDEN;

i=i+1;
Knames{i} = 'HgOH_to_Hg0';
krx(:,i) = GCARR_ac(Met, 1.22e-9, -5720.0) * NUMDEN;

i=i+1;
Knames{i} = 'HgOH_NO2_to_HgOHNO2';
krx(:,i) = GCJPLPR_abab(Met, 4.1e-30, 5.9, 1.2e-10, 1.90, 0.6);

i=i+1;
Knames{i} = 'HgOH_HO2_to_HgOHHO2';
krx(:,i) = GCJPLPR_abab(Met, 4.1e-30, 5.9, 6.9e-11, 2.40, 0.6);

i=i+1;
Knames{i} = 'HgOH_ClO_to_HgOHClO';
krx(:,i) = GCJPLPR_abab(Met, 4.1e-30, 5.9, 6.9e-11, 2.40, 0.6);

i=i+1;
Knames{i} = 'HgOH_BrO_to_HgOHBrO';
krx(:,i) = GCJPLPR_abab(Met, 4.1e-30, 5.9, 6.9e-11, 2.40, 0.6);

i=i+1;
Knames{i} = 'HgOH_Br_to_HgBrOH';
krx(:,i) = GCJPLPR_abab(Met, 4.1e-30, 5.9, 6.9e-11, 2.40, 0.6);

i=i+1;
Knames{i} = 'HgOH_OH_to_HgOHOH';
krx(:,i) = GCJPLPR_abab(Met, 4.1e-30, 5.9, 6.9e-11, 2.40, 0.6);

i=i+1;
Knames{i} = 'HgOHO_CH4_to_HgOHOH';
krx(:,i) = GCARR_ac(Met, 4.1e-12, -856.0);

i=i+1;
Knames{i} = 'HgOHO_CO_to_HgOH';
krx(:,i) = GCARR_ac(Met, 6.0e-11, -550.0);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%    ACCUMULATE INTO STRUCT    %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check that lens are consistent before putting names/values into the output 
% struct, K. 
if length(Knames) ~= length(krx(1,:))
    error(['In GEOSCHEM_K.m, the vars "Knames" and "krx" are ', ...
           'not the same length! Will lead to mismatch ', ...
           'in rates used!']);
end

% Accumulate all the rate names/values into our wanna-be rate dictionary, the
% output struct, K.

K = struct;
for i=1:length(Knames)
    K.(Knames{i}) = krx(:,i);
end

function  [TEMP, PRESS, NUMDEN, H2O,TEMP_OVER_K300, K300_OVER_TEMP,INV_TEMP,SR_TEMP,RELHUM] = extract_MetVars(Met) 
    % This function is used to extract invidividual met vars from the structure Met and is 
    % called at the top of all the rate functions defined in this file so they have access 
    % to these variables to use. These vars defined as "globals" in GEOS-Chem and don't appear 
    % as function input vars, so we just pass the struct "MET" as input to each function and 
    % use this one to then extract them out in each function ... 
	
    TEMP= Met.T;                % Temperature (K)
    PRESS= Met.P;               % Pressure (mbar / hPa) (GC uses hPa)
    NUMDEN= Met.M;              % Number density (molecules/cm3)
    H2O= Met.H2O;               % Water vapor number density (molecules/cm3)
    TEMP_OVER_K300= Met.T./300; % Temp (K) divided by 300 K. 
    K300_OVER_TEMP= 300./Met.T;  % 300K divided by Temperature(K)
    INV_TEMP= 1./Met.T;         % Inverse Temperature (1/K) 
    SR_TEMP = sqrt(Met.T);      % Square root of Temperature (K^0.5)
    RELHUM = Met.RH;            % Relative Humididty as a percentage 
end

end