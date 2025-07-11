function K = GEOSCHEM_K_v1452(Met)
% This function calculates the rate constants for the F0AM-compliant
% GEOS-Chem Mechanism version 14.5.2-Hg. This file was generated
% by the GEOS-Chem Emulator on 2025-07-10 from the following KPP file:
%    /uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/GEOSChem_Emulator/Input_Files/GC_Mech_Files/v14.5.2-Hg/Hg.eqn 
% This file contains (only) the definitions of not-commented-out, named
% reactions rates that are referenced within:
%   /uufs/chpc.utah.edu/common/home/u6044586/MATLAB/F0AM_v4.4.2/Chem/GEOSChem/v14.5.2-Hg-reg/GEOSCHEM_HetRxns.m
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
%   required_functions={Het_HgIIP_Org(),SR_MW(),Het_HgIIP_Inorg()};


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%    Import Functions Used in Rate Constants    %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import the MATALB dictionary containing a cell array with the handles to
% individual rate constant functions. In this dictionary, the keys are the
% function names and the  values are the function handles (which are used
% to set rate constants).
rate_dict = import_GC_rates_v1452({'Het_HgIIP_Org','SR_MW','Het_HgIIP_Inorg'});

% Unpack function handles from the dictionary we just returned to use herein.
[Het_HgIIP_Org,SR_MW,Het_HgIIP_Inorg]= rate_dict{:};

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%    Initialize Rate Vars/Functions    %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the total number of rate constants
nk = 35;

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
Knames{i} = 'HgBrNO2_to_Hg2ORGP';
krx(:,i) = Het_HgIIP_Org(Met, State_Het, 0.1, SR_MW(ind_HgBrNO2));

i=i+1;
Knames{i} = 'HgBrHO2_to_Hg2ORGP';
krx(:,i) = Het_HgIIP_Org(Met, State_Het, 0.1, SR_MW(ind_HgBrHO2));

i=i+1;
Knames{i} = 'HgBrOH_to_Hg2ORGP';
krx(:,i) = Het_HgIIP_Org(Met, State_Het, 0.1, SR_MW(ind_HgBrOH ));

i=i+1;
Knames{i} = 'HgBrBrO_to_Hg2ORGP';
krx(:,i) = Het_HgIIP_Org(Met, State_Het, 0.1, SR_MW(ind_HgBrBrO));

i=i+1;
Knames{i} = 'HgBrClO_to_Hg2ORGP';
krx(:,i) = Het_HgIIP_Org(Met, State_Het, 0.1, SR_MW(ind_HgBrClO));

i=i+1;
Knames{i} = 'HgBr2_to_Hg2ORGP';
krx(:,i) = Het_HgIIP_Org(Met, State_Het, 0.1, SR_MW(ind_HgBr2  ));

i=i+1;
Knames{i} = 'HgClNO2_to_Hg2ORGP';
krx(:,i) = Het_HgIIP_Org(Met, State_Het, 0.1, SR_MW(ind_HgClNO2));

i=i+1;
Knames{i} = 'HgClHO2_to_Hg2ORGP';
krx(:,i) = Het_HgIIP_Org(Met, State_Het, 0.1, SR_MW(ind_HgClHO2));

i=i+1;
Knames{i} = 'HgClOH_to_Hg2ORGP';
krx(:,i) = Het_HgIIP_Org(Met, State_Het, 0.1, SR_MW(ind_HgClOH ));

i=i+1;
Knames{i} = 'HgClBrO_to_Hg2ORGP';
krx(:,i) = Het_HgIIP_Org(Met, State_Het, 0.1, SR_MW(ind_HgClBrO));

i=i+1;
Knames{i} = 'HgClClO_to_Hg2ORGP';
krx(:,i) = Het_HgIIP_Org(Met, State_Het, 0.1, SR_MW(ind_HgClClO));

i=i+1;
Knames{i} = 'HgClBr_to_Hg2ORGP';
krx(:,i) = Het_HgIIP_Org(Met, State_Het, 0.1, SR_MW(ind_HgClBr ));

i=i+1;
Knames{i} = 'HgCl2_to_Hg2ORGP';
krx(:,i) = Het_HgIIP_Org(Met, State_Het, 0.1, SR_MW(ind_HgCl2  ));

i=i+1;
Knames{i} = 'HgOHNO2_to_Hg2ORGP';
krx(:,i) = Het_HgIIP_Org(Met, State_Het, 0.1, SR_MW(ind_HgOHNO2));

i=i+1;
Knames{i} = 'HgOHHO2_to_Hg2ORGP';
krx(:,i) = Het_HgIIP_Org(Met, State_Het, 0.1, SR_MW(ind_HgOHHO2));

i=i+1;
Knames{i} = 'HgOHOH_to_Hg2ORGP';
krx(:,i) = Het_HgIIP_Org(Met, State_Het, 0.1, SR_MW(ind_HgOHOH ));

i=i+1;
Knames{i} = 'HgOHBrO_to_Hg2ORGP';
krx(:,i) = Het_HgIIP_Org(Met, State_Het, 0.1, SR_MW(ind_HgOHBrO));

i=i+1;
Knames{i} = 'HgOHClO_to_Hg2ORGP';
krx(:,i) = Het_HgIIP_Org(Met, State_Het, 0.1, SR_MW(ind_HgOHClO));

i=i+1;
Knames{i} = 'HgBrNO2_to_Hg2ClP';
krx(:,i) = Het_HgIIP_Inorg(Met, State_Het, 0.1, SR_MW(ind_HgBrNO2));

i=i+1;
Knames{i} = 'HgBrHO2_to_Hg2ClP';
krx(:,i) = Het_HgIIP_Inorg(Met, State_Het, 0.1, SR_MW(ind_HgBrHO2));

i=i+1;
Knames{i} = 'HgBrOH_to_Hg2ClP';
krx(:,i) = Het_HgIIP_Inorg(Met, State_Het, 0.1, SR_MW(ind_HgBrOH ));

i=i+1;
Knames{i} = 'HgBrBrO_to_Hg2ClP';
krx(:,i) = Het_HgIIP_Inorg(Met, State_Het, 0.1, SR_MW(ind_HgBrBrO));

i=i+1;
Knames{i} = 'HgBrClO_to_Hg2ClP';
krx(:,i) = Het_HgIIP_Inorg(Met, State_Het, 0.1, SR_MW(ind_HgBrClO));

i=i+1;
Knames{i} = 'HgBr2_to_Hg2ClP';
krx(:,i) = Het_HgIIP_Inorg(Met, State_Het, 0.1, SR_MW(ind_HgBr2  ));

i=i+1;
Knames{i} = 'HgClNO2_to_Hg2ClP';
krx(:,i) = Het_HgIIP_Inorg(Met, State_Het, 0.1, SR_MW(ind_HgClNO2));

i=i+1;
Knames{i} = 'HgClHO2_to_Hg2ClP';
krx(:,i) = Het_HgIIP_Inorg(Met, State_Het, 0.1, SR_MW(ind_HgClHO2));

i=i+1;
Knames{i} = 'HgClOH_to_Hg2ClP';
krx(:,i) = Het_HgIIP_Inorg(Met, State_Het, 0.1, SR_MW(ind_HgClOH ));

i=i+1;
Knames{i} = 'HgClBrO_to_Hg2ClP';
krx(:,i) = Het_HgIIP_Inorg(Met, State_Het, 0.1, SR_MW(ind_HgClBrO));

i=i+1;
Knames{i} = 'HgClClO_to_Hg2ClP';
krx(:,i) = Het_HgIIP_Inorg(Met, State_Het, 0.1, SR_MW(ind_HgClClO));

i=i+1;
Knames{i} = 'HgClBr_to_Hg2ClP';
krx(:,i) = Het_HgIIP_Inorg(Met, State_Het, 0.1, SR_MW(ind_HgClBr ));

i=i+1;
Knames{i} = 'HgOHNO2_to_Hg2ClP';
krx(:,i) = Het_HgIIP_Inorg(Met, State_Het, 0.1, SR_MW(ind_HgOHNO2));

i=i+1;
Knames{i} = 'HgOHHO2_to_Hg2ClP';
krx(:,i) = Het_HgIIP_Inorg(Met, State_Het, 0.1, SR_MW(ind_HgOHHO2));

i=i+1;
Knames{i} = 'HgOHOH_to_Hg2ClP';
krx(:,i) = Het_HgIIP_Inorg(Met, State_Het, 0.1, SR_MW(ind_HgOHOH ));

i=i+1;
Knames{i} = 'HgOHBrO_to_Hg2ClP';
krx(:,i) = Het_HgIIP_Inorg(Met, State_Het, 0.1, SR_MW(ind_HgOHBrO));

i=i+1;
Knames{i} = 'HgOHClO_to_Hg2ClP';
krx(:,i) = Het_HgIIP_Inorg(Met, State_Het, 0.1, SR_MW(ind_HgOHClO));


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