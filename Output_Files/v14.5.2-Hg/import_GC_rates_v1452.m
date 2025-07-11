function out =  import_GC_rates_v1452(cell_of_rates_to_get) 

    % Function to make these local functions available outside this file.
    % This "main function" returns function handles for each of the local functions.
    % that are defined in this file. So, this is just a hack to not have to save each of these
    % functions in a different file. It also lets me load all of them at once in a one liner
    % in the master mechanism file.
    %
    % Inputs: cell_of_rates_to_get = cell array of rates you'd like to
    %           retrieve or cell array equal to {'all'} if you want them all. 
    %
    % Outputs: cell array with values corresponding to a cell of len 1 containing
    %         function handles, returned in the same order as
    %         cell_of_rates_to_get, or alphabetically if passed "all". 
    %         See Examples to see how to unpack these appropariately in
    %         script. 
    % 
    % Change Log:
    %          JDH 8/16/2022 Created (github:{@jhaskinsPhD, jhaskins@alum.mit.edu)
    %          JDH 11/17/2022 Modified to use new MATLAB dictionary, for easier usage.
    %                         Requries MATLAB 2022b or higher b/c that's when
    %                         dictionaries were first implemented in MATALB. 
    %
    % Examples:
    % 
    %    (1) Either retrieve specific rates by passing their names as a cell in "cell_of_rates_to_get" as follows:
    %
    %    	out = import_GC_rates({'GCARR_ac','GC_OHHNO3_acacac'});
    %    	[GCARR_ac,GC_OHHNO3_acacac]=out{:};
    %    
    %    (2) Or use "all" in "cell_of_rates_to_get" to retrieve all defined functions in alphabetical order as follows:
    %
    %     	out = import_GC_rates({'all'});
    %     	[GCARR_ab,GCARR_abc,GCJPLPR_abab,GCARR_ac]= out{:}; 

    % Define all keys for dictionary 
    rate_keys=cell({'GCARR_ab','GCARR_abc','GCJPLPR_abab','GCARR_ac'}); 

    % Define all values for dictionary ... must be in SAME ORDER!!!!!!!!!!!
    rate_values=cell({@GCARR_ab,@GCARR_abc,@GCJPLPR_abab,@GCARR_ac});

    % Use keys and values to define a dictionary. 
    rate_dict=dictionary(rate_keys,rate_values);

    % If the user wants to import all the rates, then grab all keys in rate_dict.  
    if any(strcmp('all', cell_of_rates_to_get))
        cell_of_rates_to_get=rate_keys;
    end
    
    % Initialize all output arrays: 
    got_it=zeros([length(cell_of_rates_to_get),1],'double');
    baddie=repmat({''},[length(cell_of_rates_to_get),1]);
    out=repmat({''},[length(cell_of_rates_to_get),1]);
    
    % Loop over cell_of_rates_to_get and add values from dict to cell, "out"
    for i=1:length(cell_of_rates_to_get)
        rate_i=cell_of_rates_to_get(i); 
        if any(strcmp(rate_i, rate_keys))
            % If this rate is defined in this function... then grab it! 
            out(i)=rate_dict(rate_i); got_it(i)=1;
        else
            % If you don't have it, point to function "None"... prints error! 
            out(i)={@None}; got_it(i)=0; baddie(i)=rate_i; 
        end
    end
	
    % If nothing was ever added, or its all equal to @None, then print
    % error & pass back an empty cell, out. 
    if isempty(got_it(got_it==1))
        disp(['%%%%%%%%%  ERROR in import_GC_rates()  %%%%%%%%%%%%%%%%%%%%%%%%%%%%' ,newline, ...
            '     Could not find any valid names of functions in items passed. Returning empty!'])
        out=cell({}); 
    elseif ~isempty(got_it(got_it==0))
        bad='     '; bind=find(got_it==0); 
        for b =1:numel(bind); bad=strcat(bad,baddie{bind(b)},","); end  
        bad=char(bad); bad=bad(1:length(bad)-1);
        disp('%%%%%% WARNING %%%%%% in import_GC_rates(). Could not the following rates: ')
        disp(['     ',bad])
    end

end

function  [temp, press, numden, h2o,temp_over_k300, k300_over_temp,inv_temp,sr_temp,relhum] = extract_MetVars(Met) 
    % This function is used to extract invidividual met vars from the structure Met and is 
    % called at the top of all the rate functions defined in this file so they have access 
    % to these variables to use. These vars defined as "globals" in GEOS-Chem and don't appear 
    % as function input vars, so we just pass the struct "MET" as input to each function and 
    % use this one to then extract them out in each function ... 
	
    temp= Met.T;                % Temperature (K)
    press= Met.P;               % Pressure (mbar / hPa) (GC uses hPa)
    numden= Met.M;              % Number density (molecules/cm3)
    h2o= Met.H2O;               % Water vapor number density (molecules/cm3)
    temp_over_k300= Met.T./300; % Temp (K) divided by 300 K. 
    k300_over_temp= 300./Met.T;  % 300K divided by Temperature(K)
    inv_temp= 1./Met.T;         % Inverse Temperature (1/K) 
    sr_temp = sqrt(Met.T);      % Square root of Temperature (K^0.5)
    relhum = Met.RH;            % Relative Humididty as a percentage 
	
end

function [k] = None() 
    disp('Function pointed to "None", meaning you passed an arg to import_GC_rates() that it did not understand.')
    k=[];
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      BEGIN FORTRAN-->MATLAB Converted Functions from gcKPP_Rates.F90
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [k] = GCJPLPR_abab(Met,a1,b1,a2,b2,fv)
    % Pass input struct, "Met" to function "extract_MetVars()" to get vars this function may depend on...
    [temp, press, numden, h2o,temp_over_k300, k300_over_temp,inv_temp,sr_temp,relhum]= extract_MetVars(Met);

    % Third body effect for pressure dependence of rate coefficients.
    % a1, b1 are the Arrhenius parameters for the lower-limit rate.
    % a2, b2 are the Arrhenius parameters for the upper-limit rate.
    % fv     is the falloff curve paramter, (see ATKINSON ET. AL (1992)
    % J. Phys. Chem. Ref. Data 21, P. 1145). Usually fv = 0.6.
    % For these reactions, these Arrhenius law terms evaluate to 1:
    % exp(c1./T)
    % exp(c2./T)
    % because c1 = c2 = 0.  Therefore we can skip computing these
    % terms.  Also, fct1 = fct2 = 0, so we will skip computing these
    % terms as well.  This is more computationally efficient.
    % (bmy, 06 Jan 2022)
    % REAL(dp), INTENT(IN) :: a1,   b1,    a2,    b2,   fv
    % REAL(dp)             :: rlow, rhigh, xyrat, blog, fexp, k
    rlow  = a1.*( k300_over_temp.^b1 ).*numden;
    rhigh = a2.*( k300_over_temp.^b2 );
    xyrat = rlow./rhigh;
    blog  = log10( xyrat );
    fexp  = 1.0./( 1.0 + ( blog.*blog ) );
    k     = rlow.*( fv.^fexp )./( 1.0 + xyrat );
end
