% F0AM-compliant GEOS-Chem Mechanism for version 14.5.2-Hg generated 
% by the GEOS-Chem Emulator on 2025-07-10 from the following KPP file: 
%    /uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/GEOSChem_Emulator/Input_Files/GC_Mech_Files/v14.5.2-Hg/Hg.eqn
%
% Only heterogenous reactions in the GEOS-Chem mechanism are included in this
% file. Thus, this is a sub-selection of the mechanism in the KPP file above 
% intended to be used with the reactions in "GEOSCHEM_GasRxns.m" / not on its own.
%
% Thus, this mechanism includes: 
%     # of ADDITIONAL species  = 1        
%     # of ADDITIONAL reactions = 35  
% including ONLY those species/ reactions that are NOT declared in 
% "GEOSChem_GasRxns.m".
%
% The following named reaction rates are referenced in this file & must be 
% defined in "GEOSCHEM_K.m" (as named rxns) with the functions defined in 
% the "import_GC_rates.m" file. 
%    Het_HgIIP_Org(),SR_MW(),Het_HgIIP_Inorg()

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%    BEGIN F0AM MECHANISM     %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SpeciesToAdd = {'Hg2ClP'};

RO2ToAdd = {}; % In a GEOS-Chem Mechanism, this list should ALWAYS be empty. 

AddSpecies

i=i+1;
Rnames{i}='HgBrNO2=Hg2ORGP';
k(:,i) = HgBrNO2_to_Hg2ORGP;
Gstr{i,1}='HgBrNO2'; 
fHgBrNO2(i)=fHgBrNO2(i)-1.0; fHg2ORGP(i)=fHg2ORGP(i)+1.0;

i=i+1;
Rnames{i}='HgBrHO2=Hg2ORGP';
k(:,i) = HgBrHO2_to_Hg2ORGP;
Gstr{i,1}='HgBrHO2'; 
fHgBrHO2(i)=fHgBrHO2(i)-1.0; fHg2ORGP(i)=fHg2ORGP(i)+1.0;

i=i+1;
Rnames{i}='HgBrOH=Hg2ORGP';
k(:,i) = HgBrOH_to_Hg2ORGP;
Gstr{i,1}='HgBrOH'; 
fHgBrOH(i)=fHgBrOH(i)-1.0; fHg2ORGP(i)=fHg2ORGP(i)+1.0;

i=i+1;
Rnames{i}='HgBrBrO=Hg2ORGP';
k(:,i) = HgBrBrO_to_Hg2ORGP;
Gstr{i,1}='HgBrBrO'; 
fHgBrBrO(i)=fHgBrBrO(i)-1.0; fHg2ORGP(i)=fHg2ORGP(i)+1.0;

i=i+1;
Rnames{i}='HgBrClO=Hg2ORGP';
k(:,i) = HgBrClO_to_Hg2ORGP;
Gstr{i,1}='HgBrClO'; 
fHgBrClO(i)=fHgBrClO(i)-1.0; fHg2ORGP(i)=fHg2ORGP(i)+1.0;

i=i+1;
Rnames{i}='HgBr2=Hg2ORGP';
k(:,i) = HgBr2_to_Hg2ORGP;
Gstr{i,1}='HgBr2'; 
fHgBr2(i)=fHgBr2(i)-1.0; fHg2ORGP(i)=fHg2ORGP(i)+1.0;

i=i+1;
Rnames{i}='HgClNO2=Hg2ORGP';
k(:,i) = HgClNO2_to_Hg2ORGP;
Gstr{i,1}='HgClNO2'; 
fHgClNO2(i)=fHgClNO2(i)-1.0; fHg2ORGP(i)=fHg2ORGP(i)+1.0;

i=i+1;
Rnames{i}='HgClHO2=Hg2ORGP';
k(:,i) = HgClHO2_to_Hg2ORGP;
Gstr{i,1}='HgClHO2'; 
fHgClHO2(i)=fHgClHO2(i)-1.0; fHg2ORGP(i)=fHg2ORGP(i)+1.0;

i=i+1;
Rnames{i}='HgClOH=Hg2ORGP';
k(:,i) = HgClOH_to_Hg2ORGP;
Gstr{i,1}='HgClOH'; 
fHgClOH(i)=fHgClOH(i)-1.0; fHg2ORGP(i)=fHg2ORGP(i)+1.0;

i=i+1;
Rnames{i}='HgClBrO=Hg2ORGP';
k(:,i) = HgClBrO_to_Hg2ORGP;
Gstr{i,1}='HgClBrO'; 
fHgClBrO(i)=fHgClBrO(i)-1.0; fHg2ORGP(i)=fHg2ORGP(i)+1.0;

i=i+1;
Rnames{i}='HgClClO=Hg2ORGP';
k(:,i) = HgClClO_to_Hg2ORGP;
Gstr{i,1}='HgClClO'; 
fHgClClO(i)=fHgClClO(i)-1.0; fHg2ORGP(i)=fHg2ORGP(i)+1.0;

i=i+1;
Rnames{i}='HgClBr=Hg2ORGP';
k(:,i) = HgClBr_to_Hg2ORGP;
Gstr{i,1}='HgClBr'; 
fHgClBr(i)=fHgClBr(i)-1.0; fHg2ORGP(i)=fHg2ORGP(i)+1.0;

i=i+1;
Rnames{i}='HgCl2=Hg2ORGP';
k(:,i) = HgCl2_to_Hg2ORGP;
Gstr{i,1}='HgCl2'; 
fHgCl2(i)=fHgCl2(i)-1.0; fHg2ORGP(i)=fHg2ORGP(i)+1.0;

i=i+1;
Rnames{i}='HgOHNO2=Hg2ORGP';
k(:,i) = HgOHNO2_to_Hg2ORGP;
Gstr{i,1}='HgOHNO2'; 
fHgOHNO2(i)=fHgOHNO2(i)-1.0; fHg2ORGP(i)=fHg2ORGP(i)+1.0;

i=i+1;
Rnames{i}='HgOHHO2=Hg2ORGP';
k(:,i) = HgOHHO2_to_Hg2ORGP;
Gstr{i,1}='HgOHHO2'; 
fHgOHHO2(i)=fHgOHHO2(i)-1.0; fHg2ORGP(i)=fHg2ORGP(i)+1.0;

i=i+1;
Rnames{i}='HgOHOH=Hg2ORGP';
k(:,i) = HgOHOH_to_Hg2ORGP;
Gstr{i,1}='HgOHOH'; 
fHgOHOH(i)=fHgOHOH(i)-1.0; fHg2ORGP(i)=fHg2ORGP(i)+1.0;

i=i+1;
Rnames{i}='HgOHBrO=Hg2ORGP';
k(:,i) = HgOHBrO_to_Hg2ORGP;
Gstr{i,1}='HgOHBrO'; 
fHgOHBrO(i)=fHgOHBrO(i)-1.0; fHg2ORGP(i)=fHg2ORGP(i)+1.0;

i=i+1;
Rnames{i}='HgOHClO=Hg2ORGP';
k(:,i) = HgOHClO_to_Hg2ORGP;
Gstr{i,1}='HgOHClO'; 
fHgOHClO(i)=fHgOHClO(i)-1.0; fHg2ORGP(i)=fHg2ORGP(i)+1.0;

i=i+1;
Rnames{i}='HgBrNO2=Hg2ClP';
k(:,i) = HgBrNO2_to_Hg2ClP;
Gstr{i,1}='HgBrNO2'; 
fHgBrNO2(i)=fHgBrNO2(i)-1.0; fHg2ClP(i)=fHg2ClP(i)+1.0;

i=i+1;
Rnames{i}='HgBrHO2=Hg2ClP';
k(:,i) = HgBrHO2_to_Hg2ClP;
Gstr{i,1}='HgBrHO2'; 
fHgBrHO2(i)=fHgBrHO2(i)-1.0; fHg2ClP(i)=fHg2ClP(i)+1.0;

i=i+1;
Rnames{i}='HgBrOH=Hg2ClP';
k(:,i) = HgBrOH_to_Hg2ClP;
Gstr{i,1}='HgBrOH'; 
fHgBrOH(i)=fHgBrOH(i)-1.0; fHg2ClP(i)=fHg2ClP(i)+1.0;

i=i+1;
Rnames{i}='HgBrBrO=Hg2ClP';
k(:,i) = HgBrBrO_to_Hg2ClP;
Gstr{i,1}='HgBrBrO'; 
fHgBrBrO(i)=fHgBrBrO(i)-1.0; fHg2ClP(i)=fHg2ClP(i)+1.0;

i=i+1;
Rnames{i}='HgBrClO=Hg2ClP';
k(:,i) = HgBrClO_to_Hg2ClP;
Gstr{i,1}='HgBrClO'; 
fHgBrClO(i)=fHgBrClO(i)-1.0; fHg2ClP(i)=fHg2ClP(i)+1.0;

i=i+1;
Rnames{i}='HgBr2=Hg2ClP';
k(:,i) = HgBr2_to_Hg2ClP;
Gstr{i,1}='HgBr2'; 
fHgBr2(i)=fHgBr2(i)-1.0; fHg2ClP(i)=fHg2ClP(i)+1.0;

i=i+1;
Rnames{i}='HgClNO2=Hg2ClP';
k(:,i) = HgClNO2_to_Hg2ClP;
Gstr{i,1}='HgClNO2'; 
fHgClNO2(i)=fHgClNO2(i)-1.0; fHg2ClP(i)=fHg2ClP(i)+1.0;

i=i+1;
Rnames{i}='HgClHO2=Hg2ClP';
k(:,i) = HgClHO2_to_Hg2ClP;
Gstr{i,1}='HgClHO2'; 
fHgClHO2(i)=fHgClHO2(i)-1.0; fHg2ClP(i)=fHg2ClP(i)+1.0;

i=i+1;
Rnames{i}='HgClOH=Hg2ClP';
k(:,i) = HgClOH_to_Hg2ClP;
Gstr{i,1}='HgClOH'; 
fHgClOH(i)=fHgClOH(i)-1.0; fHg2ClP(i)=fHg2ClP(i)+1.0;

i=i+1;
Rnames{i}='HgClBrO=Hg2ClP';
k(:,i) = HgClBrO_to_Hg2ClP;
Gstr{i,1}='HgClBrO'; 
fHgClBrO(i)=fHgClBrO(i)-1.0; fHg2ClP(i)=fHg2ClP(i)+1.0;

i=i+1;
Rnames{i}='HgClClO=Hg2ClP';
k(:,i) = HgClClO_to_Hg2ClP;
Gstr{i,1}='HgClClO'; 
fHgClClO(i)=fHgClClO(i)-1.0; fHg2ClP(i)=fHg2ClP(i)+1.0;

i=i+1;
Rnames{i}='HgClBr=Hg2ClP';
k(:,i) = HgClBr_to_Hg2ClP;
Gstr{i,1}='HgClBr'; 
fHgClBr(i)=fHgClBr(i)-1.0; fHg2ClP(i)=fHg2ClP(i)+1.0;

i=i+1;
Rnames{i}='HgOHNO2=Hg2ClP';
k(:,i) = HgOHNO2_to_Hg2ClP;
Gstr{i,1}='HgOHNO2'; 
fHgOHNO2(i)=fHgOHNO2(i)-1.0; fHg2ClP(i)=fHg2ClP(i)+1.0;

i=i+1;
Rnames{i}='HgOHHO2=Hg2ClP';
k(:,i) = HgOHHO2_to_Hg2ClP;
Gstr{i,1}='HgOHHO2'; 
fHgOHHO2(i)=fHgOHHO2(i)-1.0; fHg2ClP(i)=fHg2ClP(i)+1.0;

i=i+1;
Rnames{i}='HgOHOH=Hg2ClP';
k(:,i) = HgOHOH_to_Hg2ClP;
Gstr{i,1}='HgOHOH'; 
fHgOHOH(i)=fHgOHOH(i)-1.0; fHg2ClP(i)=fHg2ClP(i)+1.0;

i=i+1;
Rnames{i}='HgOHBrO=Hg2ClP';
k(:,i) = HgOHBrO_to_Hg2ClP;
Gstr{i,1}='HgOHBrO'; 
fHgOHBrO(i)=fHgOHBrO(i)-1.0; fHg2ClP(i)=fHg2ClP(i)+1.0;

i=i+1;
Rnames{i}='HgOHClO=Hg2ClP';
k(:,i) = HgOHClO_to_Hg2ClP;
Gstr{i,1}='HgOHClO'; 
fHgOHClO(i)=fHgOHClO(i)-1.0; fHg2ClP(i)=fHg2ClP(i)+1.0;

