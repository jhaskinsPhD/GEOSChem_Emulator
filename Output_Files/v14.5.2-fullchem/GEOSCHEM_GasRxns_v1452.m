% F0AM-compliant GEOS-Chem Mechanism for version 14.5.2-Hg generated
% by the GEOS-Chem Emulator on 2025-07-10 from the following KPP file(s):
%    /uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/GEOSChem_Emulator/Input_Files/GC_Mech_Files/v14.5.2-Hg/Hg.eqn
%
%    Photolysis Frequency Mapping JMAP_TYPE = CrossSec_Match
%    # of (non-het) species                 = 37
%    # of (non-het) reactions               = 60
% 
% The following (named) photolysis rates must be defined in "GEOSChem_J.m"
% for this mechanism to work:
%    required_Js = {'jNO2_b','jClO','jBrO','jMB','jMBNO2','jMBHO2',...
%    'jMBOH','jMB2','jMBBRO','jMBCLO','jMC2','jMO','jNO2_a'};
%
% The following rate functions used in the named reaction rates used here/
% called within "GEOSCHEM_K.m" must be defined in the "import_GC_rates.m"
% file for this mechanism to work: 
%    required_rates = {GCARR_ab(),GCARR_abc(),GCJPLPR_abab(),GCARR_ac()};

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%   BEGIN F0AM  MECHANISM     %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SpeciesToAdd = {'Br';'BrO';'CH4';'CO';'Cl';'ClO';'HO2';'Hg0';...
'Hg2ORGP';'HgBr';'HgBr2';'HgBrBrO';'HgBrClO';'HgBrHO2';'HgBrNO2';'HgBrO';'HgBrOH';...
'HgCl';'HgCl2';'HgClBr';'HgClBrO';'HgClClO';'HgClHO2';'HgClNO2';'HgClO';'HgClOH';...
'HgOH';'HgOHBrO';'HgOHClO';'HgOHHO2';'HgOHNO2';'HgOHO';'HgOHOH';'NO';'NO2';'O3';...
'OH'};

RO2ToAdd = {}; % In a GEOS-Chem Mechanism, this list should ALWAYS be empty.

AddSpecies


i=i+1;
Rnames{i}='Hg0+Br=HgBr';
k(:,i) = Hg0_Br_to_HgBr;
Gstr{i,1}='Hg0';  Gstr{i,2}='Br'; 
fHg0(i)=fHg0(i)-1.0; fBr(i)=fBr(i)-1.0; fHgBr(i)=fHgBr(i)+1.0;

i=i+1;
Rnames{i}='HgBr=Hg0';
k(:,i) = HgBr_to_Hg0;
Gstr{i,1}='HgBr'; 
fHgBr(i)=fHgBr(i)-1.0; fHg0(i)=fHg0(i)+1.0;

i=i+1;
Rnames{i}='HgBr+Br=Hg0';
k(:,i) = 3.9e-11;
Gstr{i,1}='HgBr';  Gstr{i,2}='Br'; 
fHgBr(i)=fHgBr(i)-1.0; fBr(i)=fBr(i)-1.0; fHg0(i)=fHg0(i)+1.0;

i=i+1;
Rnames{i}='HgBr+NO2=Hg0';
k(:,i) = 3.0e-12;
Gstr{i,1}='HgBr';  Gstr{i,2}='NO2'; 
fHgBr(i)=fHgBr(i)-1.0; fNO2(i)=fNO2(i)-1.0; fHg0(i)=fHg0(i)+1.0;

i=i+1;
Rnames{i}='HgBr+NO2=HgBrNO2';
k(:,i) = HgBr_NO2_to_HgBrNO2;
Gstr{i,1}='HgBr';  Gstr{i,2}='NO2'; 
fHgBr(i)=fHgBr(i)-1.0; fNO2(i)=fNO2(i)-1.0; fHgBrNO2(i)=fHgBrNO2(i)+1.0;

i=i+1;
Rnames{i}='HgBr+HO2=HgBrHO2';
k(:,i) = HgBr_HO2_to_HgBrHO2;
Gstr{i,1}='HgBr';  Gstr{i,2}='HO2'; 
fHgBr(i)=fHgBr(i)-1.0; fHO2(i)=fHO2(i)-1.0; fHgBrHO2(i)=fHgBrHO2(i)+1.0;

i=i+1;
Rnames{i}='HgBr+ClO=HgBrClO';
k(:,i) = HgBr_ClO_to_HgBrClO;
Gstr{i,1}='HgBr';  Gstr{i,2}='ClO'; 
fHgBr(i)=fHgBr(i)-1.0; fClO(i)=fClO(i)-1.0; fHgBrClO(i)=fHgBrClO(i)+1.0;

i=i+1;
Rnames{i}='HgBr+BrO=HgBrBrO';
k(:,i) = HgBr_BrO_to_HgBrBrO;
Gstr{i,1}='HgBr';  Gstr{i,2}='BrO'; 
fHgBr(i)=fHgBr(i)-1.0; fBrO(i)=fBrO(i)-1.0; fHgBrBrO(i)=fHgBrBrO(i)+1.0;

i=i+1;
Rnames{i}='HgBr+OH=HgBrOH';
k(:,i) = HgBr_OH_to_HgBrOH;
Gstr{i,1}='HgBr';  Gstr{i,2}='OH'; 
fHgBr(i)=fHgBr(i)-1.0; fOH(i)=fOH(i)-1.0; fHgBrOH(i)=fHgBrOH(i)+1.0;

i=i+1;
Rnames{i}='HgBr+Br=HgBr2';
k(:,i) = 3.0e-11;
Gstr{i,1}='HgBr';  Gstr{i,2}='Br'; 
fHgBr(i)=fHgBr(i)-1.0; fBr(i)=fBr(i)-1.0; fHgBr2(i)=fHgBr2(i)+1.0;

i=i+1;
Rnames{i}='HgBrO+CH4=HgBrOH';
k(:,i) = HgBrO_CH4_to_HgBrOH;
Gstr{i,1}='HgBrO';  Gstr{i,2}='CH4'; 
fHgBrO(i)=fHgBrO(i)-1.0; fCH4(i)=fCH4(i)-1.0; fHgBrOH(i)=fHgBrOH(i)+1.0;

i=i+1;
Rnames{i}='HgBrO+CO=HgBr';
k(:,i) = HgBrO_CO_to_HgBr;
Gstr{i,1}='HgBrO';  Gstr{i,2}='CO'; 
fHgBrO(i)=fHgBrO(i)-1.0; fCO(i)=fCO(i)-1.0; fHgBr(i)=fHgBr(i)+1.0;

i=i+1;
Rnames{i}='HgBr+O3=HgBrO';
k(:,i) = 3.0e-11;
Gstr{i,1}='HgBr';  Gstr{i,2}='O3'; 
fHgBr(i)=fHgBr(i)-1.0; fO3(i)=fO3(i)-1.0; fHgBrO(i)=fHgBrO(i)+1.0;

i=i+1;
Rnames{i}='Hg0+Cl=HgCl';
k(:,i) = Hg0_Cl_to_HgCl;
Gstr{i,1}='Hg0';  Gstr{i,2}='Cl'; 
fHg0(i)=fHg0(i)-1.0; fCl(i)=fCl(i)-1.0; fHgCl(i)=fHgCl(i)+1.0;

i=i+1;
Rnames{i}='HgCl+Cl=Hg0';
k(:,i) = 3.9e-11;
Gstr{i,1}='HgCl';  Gstr{i,2}='Cl'; 
fHgCl(i)=fHgCl(i)-1.0; fCl(i)=fCl(i)-1.0; fHg0(i)=fHg0(i)+1.0;

i=i+1;
Rnames{i}='HgCl+NO2=Hg0';
k(:,i) = 3.0e-12;
Gstr{i,1}='HgCl';  Gstr{i,2}='NO2'; 
fHgCl(i)=fHgCl(i)-1.0; fNO2(i)=fNO2(i)-1.0; fHg0(i)=fHg0(i)+1.0;

i=i+1;
Rnames{i}='HgCl+NO2=HgClNO2';
k(:,i) = HgCl_NO2_to_HgClNO2;
Gstr{i,1}='HgCl';  Gstr{i,2}='NO2'; 
fHgCl(i)=fHgCl(i)-1.0; fNO2(i)=fNO2(i)-1.0; fHgClNO2(i)=fHgClNO2(i)+1.0;

i=i+1;
Rnames{i}='HgCl+HO2=HgClHO2';
k(:,i) = HgCl_HO2_to_HgClHO2;
Gstr{i,1}='HgCl';  Gstr{i,2}='HO2'; 
fHgCl(i)=fHgCl(i)-1.0; fHO2(i)=fHO2(i)-1.0; fHgClHO2(i)=fHgClHO2(i)+1.0;

i=i+1;
Rnames{i}='HgCl+ClO=HgClClO';
k(:,i) = HgCl_ClO_to_HgClClO;
Gstr{i,1}='HgCl';  Gstr{i,2}='ClO'; 
fHgCl(i)=fHgCl(i)-1.0; fClO(i)=fClO(i)-1.0; fHgClClO(i)=fHgClClO(i)+1.0;

i=i+1;
Rnames{i}='HgCl+BrO=HgClBrO';
k(:,i) = HgCl_BrO_to_HgClBrO;
Gstr{i,1}='HgCl';  Gstr{i,2}='BrO'; 
fHgCl(i)=fHgCl(i)-1.0; fBrO(i)=fBrO(i)-1.0; fHgClBrO(i)=fHgClBrO(i)+1.0;

i=i+1;
Rnames{i}='HgCl+Br=HgClBr';
k(:,i) = 3.0e-11;
Gstr{i,1}='HgCl';  Gstr{i,2}='Br'; 
fHgCl(i)=fHgCl(i)-1.0; fBr(i)=fBr(i)-1.0; fHgClBr(i)=fHgClBr(i)+1.0;

i=i+1;
Rnames{i}='HgCl+OH=HgClOH';
k(:,i) = 3.0e-11;
Gstr{i,1}='HgCl';  Gstr{i,2}='OH'; 
fHgCl(i)=fHgCl(i)-1.0; fOH(i)=fOH(i)-1.0; fHgClOH(i)=fHgClOH(i)+1.0;

i=i+1;
Rnames{i}='HgCl+O3=HgClO';
k(:,i) = 3.0e-11;
Gstr{i,1}='HgCl';  Gstr{i,2}='O3'; 
fHgCl(i)=fHgCl(i)-1.0; fO3(i)=fO3(i)-1.0; fHgClO(i)=fHgClO(i)+1.0;

i=i+1;
Rnames{i}='HgClO+CH4=HgClOH';
k(:,i) = HgClO_CH4_to_HgClOH;
Gstr{i,1}='HgClO';  Gstr{i,2}='CH4'; 
fHgClO(i)=fHgClO(i)-1.0; fCH4(i)=fCH4(i)-1.0; fHgClOH(i)=fHgClOH(i)+1.0;

i=i+1;
Rnames{i}='HgClO+CO=HgCl';
k(:,i) = HgClO_CO_to_HgCl;
Gstr{i,1}='HgClO';  Gstr{i,2}='CO'; 
fHgClO(i)=fHgClO(i)-1.0; fCO(i)=fCO(i)-1.0; fHgCl(i)=fHgCl(i)+1.0;

i=i+1;
Rnames{i}='Hg0+OH=HgOH';
k(:,i) = Hg0_OH_to_HgOH;
Gstr{i,1}='Hg0';  Gstr{i,2}='OH'; 
fHg0(i)=fHg0(i)-1.0; fOH(i)=fOH(i)-1.0; fHgOH(i)=fHgOH(i)+1.0;

i=i+1;
Rnames{i}='HgOH=Hg0';
k(:,i) = HgOH_to_Hg0;
Gstr{i,1}='HgOH'; 
fHgOH(i)=fHgOH(i)-1.0; fHg0(i)=fHg0(i)+1.0;

i=i+1;
Rnames{i}='HgOH+NO2=HgOHNO2';
k(:,i) = HgOH_NO2_to_HgOHNO2;
Gstr{i,1}='HgOH';  Gstr{i,2}='NO2'; 
fHgOH(i)=fHgOH(i)-1.0; fNO2(i)=fNO2(i)-1.0; fHgOHNO2(i)=fHgOHNO2(i)+1.0;

i=i+1;
Rnames{i}='HgOH+HO2=HgOHHO2';
k(:,i) = HgOH_HO2_to_HgOHHO2;
Gstr{i,1}='HgOH';  Gstr{i,2}='HO2'; 
fHgOH(i)=fHgOH(i)-1.0; fHO2(i)=fHO2(i)-1.0; fHgOHHO2(i)=fHgOHHO2(i)+1.0;

i=i+1;
Rnames{i}='HgOH+ClO=HgOHClO';
k(:,i) = HgOH_ClO_to_HgOHClO;
Gstr{i,1}='HgOH';  Gstr{i,2}='ClO'; 
fHgOH(i)=fHgOH(i)-1.0; fClO(i)=fClO(i)-1.0; fHgOHClO(i)=fHgOHClO(i)+1.0;

i=i+1;
Rnames{i}='HgOH+BrO=HgOHBrO';
k(:,i) = HgOH_BrO_to_HgOHBrO;
Gstr{i,1}='HgOH';  Gstr{i,2}='BrO'; 
fHgOH(i)=fHgOH(i)-1.0; fBrO(i)=fBrO(i)-1.0; fHgOHBrO(i)=fHgOHBrO(i)+1.0;

i=i+1;
Rnames{i}='HgOH+Br=HgBrOH';
k(:,i) = HgOH_Br_to_HgBrOH;
Gstr{i,1}='HgOH';  Gstr{i,2}='Br'; 
fHgOH(i)=fHgOH(i)-1.0; fBr(i)=fBr(i)-1.0; fHgBrOH(i)=fHgBrOH(i)+1.0;

i=i+1;
Rnames{i}='HgOH+OH=HgOHOH';
k(:,i) = HgOH_OH_to_HgOHOH;
Gstr{i,1}='HgOH';  Gstr{i,2}='OH'; 
fHgOH(i)=fHgOH(i)-1.0; fOH(i)=fOH(i)-1.0; fHgOHOH(i)=fHgOHOH(i)+1.0;

i=i+1;
Rnames{i}='HgOH+O3=HgOHO';
k(:,i) = 3.0e-11;
Gstr{i,1}='HgOH';  Gstr{i,2}='O3'; 
fHgOH(i)=fHgOH(i)-1.0; fO3(i)=fO3(i)-1.0; fHgOHO(i)=fHgOHO(i)+1.0;

i=i+1;
Rnames{i}='HgOHO+CH4=HgOHOH';
k(:,i) = HgOHO_CH4_to_HgOHOH;
Gstr{i,1}='HgOHO';  Gstr{i,2}='CH4'; 
fHgOHO(i)=fHgOHO(i)-1.0; fCH4(i)=fCH4(i)-1.0; fHgOHOH(i)=fHgOHOH(i)+1.0;

i=i+1;
Rnames{i}='HgOHO+CO=HgOH';
k(:,i) = HgOHO_CO_to_HgOH;
Gstr{i,1}='HgOHO';  Gstr{i,2}='CO'; 
fHgOHO(i)=fHgOHO(i)-1.0; fCO(i)=fCO(i)-1.0; fHgOH(i)=fHgOH(i)+1.0;

i=i+1;
Rnames{i}='NO2+hv=NO+O3';
k(:,i)=jNO2_b; % PHOTOL(25) = jNO2_b & used for NO2+hv->NO+O
Gstr{i,1}='NO2'; 
fNO2(i)=fNO2(i)-1.0; fNO(i)=fNO(i)+1.0; fO3(i)=fO3(i)+1.0;

i=i+1;
Rnames{i}='BrO+hv=Br+O3';
k(:,i)=jClO; % PHOTOL(26) = jClO & used for ClO+hv->Cl+O
Gstr{i,1}='BrO'; 
fBrO(i)=fBrO(i)-1.0; fBr(i)=fBr(i)+1.0; fO3(i)=fO3(i)+1.0;

i=i+1;
Rnames{i}='ClO+hv=Cl+O3';
k(:,i)=jBrO; % PHOTOL(27) = jBrO & used for BrO+hv->Br+O
Gstr{i,1}='ClO'; 
fClO(i)=fClO(i)-1.0; fCl(i)=fCl(i)+1.0; fO3(i)=fO3(i)+1.0;

i=i+1;
Rnames{i}='HgBr+hv=Hg0+Br';
k(:,i)=jMB; % PHOTOL(4) = jMB & used for HGBR+hv->HG0+BR
Gstr{i,1}='HgBr'; 
fHgBr(i)=fHgBr(i)-1.0; fHg0(i)=fHg0(i)+1.0; fBr(i)=fBr(i)+1.0;

i=i+1;
Rnames{i}='HgBrNO2+hv=0.9HgBrO+0.1HgBr+0.9NO+0.1NO2';
k(:,i)=jMBNO2; % PHOTOL(20) = jMBNO2 & used for HGOHNO2+hv->HGOH+PRODUCTS
Gstr{i,1}='HgBrNO2'; 
fHgBrNO2(i)=fHgBrNO2(i)-1.0; fHgBrO(i)=fHgBrO(i)+0.9; fHgBr(i)=fHgBr(i)+0.1; fNO(i)=fNO(i)+0.9; fNO2(i)=fNO2(i)+0.1;

i=i+1;
Rnames{i}='HgBrHO2+hv=0.25HgBrO+0.67Hg0+0.07HgBr+0.01HgBrOH+0.67Br+0.74HO2+0.26OH';
k(:,i)=jMBHO2; % PHOTOL(21) = jMBHO2 & used for HGOHHO2+hv->HGOH+PRODUCTS
Gstr{i,1}='HgBrHO2'; 
fHgBrHO2(i)=fHgBrHO2(i)-1.0; fHgBrO(i)=fHgBrO(i)+0.25; fHg0(i)=fHg0(i)+0.67; fHgBr(i)=fHgBr(i)+0.07; fHgBrOH(i)=fHgBrOH(i)+0.01; fBr(i)=fBr(i)+0.67; fHO2(i)=fHO2(i)+0.74; fOH(i)=fOH(i)+0.26;

i=i+1;
Rnames{i}='HgBrOH+hv=0.49Hg0+0.35HgOH+0.15HgBr+0.01HgBrOH+0.85Br+0.65OH';
k(:,i)=jMBOH; % PHOTOL(16) = jMBOH & used for HGCLOH+hv->HG0+PRODUCTS
Gstr{i,1}='HgBrOH'; 
fHgBrOH(i)=fHgBrOH(i)-1.0; fHg0(i)=fHg0(i)+0.49; fHgOH(i)=fHgOH(i)+0.35; fHgBr(i)=fHgBr(i)+0.15; fHgBrOH(i)=fHgBrOH(i)+0.01; fBr(i)=fBr(i)+0.85; fOH(i)=fOH(i)+0.65;

i=i+1;
Rnames{i}='HgBr2+hv=0.4Hg0+0.6HgBr+1.4Br';
k(:,i)=jMB2; % PHOTOL(17) = jMB2 & used for HGCLBR+hv->HG0+PRODUCTS
Gstr{i,1}='HgBr2'; 
fHgBr2(i)=fHgBr2(i)-1.0; fHg0(i)=fHg0(i)+0.4; fHgBr(i)=fHgBr(i)+0.6; fBr(i)=fBr(i)+1.4;

i=i+1;
Rnames{i}='HgBrBrO+hv=Hg0+BrO';
k(:,i)=jMBBRO; % PHOTOL(23) = jMBBRO & used for HGOHBRO+hv->HGOH+PRODUCTS
Gstr{i,1}='HgBrBrO'; 
fHgBrBrO(i)=fHgBrBrO(i)-1.0; fHg0(i)=fHg0(i)+1.0; fBrO(i)=fBrO(i)+1.0;

i=i+1;
Rnames{i}='HgBrClO+hv=Hg0+Cl+Br';
k(:,i)=jMBCLO; % PHOTOL(24) = jMBCLO & used for HGOHCLO+hv->HGOH+PRODUCTS
Gstr{i,1}='HgBrClO'; 
fHgBrClO(i)=fHgBrClO(i)-1.0; fHg0(i)=fHg0(i)+1.0; fCl(i)=fCl(i)+1.0; fBr(i)=fBr(i)+1.0;

i=i+1;
Rnames{i}='HgCl2+hv=Hg0+2Cl';
k(:,i)=jMC2; % PHOTOL(22) = jMC2 & used for HGOHOH+hv->HGOH+PRODUCTS
Gstr{i,1}='HgCl2'; 
fHgCl2(i)=fHgCl2(i)-1.0; fHg0(i)=fHg0(i)+1.0; fCl(i)=fCl(i)+2.0;

i=i+1;
Rnames{i}='HgClNO2+hv=0.9HgClO+0.1HgCl+0.9NO+0.1NO2';
k(:,i)=jMBNO2; % PHOTOL(20) = jMBNO2 & used for HGOHNO2+hv->HGOH+PRODUCTS
Gstr{i,1}='HgClNO2'; 
fHgClNO2(i)=fHgClNO2(i)-1.0; fHgClO(i)=fHgClO(i)+0.9; fHgCl(i)=fHgCl(i)+0.1; fNO(i)=fNO(i)+0.9; fNO2(i)=fNO2(i)+0.1;

i=i+1;
Rnames{i}='HgClHO2+hv=0.25HgClO+0.67Hg0+0.07HgCl+0.01HgClOH+0.67Cl+0.74HO2+0.26OH';
k(:,i)=jMBHO2; % PHOTOL(21) = jMBHO2 & used for HGOHHO2+hv->HGOH+PRODUCTS
Gstr{i,1}='HgClHO2'; 
fHgClHO2(i)=fHgClHO2(i)-1.0; fHgClO(i)=fHgClO(i)+0.25; fHg0(i)=fHg0(i)+0.67; fHgCl(i)=fHgCl(i)+0.07; fHgClOH(i)=fHgClOH(i)+0.01; fCl(i)=fCl(i)+0.67; fHO2(i)=fHO2(i)+0.74; fOH(i)=fOH(i)+0.26;

i=i+1;
Rnames{i}='HgClOH+hv=0.49Hg0+0.35HgOH+0.15HgCl+0.01HgClOH+0.85Cl+0.65OH';
k(:,i)=jMBOH; % PHOTOL(16) = jMBOH & used for HGCLOH+hv->HG0+PRODUCTS
Gstr{i,1}='HgClOH'; 
fHgClOH(i)=fHgClOH(i)-1.0; fHg0(i)=fHg0(i)+0.49; fHgOH(i)=fHgOH(i)+0.35; fHgCl(i)=fHgCl(i)+0.15; fHgClOH(i)=fHgClOH(i)+0.01; fCl(i)=fCl(i)+0.85; fOH(i)=fOH(i)+0.65;

i=i+1;
Rnames{i}='HgClBr+hv=HgCl+Br';
k(:,i)=jMB2; % PHOTOL(17) = jMB2 & used for HGCLBR+hv->HG0+PRODUCTS
Gstr{i,1}='HgClBr'; 
fHgClBr(i)=fHgClBr(i)-1.0; fHgCl(i)=fHgCl(i)+1.0; fBr(i)=fBr(i)+1.0;

i=i+1;
Rnames{i}='HgClBrO+hv=HgCl+BrO';
k(:,i)=jMBBRO; % PHOTOL(23) = jMBBRO & used for HGOHBRO+hv->HGOH+PRODUCTS
Gstr{i,1}='HgClBrO'; 
fHgClBrO(i)=fHgClBrO(i)-1.0; fHgCl(i)=fHgCl(i)+1.0; fBrO(i)=fBrO(i)+1.0;

i=i+1;
Rnames{i}='HgClClO+hv=HgCl+ClO';
k(:,i)=jMBCLO; % PHOTOL(24) = jMBCLO & used for HGOHCLO+hv->HGOH+PRODUCTS
Gstr{i,1}='HgClClO'; 
fHgClClO(i)=fHgClClO(i)-1.0; fHgCl(i)=fHgCl(i)+1.0; fClO(i)=fClO(i)+1.0;

i=i+1;
Rnames{i}='HgOH+hv=Hg0+OH';
k(:,i)=jMO; % PHOTOL(6) = jMO & used for HGOH+hv->HG0+OH
Gstr{i,1}='HgOH'; 
fHgOH(i)=fHgOH(i)-1.0; fHg0(i)=fHg0(i)+1.0; fOH(i)=fOH(i)+1.0;

i=i+1;
Rnames{i}='HgOHNO2+hv=0.9HgOHO+0.1HgOH+0.9NO+0.1NO2';
k(:,i)=jMBNO2; % PHOTOL(20) = jMBNO2 & used for HGOHNO2+hv->HGOH+PRODUCTS
Gstr{i,1}='HgOHNO2'; 
fHgOHNO2(i)=fHgOHNO2(i)-1.0; fHgOHO(i)=fHgOHO(i)+0.9; fHgOH(i)=fHgOH(i)+0.1; fNO(i)=fNO(i)+0.9; fNO2(i)=fNO2(i)+0.1;

i=i+1;
Rnames{i}='HgOHHO2+hv=0.25HgOHO+0.67Hg0+0.07HgOH+0.01HgOHOH+0.67Cl+0.74HO2+0.26OH';
k(:,i)=jMBHO2; % PHOTOL(21) = jMBHO2 & used for HGOHHO2+hv->HGOH+PRODUCTS
Gstr{i,1}='HgOHHO2'; 
fHgOHHO2(i)=fHgOHHO2(i)-1.0; fHgOHO(i)=fHgOHO(i)+0.25; fHg0(i)=fHg0(i)+0.67; fHgOH(i)=fHgOH(i)+0.07; fHgOHOH(i)=fHgOHOH(i)+0.01; fCl(i)=fCl(i)+0.67; fHO2(i)=fHO2(i)+0.74; fOH(i)=fOH(i)+0.26;

i=i+1;
Rnames{i}='HgOHOH+hv=Hg0+2OH';
k(:,i)=jMC2; % PHOTOL(22) = jMC2 & used for HGOHOH+hv->HGOH+PRODUCTS
Gstr{i,1}='HgOHOH'; 
fHgOHOH(i)=fHgOHOH(i)-1.0; fHg0(i)=fHg0(i)+1.0; fOH(i)=fOH(i)+2.0;

i=i+1;
Rnames{i}='HgOHBrO+hv=HgOH+BrO';
k(:,i)=jMBBRO; % PHOTOL(23) = jMBBRO & used for HGOHBRO+hv->HGOH+PRODUCTS
Gstr{i,1}='HgOHBrO'; 
fHgOHBrO(i)=fHgOHBrO(i)-1.0; fHgOH(i)=fHgOH(i)+1.0; fBrO(i)=fBrO(i)+1.0;

i=i+1;
Rnames{i}='HgOHClO+hv=HgOH+ClO';
k(:,i)=jMBCLO; % PHOTOL(24) = jMBCLO & used for HGOHCLO+hv->HGOH+PRODUCTS
Gstr{i,1}='HgOHClO'; 
fHgOHClO(i)=fHgOHClO(i)-1.0; fHgOH(i)=fHgOH(i)+1.0; fClO(i)=fClO(i)+1.0;

i=i+1;
Rnames{i}='Hg2ORGP+hv=Hg0';
k(:,i)=jNO2_a; % PHOTOL(13) = jNO2_a & used for HG2ORGP+hv->HG0+PRODUCTS
Gstr{i,1}='Hg2ORGP'; 
fHg2ORGP(i)=fHg2ORGP(i)-1.0; fHg0(i)=fHg0(i)+1.0;

