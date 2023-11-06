WHAT IS 'JMAP_TYPE'? 

	The string you enter for the input 'jmap_type' determines how we map the
	photolysis frequencies that are calculated in the FJX module of GEOS-Chem 
	to a photolysis frequency already defined in F0AM. Most of the photolysis frequencies 
	currently defined in F0AM are MCM photolysis rates, but there is additional info on other
	rates (including for halogens). 

WHAT ARE THE ACCEPTED 'JMAP_TYPE' INPUT ARGUEMENTS? 

	All 'jmap_type' inputs are case insensitive. The allowed options are as follows: 

	(1) 'CrossSec_Match'  - Matches j-values from GEOS-Chem to MCM using the best 
							match to the cross-section that defines that j-value in
							the GEOS-Chem FJX version correponding to your GC_version. 
							This is the DEFAULT that is used if no input for 'jamp_type'
							is given because it will give you the MOST analaogous 
							chemistry in F0AM compared to what is in GEOS-Chem. 
										 
	(2) 'MCM'             - Matches j-values from GEOS-Chem to MCM using their BEST match 
							in the MCM (based on the compound), regardless of whether FJX & 
							GEOS-Chem use that cross-section to define that j-value. For example, 
							the photolysis frequencies of XXX compound in GEOS-Chem might be 
							defined using the cross section of YYY compound owing to a lack of 
							data about that specific compound's photolysis, but the MCM might 
							contain a specific photolysis frequency JUST for compound XXX since it 
							is generally more detailed than GEOS-Chem is, chemically. 							
							So, you might want to set 'jmap_type' to 'MCM' in order to use the
							MCM photolysis frequency since its more specific to compound XXX 
							than just assigning it the MCM analoge of the j-value for compound YYY. 
				   
	(3) 'RCIM'            - Matches j-values from GEOS-Chem to MCM following the
							methodology that the Reduced Caltech Isoprene Mechanism (RCIM)
							uses, which uses "more realistic" J-Values that are pre-defined 
							in F0AM, and does differ from those defined in the MCM. For example, 
							the photolysis of MVK might be defined in the MCM, but have a newer, 
							more updated photolysis frequency in the RCIM. So you might choose to 
							set 'jmap_type' equal to 'RCIM' in order to have a more realistic/updated 
							version of the phototlysis of MVK to use. See the Wenneberg et al., 2019 
							Supplement for more details. 
              
	NOTE: Both 'RCIM' and 'MCM' options will not replicate perfectly GEOS-Chem's 
		  chemical mechanism, as they contain more detail for the j-values than GEOS-Chem does, which is why 
		  the DEFAULT is to do the 'CrossSec_Match' match based off the cross-section used for each j-value 
		  in the FJX file corresponding to your GEOS-Chem version.  