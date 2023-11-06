function [set, get]=  import_global_functs()
    % Wrapper Function to import handles for setting and getting global met vars.  
    set=@setGlobal_MetVars; 
    get=@getGlobal_MetVars; 
end 


function setGlobal_MetVars(Met)
    struct2var(Met)
    global TEMP PRESS NUMDEN TEMP_OVER_K300 K300_OVER_TEMP

    TEMP=T;   % [molecules./cm3]
    PRESS=P;  % [mbar]
    NUMDEN=M; % [molecules./cm3]
    TEMP_OVER_K300 = T./300;
    K300_OVER_TEMP= 300./T;
end

function [temp, press,numden, temp_over_k300, k300_over_temp]= getGlobal_MetVars()

    global TEMP PRESS NUMDEN TEMP_OVER_K300 K300_OVER_TEMP
    temp = TEMP;
    press = PRESS;
    numden = NUMDEN;
    temp_over_k300 = TEMP_OVER_K300;
    k300_over_temp = K300_OVER_TEMP;
end




