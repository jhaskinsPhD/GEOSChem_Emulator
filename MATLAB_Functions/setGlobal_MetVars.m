function setGlobal_MetVars(Met)
    struct2var(Met)
    global TEMP PRESS NUMDEN TEMP_OVER_K300 K300_OVER_TEMP

    TEMP=T;   % [molecules./cm3]
    PRESS=P;  % [mbar]
    NUMDEN=M; % [molecules./cm3]
    TEMP_OVER_K300 = T./300;
    K300_OVER_TEMP= 300./T;
end




