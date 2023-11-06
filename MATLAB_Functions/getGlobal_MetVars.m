function [temp, press,numden, temp_over_k300, k300_over_temp]= getGlobal_MetVars()

    global TEMP PRESS NUMDEN TEMP_OVER_K300 K300_OVER_TEMP
    temp = TEMP;
    press = PRESS;
    numden = NUMDEN;
    temp_over_k300 = TEMP_OVER_K300;
    k300_over_temp = K300_OVER_TEMP;
end