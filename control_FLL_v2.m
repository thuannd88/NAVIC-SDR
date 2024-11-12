
% Control FLL.
% This program control the FLL status, checking the freqeuncy of the local
% carrier frequency, passed from the FLL program.
% The check is based on a FIFO structure and is performed at a predefined
% control time, which could be for example every 0.25 s.

function [FLL_status,freq_control] = control_FLL_v2(act_freq,loop,thres);


global Freq_sum;
global Old_Freq_sum;

Freq_sum = (1/20)*act_freq + ((20-1)/20)*Freq_sum;

FLL_status = 0;
freq_control = 0;

if (rem(loop,100)==0)
    
    if (abs((Freq_sum - Old_Freq_sum))< 1)
        FLL_status = 1;
        freq_control = Freq_sum;
    end
   % figure(100),hold on, plot(loop, Freq_sum,'o');
    Old_Freq_sum = Freq_sum;
end


if(loop >= 5000)       % After 5 seconds (5000 msec, since the PDI used in the FLL is equal to 1 msec), the time is out! FLL not locked
    FLL_status = -1;
end

return