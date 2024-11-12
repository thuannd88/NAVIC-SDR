% Marco Pini
% GPS Receiver Architectures
% Assignment 1 - 01/18/05

function [mysamples] = cacode(svnum,fs,coderate,numsamp,offset)

% Input Arguments:
% svnum - the Satellite's PRN number
% fs - sampling frequency
% coderate - PRN code rate
% numsamp - number of samples in the output vector
% offset - offset beetwen the code and the sampling frequency
% 
% Output Argument
% mysamples - a 1023 element vector containing the desired output sequence
%
% D. Akos - UCBoulder - ASEN5519 - Spring 2005

% the g2s vector holds the appropriate shift of the g2 code to generate
% the C/A code (ex. for SV#19 - use a G2 shift of g2shift(19,1)=471)
g2s = [5;6;7;8;17;18;139;140;141;251;252;254;255;256;257;258;469;470;471; ...
       472;473;474;509;512;513;514;515;516;859;860;861;862];

g2shift=g2s(svnum,1);
           
% Generate G1 code
     %   load shift register
          reg = -1*ones(1,10);
     %
     for i = 1:1023,
         g1(i) = reg(10);
         save1 = reg(3)*reg(10);
         reg(1,2:10) = reg(1:1:9);
         reg(1) = save1;
     end
%
% Generate G2 code
%
     %   load shift register
           reg = -1*ones(1,10);
     %
     for i = 1:1023,
         g2(i) = reg(10);
         save2 = reg(2)*reg(3)*reg(6)*reg(8)*reg(9)*reg(10);
         reg(1,2:10) = reg(1:1:9);
         reg(1) = save2;
     end
     %
     %    Shift G2 code
     %
     g2tmp(1,1:g2shift)=g2(1,1023-g2shift+1:1023);
     g2tmp(1,g2shift+1:1023)=g2(1,1:1023-g2shift);
     %
     g2 = g2tmp;
%
%  Form single sample C/A code by multiplying G1 and G2 point by point
%
% ca = g1.*g2;
ca = GNSScodegen(svnum,'L5I',0);



% ----- Inizialization Section
mysamples = [];                 % output vector
Ts = 1/fs;                      % duration of the sampling period
Tc = 1/coderate;                % duration of one chip

Tc_norm = 1;                    % normalized duration of one chip
Ts_norm = (Tc_norm*Ts)/Tc;      % normalized sampling time
TC_flag_norm = 1;
% --------------------
% ------ Offset Control Section
if (offset < 0),
    indexcode = 1023
    epoc = 1 - abs(offset);
else
indexcode = 1;                 % C/A code index
epoc = offset;                 % sampling instant
end
% -------------------- 
% ------ Core Section
for (index=1:numsamp),    
    
    if epoc >= TC_flag_norm,               % a new chip must be sampled
        indexcode = indexcode + 1;
        TC_flag_norm = TC_flag_norm + 1;
            if indexcode > 1023,           % the new chip is the first of the C/A code frame
              indexcode = 1;
            end
    end   

mysamples(index) = ca(indexcode);         % sampling process
epoc = epoc + Ts_norm;                    % future sampling instant

end
return;
% ----------------