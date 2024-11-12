clear all;

global f_sampling; f_sampling = 5e6; % sampling frequency [Hz]
global nominalfreq; nominalfreq = 0; % IF frequency [Hz]
global sampletype; sampletype = 2; 
global code_rate; code_rate =  1.023e6;
global code_length; code_length = 1023;

fid1=fopen('C:\teleorbit\simout1.o','rb');
fid2=fopen('C:\teleorbit\simout2.o','rb');

fidout=fopen('C:\teleorbit\out.bin','w');
fclose(fidout);

fidout=fopen('C:\teleorbit\out.bin','a');

n_code_extr = 1; %s

while (~feof(fid1) & ~feof(fid2))
[gpsdata1,scount1] = fread(fid1,n_code_extr*sampletype*f_sampling,'int16');
[gpsdata2,scount2] = fread(fid2,n_code_extr*sampletype*f_sampling,'int16');

% fwrite(fidout,0.96*gpsdata1+0.04*gpsdata2,'int16');
fwrite(fidout,0.5*gpsdata1+0.5*gpsdata2,'int16');
end
fclose('all');