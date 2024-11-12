function data_saving_for_demodulation(PRN_inview)
% store data to channels' structure
% channels(sat_i).PRN
%                .symbol_vect: data symbols
%                .absoluteSample: absolute location in input file (byte) equivalent to a
%                                 data symbol
channels = [];
for iii = 1:length(PRN_inview)
    PRN = PRN_inview(iii);
    mat_file = ['./GPS_TrackPLL_PRN',num2str(PRN)];
    load(mat_file);    
    channels(iii).PRN = PRN;
    symbol_vect(p_i > 0)  =  1;
    symbol_vect(p_i <= 0) = -1;
    channels(iii).symbol_vect = symbol_vect;
    channels(iii).absoluteSample = absoluteSample;
end
save('./GPS_channels', 'channels')