
function [CNo_Emanuela, CNo_bluebook, CNo_SNV, CNo_MM, CNo_Bea] = trackcarrPLL_DLL_testjitter_10ms(file_in,f_sampling,prn,skp_msec,phaseFLL,mycarrfreq,msec,plotme);

% Code/Carrier Tracking 
% M. Pini on the basis of D. Akos - UCBoulder - ASEN5519 - Spring 2005; (c) D. Akos
%
%Inputs
%  prn - the prn number of the SV to track
%  mycodephase - the starting sample number in the data file to start tracking at (from acq)
%  mycarrfreq - the starting carrier frequency for the SV to track
%  msec - the number of msec to run for (be sure there is enough data for this)
%  plotme - does a plot of the output if "1", otherwise no plot is drawn
%
%Outputs
%  e_i,e_q,p_i,p_q,l_i,l_q - the 6 output accumulator values recorded once each code period corresponding
%                            to early, late, and prompt (e, l, p) and inphase and quadrature (i, q)
%
%Requires
% m-file: cacode.m (to generate the C/A code, 1 sample/chip)
% input data file named "gpsdata.bin" of signed char type data format  "gpsdata_long.bin" of signed char type data format
% "numsamp.bin" of short type data format

% Open this file only if you want to save the Prompt correlator output
% fid_1=fopen('C:\Marco_Pini\ISMB_2007\IRGAL\Tech\CNo Estimation\CNo_Data_sets_08Gennaio2008\CNo_test_IP.bin','w');
% fid_2=fopen('C:\Marco_Pini\ISMB_2007\IRGAL\Tech\CNo Estimation\CNo_Data_sets_08Gennaio2008\CNo_test_QP.bin','w');

%define sampling frequency
fs=f_sampling;


%skip through that data file to start at the appropriate sample (corresponding to code phase)
skp_msec = rem(skp_msec,20);
fseek(file_in,((phaseFLL-1)+round((20-skp_msec)*f_sampling*1e-3)),'bof');
% IMPORTANT: mycodephase: number of samples skiped during the acquisition  and processed by the FLL from the beginning of file;
%            round((20-cnt_skp)*f_sampling*1e-3) is the number of samples to skip, in order to start the phase tracking at the beginning of the next data bit.

%get a vector with the C/A code sampled 1x/chip
ca=cacode(prn);
%then make 1st value be the last chip and the last value the very first chip
ca=[ca(1023) ca ca ca ca ca ca ca ca ca ca ca(1)];

% -------------------------------------------------------------------------
% perform various initializations
%loop counter for times through the loop
loopcnt=1;

% define initial code frequency and code frequency basis of NCO in chips/sec
codefreq = 1.0230e6;
% NCO values initialization
codefreq_basis = 1.023e6;
%define number of chips in a code period
numchips=1023*10;
%define code phase (in chips) which is "left over"
remcodephase = 0.0;
%define early-late offset (in chips)
earlylate = 0.5;

%define carrier frequency which is used over whole tracking period (ok if apprx as noncoherent DLL)
carrfreq = mycarrfreq;
carrierfreq_basis = mycarrfreq; 
%define how much carrier phase is "left over" at the end of each code period to reset trig argument
remcarrphase = 0.0;
%--------------------------------------------------------------------------


% CARRIER TRACKING LOOP
zeta_carr = 1/sqrt(2);      % Butterworth 
Bl_carr =25;
wn_carr = Bl_carr/0.53;            % Natural Frequency;
k_carr =10;                 % Gain of the overall loop
% loop filter coefficients
t1_carr = k_carr/(wn_carr*wn_carr);
t2_carr = 2*zeta_carr/wn_carr; 
% filter values initialization
t2_div_t1_carr = t2_carr/t1_carr;
olderrort_carr = 0;
oldcarrier_nco = 0;
delta_carr = (10e-3)*1;               % integration time
delta_div_t1_carr = delta_carr/t1_carr;
vec_freq_carr = [];
vec_errort_carr = [];
vec_nco_carr =[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CODE TRACKING LOOP
zeta = 1/sqrt(2);      % Butterworth 
Bl = 2;
wn=Bl/0.53;              % Natural Frequency;
k = 1;                 % Gain of the overall loop

% loop filter coefficients
t1 = k/(wn*wn);
t2 = 2*zeta/wn; 
% filter values initialization
t2_div_t1 = t2/t1;
olderrort_code = 0;
oldcode_nco = 0;
delta = (10e-3)*1;               % integration time
delta_div_t1 = delta/t1;
blksize_vect = [];
Test_vect = [];
Discrim_out = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DEBUG %
global CNo_Emanuela;
CNo_Emanuela = [];
global CNo_bluebook;
CNo_bluebook = [];
global CNo_SNV;
CNo_SNV = [];
global CNo_MM;
CNo_MM = [];
global CNo_Bea;
CNo_Bea = [];
%start a timer
tic;

% while less than the number of specified code periods, continue processing
while (loopcnt <=  msec/10)
    
    %since it can be time consuming, do a periodic display of current stage
    %also a good point to display debugging information
    if (rem(loopcnt,10)==0)
        disp(['   Currently on interation:',int2str(loopcnt*10),' of:',int2str(msec)])
    end
    
    %update the phasestep based on code freq (variable) and sampling frequency (fixed)
    codephasestep=codefreq/fs;
    vet_freq(loopcnt) = codefreq;
    
    %find the size of a "block" or code period in whole samples
    %note - not starting from zero, but rather where you left off from last code period
    %this could be done using the CA code generator with sampling freq, code freq, and offset inputs
    %but this will be faster
    blksize=ceil((numchips-remcodephase)/codephasestep);
    blksize_vect(loopcnt) = blksize;
    %read in the appropriate number of samples to process this interation
    [rawdata,scount] = fread(file_in,blksize,'schar');
    %rawdata=rawdata'-128;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                        BLANKING
%  index = find((rawdata>2) | (rawdata>-3));
%  rawdata(index)=0;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %if did not read in enough samples, then could be out of data - better exit
    if (scount ~= blksize)
        disp('Not able to read the specified number of samples, exiting...')
        e_i=0.0;e_q=0.0;p_i=0.0;p_q=0.0;l_i=0.0;l_q=0.0;
        fclose(file_in);
        return
    end
    %define index into early code vector
    tcode=(remcodephase-earlylate):codephasestep:((blksize-1)*codephasestep+remcodephase-earlylate);
    tcode2=ceil(tcode)+1;
    earlycode=ca(tcode2);
    %define index into late code vector
    tcode=(remcodephase+earlylate):codephasestep:((blksize-1)*codephasestep+remcodephase+earlylate);
    tcode2=ceil(tcode)+1;
    latecode=ca(tcode2);
    %define index into prompt code vector
    tcode=remcodephase:codephasestep:((blksize-1)*codephasestep+remcodephase);
    tcode2=ceil(tcode)+1;
    promptcode=ca(tcode2);
    %now compute the remainder for next time around
    remcodephase_old = remcodephase;
    
    remcodephase = (tcode(blksize) + codephasestep) - 1023.0*10;
    
    Test_vect(loopcnt) = (blksize-1)*(1/fs)  + (codephasestep-remcodephase)*(1/codefreq) + remcodephase_old*(1/codefreq);
    %Test_vect(loopcnt) = (blksize-1)*(1/fs)  - remcodephase*(1/codefreq);
    
    %generate the carrier frequency to mix the signal to baseband
    time=(0:blksize) ./ fs;
    %get the argument to sin/cos functions
    trigarg = ((carrfreq * 2.0 * pi) .* time) + remcarrphase;
    %compute the "leftover" on sample after the last to start there next time
    remcarrphase=rem(trigarg(blksize+1),(2 * pi));
    %finally compute the signal to mix the collected data to bandband
    carrcos=cos(trigarg(1:blksize));
    carrsin=sin(trigarg(1:blksize));
    
    %generate the six standard accumulated values
    %first mix to baseband
    tempdatacos = carrcos .* rawdata';
    tempdatasin = carrsin .* rawdata';
    
    %now get early, late, and prompt values for each
    e_i(loopcnt) = sum(earlycode .* tempdatasin);
    e_q(loopcnt) = sum(earlycode .* tempdatacos);
    p_i(loopcnt) = sum(promptcode .* tempdatasin);
    p_q(loopcnt) = sum(promptcode .* tempdatacos);
    l_i(loopcnt) = sum(latecode .* tempdatasin);
    l_q(loopcnt) = sum(latecode .* tempdatacos);
    
    
    % PLL discriminator and feed-back 
    %implement carrier loop discriminator (phase detector)
    errort_carr = atan ( p_q(loopcnt)/p_i(loopcnt) );    % discriminator as arc-tan
    
    vec_errort_carr(loopcnt) = errort_carr;        
    %implement code loop filter and generate NCO command
    carrier_nco = oldcarrier_nco + (t2_div_t1_carr * (errort_carr - olderrort_carr)) + errort_carr*delta_div_t1_carr;
    %modify code freq based on NCO command
    %delta_freq = code_nco - oldcode_nco;
    carrfreq =  mycarrfreq + carrier_nco;
    % N.B! freq = d(phase)/dt; 
    vec_nco_carr(loopcnt) = carrier_nco; 
    vec_freq_carr(loopcnt) = carrfreq;
    olderrort_carr = errort_carr;
    oldcarrier_nco = carrier_nco; 
    
    %-------------------------------
    
    % DLL discriminator and feed-back
    errort_code =0.5*((( e_i(loopcnt)^2 + e_q(loopcnt)^2 ) - ( l_i(loopcnt)^2 + l_q(loopcnt)^2 ))/(( e_i(loopcnt)^2 + e_q(loopcnt)^2 ) + (l_i(loopcnt)^2 + l_q(loopcnt)^2 )));
    %errort_code = ((e_i(loopcnt) - l_i(loopcnt))*sign(p_i(loopcnt)))/(p_i(loopcnt));
    % N.B: Normalized discriminator, the gain of the overall loop is set to 1 (k=1) 
    Discrim_out(loopcnt) = errort_code;
    
    %implement code loop filter and generate NCO command
    code_nco = oldcode_nco + (t2_div_t1 * (errort_code - olderrort_code)) + errort_code*delta_div_t1;
    
    %modify code freq based on NCO command
    codefreq = codefreq_basis - code_nco;
    
    olderrort_code = errort_code;
    oldcode_nco = code_nco; 
    vec_errort_code(loopcnt) = errort_code;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     fwrite(fid_1,p_i(loopcnt),'double');
%     fwrite(fid_2,p_q(loopcnt),'double');
    
%    PLL_st = control_PLL_v2_EF(p_i(loopcnt),p_q(loopcnt),loopcnt,35);
    
%     if (PLL_st ~= 0) 
%         disp('PLL not locked, C/No below the threshold, the channel might be reassigned to another satellite, or try with a new fast acquisition...')
%         break;
%     end;
    
    %increment loop counter
    loopcnt=loopcnt+1;
    
end
fclose(file_in);

disp(['Time required to process ',int2str(msec),' msec of data is:',num2str(toc),'sec']);

fclose('all')

if (plotme == 1)
    
    x_axe_second = 0.01:0.01:length(p_i)/100;
    
    
    subplot(121),plot(p_i,'o-')
    grid
    xlabel('milliseconds')
    ylabel('inphase amplitude')
    title('Prompt Inphase Correlation Results')
    subplot(122),plot(p_q,'o-')
    grid   
    xlabel('milliseconds')
    ylabel('quadrature amplitude')
    title('Prompt Quadrature Correlation Results')
    
    figure,
    plot(p_i,p_q,'.'),grid on;
    
    for (i=1:length(vec_freq_carr))
        if (i<=500),
            vect_local = vec_freq_carr(1:i);
            vec_freq_mean(i) = sum(vect_local)/length(vect_local);
        else
            vect_local = vec_freq_carr(i-499:i);
            vec_freq_mean(i) = sum(vect_local)/length(vect_local);
        end
    end     
    figure
    plot (vec_freq_carr - mycarrfreq,'r'),grid on;
    hold on
    plot (vec_freq_mean - mycarrfreq,'g.-');
    xlabel('milliseconds');
    ylabel('Frequency [Hz]')
    title('Increment of the local carrier frequency with respect to the nominal value')
    %    figure
    %    plot (vec_errort,'k'),grid on;
    disp(['Time required to process ',int2str(msec),' msec of data is:',num2str(toc),'sec']);
    mean(vet_freq);
    
    figure(100)
    subplot(212),plot(x_axe_second, (p_i .^2 + p_q .^ 2), 'g.-'),grid on;
    hold on
    grid
    subplot(212),plot(x_axe_second, (e_i .^2 + e_q .^ 2), 'bx-');
    subplot(212),plot(x_axe_second, (l_i .^2 + l_q .^ 2), 'r+-');
    xlabel('milliseconds');
    ylabel('amplitude');
    title('Correlation Results');
    legend('prompt','early','late');
    
    x_axe_second_CNo = 0.5:0.5:length(CNo_bluebook)/2;
    subplot(211), plot(x_axe_second_CNo,CNo_bluebook,'k.-'),grid on;
    ylabel('C/No[dBHz]');
    
    figure(101),
    subplot(211),plot(x_axe_second,vec_errort_carr,'k*-'),hold on,grid on
    title('Phase residual error');
    subplot(212),plot(x_axe_second,(e_i-l_i)./(2*p_i),'k-'),hold on,grid on
    title('Pseudorange tracking error [chip]');
    
    figure
    plot(vet_freq,'k.-'), grid on;
    xlabel('milliseconds')
    ylabel('Frequency [Hz]')
    title('Local Code Frequency')
    mean_freq = mean(vet_freq);
    
        
    figure
    plot(vec_errort_carr,'r*-'),hold on,grid on
    title('Phase residual error')
end



