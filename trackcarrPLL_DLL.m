function trackcarrPLL_DLL(file_in,prn,skp_msec,phaseFLL,mycarrfreq,msec,plotme);
global sampletype;
% Code/Carrier Tracking
% M. Pini on the basis of D. Akos - UCBoulder - ASEN5519 - Spring 2005; (c) D. Akos

%Inputs
%  prn - the prn number of the SV to track
%  mycodephase - the starting sample number in the data file to start tracking at (from acq)
%  mycarrfreq - the starting carrier frequency for the SV to track
%  msec - the number of msec to run for (be sure there is enough data for this)
%  plotme - does a plot of the output if "1", otherwise no plot is drawn

%Outputs
%  e_i,e_q,p_i,p_q,l_i,l_q - the 6 output accumulator values recorded once each code period corresponding
%                            to early, late, and prompt (e, l, p) and inphase and quadrature (i, q)

%Requires
% m-file: cacode.m (to generate the C/A code, 1 sample/chip)
% input data file named "gpsdata.bin" of signed char type data format  "gpsdata_long.bin" of signed char type data format
% "numsamp.bin" of short type data format

global isNewRun;
global f_sampling;
if isNewRun
    fseek(file_in,0,-1);

    %skip through that data file to start at the appropriate sample (corresponding to code phase)
    skp_msec = rem(skp_msec,20);
    fseek(file_in,((phaseFLL-1)+round((20-skp_msec)*f_sampling*1e-3))*sampletype,'bof');
    % IMPORTANT: mycodephase: number of samples skiped during the acquisition  and processed by the FLL from the beginning of file;
    %            round((20-skp_msec)*f_sampling*1e-3) is the number of samples to skip, in order to start the phase tracking at the beginning of the next data bit.


    %skip through that data file to start at the appropriate sample (corresponding to code phase)
    %fseek(file_in,mycodephase-1,'bof');

    %get a vector with the C/A code sampled 1x/chip
%     ca=generateCAcode(prn);
ca = (double(gnssCACode(prn,"NavIC L5-SPS")-0.5))';
    %then make 1st value be the last chip and the last value the very first chip
    ca=[ca(1023) ca ca(1)];

    % -------------------------------------------------------------------------
    % perform various initializations
    %loop counter for times through the loop
    loopcnt=1;
    %define sampling frequency
    fs=f_sampling;

    % define initial code frequency and code frequency basis of NCO in chips/sec
    codefreq = 1.0230e6;
    % NCO values initialization
    codefreq_basis = 1.023e6;
    %define number of chips in a code period
    numchips=1023;
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
    Bl_carr =15;
    wn_carr = Bl_carr/0.53;            % Natural Frequency;
    k_carr = 1;                 % Gain of the overall loop
    % loop filter coefficients
    t1_carr = k_carr/(wn_carr*wn_carr);
    t2_carr = 2*zeta_carr/wn_carr;
    % filter values initialization
    t2_div_t1_carr = t2_carr/t1_carr;
    olderrort_carr = 0;
    oldcarrier_nco = 0;
    delta_carr = 1e-3;               % integration time
    delta_div_t1_carr = delta_carr/t1_carr;
    vec_freq_carr = [];
    vec_errort_carr = [];
    vec_nco_carr =[];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % CODE TRACKING LOOP
    zeta = 1/sqrt(2);      % Butterworth
    Bl=2;
    wn=Bl/0.53;              % Natural Frequency;
    k = 1;                 % Gain of the overall loop

    % loop filter coefficients
    t1 = k/(wn*wn);
    t2 = 2*zeta/wn;
    % filter values initialization
    t2_div_t1 = t2/t1;
    olderrort_code = 0;
    oldcode_nco = 0;
    delta = 1e-3;               % integration time
    delta_div_t1 = delta/t1;
    errort_code = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % while less than the number of specified code periods, continue processing
    while (loopcnt <=  msec)

        %since it can be time consuming, do a periodic display of current stage
        %also a good point to display debugging information
        if (rem(loopcnt,100)==0)
            disp(['   Currently on iteration:',int2str(loopcnt),' of:',int2str(msec)])
        end

        %    ? update on the basis of the code loop
        %    %read in the blksize from the previous generated code tracking
        %    blksize = fread(fid2,1,'short');

        %read in the appropriate number of samples of data and prompt PRN code to process this interation
        %    [rawdata,scount] = fread(file_in,blksize,'schar');
        %    rawdata=rawdata';
        %
        %update the phasestep based on code freq (variable) and sampling frequency (fixed)
        codephasestep=codefreq/fs;
        vet_freq(loopcnt) = codefreq;

        %find the size of a "block" or code period in whole samples
        %note - not starting from zero, but rather where you left off from last code period
        %this could be done using the CA code generator with sampling freq, code freq, and offset inputs
        %but this will be faster
        blksize=ceil((numchips-remcodephase)/codephasestep);
        %read in the appropriate number of samples to process this interation
%         [rawdata,scount] = fread(file_in,sampletype*blksize,'schar');
        [rawdata,scount] = fread(file_in,sampletype*blksize,'int16');
        if (sampletype==2)
            data1=rawdata(1:2:end);    
            data2=rawdata(2:2:end);    
            rawdata=data1 + 1i*data2;
        end 
        %rawdata=rawdata'-127;
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % %                        BLANKING
        % aa = find(rawdata > 8);
        % rawdata(aa)= 0;
        % bb = find(rawdata < -7);
        % rawdata(bb)= 0;

        % %                      PARTIAL  BLANKING
        % segno = sign(rawdata);
        % data_bin = dec2bin(abs(rawdata),7);
        % data_bin(:,1:4) = '0';
        % rawdata = segno.*(bin2dec(data_bin))';

        % %                      PARTIAL  BLANKING 1
        aa = find(rawdata > 8);
        samples_noise = round(sqrt(5)*randn(1,length(rawdata)));
        rawdata(aa)= samples_noise(aa);

        bb = find(rawdata < -7);
        samples_noise = round(sqrt(5)*randn(1,length(rawdata)));
        rawdata(bb)= samples_noise(bb);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %if did not read in enough samples, then could be out of data - better exit
        if (scount ~= blksize*sampletype)
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
        tcode=(remcodephase+earlylate):codephasestep:((blksize-1)*codephasestep+remcodephase+earlylate);;
        tcode2=ceil(tcode)+1;
        latecode=ca(tcode2);
        %define index into prompt code vector
        tcode=remcodephase:codephasestep:((blksize-1)*codephasestep+remcodephase);
        tcode2=ceil(tcode)+1;
        promptcode=ca(tcode2);

        %now compute the remainder for next time around
        remcodephase = (tcode(blksize) + codephasestep) - 1023.0;

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
        tempdata=rawdata'.*(carrcos+1i*carrsin);
        tempdatacos=real(tempdata);
        tempdatasin=imag(tempdata);
        %tempdatacos = carrcos .* rawdata';
        %tempdatasin = carrsin .* rawdata';

        %now get early, late, and prompt values for each
        e_i(loopcnt) = sum(earlycode .* tempdatasin);
        e_q(loopcnt) = sum(earlycode .* tempdatacos);
        p_i(loopcnt) = sum(promptcode .* tempdatasin);
        p_q(loopcnt) = sum(promptcode .* tempdatacos);
        l_i(loopcnt) = sum(latecode .* tempdatasin);
        l_q(loopcnt) = sum(latecode .* tempdatacos);
        absoluteSample(loopcnt) = ftell(file_in);

        % PLL discriminator and feed-back
        %implement carrier loop discriminator (phase detector)
        errort_carr = atan ( p_q(loopcnt)/p_i(loopcnt) );    % discriminator as arc-tan

        %    if (loopcnt > 500),
        %         if (sign(p_q(loopcnt)) ~= sign(p_q(loopcnt-1)) )
        %             codephase_out = loopcnt*blksize + phaseFLL;
        %             break;
        %         end
        %     end

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
        errort_code = 0.5*(( e_i(loopcnt)^2 + e_q(loopcnt)^2 ) - ( l_i(loopcnt)^2 + l_q(loopcnt)^2 ))/(( e_i(loopcnt)^2 + e_q(loopcnt)^2 ) + ( l_i(loopcnt)^2 + l_q(loopcnt)^2 ));
        %errort_code = ((e_i(loopcnt) - l_i(loopcnt))*sign(p_i(loopcnt)))/(p_i(loopcnt));
        % N.B: Normalized discriminator, the gain of the overall loop is set to 1 (k=1)


        %implement code loop filter and generate NCO command
        code_nco = oldcode_nco + (t2_div_t1 * (errort_code - olderrort_code)) + errort_code*delta_div_t1;

        %modify code freq based on NCO command
        codefreq = codefreq_basis - code_nco;

        olderrort_code = errort_code;
        oldcode_nco = code_nco;
        vec_errort_code(loopcnt) = codefreq;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Compute C
        a = abs(p_i(loopcnt))/16368;
        C(loopcnt) = (a*2/sqrt(2))^2;  % Note: p_i*2 is due to the loss I have in the cos multipliation. .../sqrt(2) is needed since I don't want
        % the peak value in the computation of C, but the effective value instead
        %increment loop counter
        loopcnt=loopcnt+1;

    end

    %disp(['Time required to process ',int2str(msec),' msec of data is:',num2str(toc),'sec'])

    if (plotme == 1)
%         close all
        subplot(121),plot(p_i)
        grid
        xlabel('milliseconds')
        ylabel('inphase amplitude')
        title('Prompt Inphase Correlation Results')
        subplot(122),plot(p_q)
        grid
        xlabel('milliseconds')
        ylabel('quadrature amplitude')
        title('Prompt Quadrature Correlation Results')

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
        title({prn; 'Increment of the local carrier frequency with respect to the nominal value'})
        %    figure
        %    plot (vec_errort,'k'),grid on;
        %disp(['Time required to process ',int2str(msec),' msec of data is:',num2str(toc),'sec']);
        mean(vet_freq)
        
       NBP = []
       WBP = []
% for k=1:msec"
for i=1:20:41
    NBP= [NBP sum(p_i(i:i+19))^2+sum(p_q(i:i+19))^2];
    WBP= [WBP sum(p_i(i:i+19).^2)+sum(p_q(i:i+19).^2)];
end
    SIraw=NBP-WBP;
    SItrend=mean(SIraw);
    SI=SIraw./SItrend;

    S4= sqrt((mean(SI.^2)-mean(SI)^2)/mean(SI)^2);

        fid = fopen('s4.txt','a');
        fprintf(fid,'%d %f \n',prn,S4);
        fclose(fid);
        
        
        
        figure
        plot(p_i .^2 + p_q .^ 2, 'g.-')
        hold on
        grid
        plot(e_i .^2 + e_q .^ 2, 'bx-');
        plot(l_i .^2 + l_q .^ 2, 'r+-');
        xlabel('milliseconds')
        ylabel('amplitude')
        title({prn; 'Correlation Results'})
        legend('prompt','early','late')
        figure
        plot(vet_freq,'k.-'), grid on;
        xlabel('milliseconds')
        ylabel('Frequency [Hz]')
        title('Local Code Frequency')
        mean_freq = mean(vet_freq);
    end
    mat_file = ['./GPS_TrackPLL_PRN',num2str(prn)];
    save (mat_file);
else
    mat_file = ['./GPS_TrackPLL_PRN',num2str(prn)];
    load(mat_file);
end


