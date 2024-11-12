function [e_i,e_q,p_i,p_q,l_i,l_q,absoluteSample] = tracking_GLO(fid,FCH,mycodephase,mycarrfreq,msec,plotme)

% Code/Carrier Tracking
% M. Pini on the basis of D. Akos - UCBoulder - ASEN5519 - Spring 2005; (c) D. Akos

%Inputs
%  FCH - the channel number of the SV to track
%  mycodephase - the starting sample number in the data file to start tracking at (from acq)
%  mycarrfreq - the starting carrier frequency for the SV to track
%  msec - the number of msec to run for (be sure there is enough data for this)
%  plotme - does a plot of the output if "1", otherwise no plot is drawn

%Outputs
%  e_i,e_q,p_i,p_q,l_i,l_q - the 6 output accumulator values recorded once each code period corresponding
%                            to early, late, and prompt (e, l, p) and inphase and quadrature (i, q)


global isNewRun;
global frf0;
global fch_step;
global f_sampling;
global nominalfreq;
global PRN;
global code_rate;
global code_length;
global seek_sec;
global sampletype;

if (isNewRun==1)

	currentSample = ((mycodephase-1)+(seek_sec*f_sampling))*floor(sampletype);
% fprintf('\n %f \n',currentSample);
    %skip through that data file to start at the appropriate sample (corresponding to code phase)
	fseek(fid,currentSample,-1);
	%skp_msec = rem(skp_msec,20);
	%fseek(fid,((mycodephase-1)+round((20-skp_msec)*f_sampling*1e-3)),'bof');
    %IMPORTANT: round((20-skp_msec)*f_sampling*1e-3) is the number of samples to skip, in order to start the phase tracking at the beginning of the next data bit.
    %fseek(fid,mycodephase-1,'bof');

    fif=nominalfreq+FCH*0.5625e6;
    %then make 1st value be the last chip and the last value the very first chip
    ca=[PRN(end) PRN PRN(1)];

    % -------------------------------------------------------------------------
    % perform various initializations
    %loop counter for times through the loop
    loopcnt=1;
    %define sampling frequency
    fs=f_sampling;

    % define initial code frequency and code frequency basis of NCO in chips/sec
    codefreq = code_rate;
    %define code phase (in chips) which is "left over"
    remcodephase = 0.0;
    %define early-late offset (in chips)
    earlylate = 0.05;
    %define carrier frequency which is used over whole tracking period (ok if apprx as noncoherent DLL)
    carrfreq = mycarrfreq;
    %define how much carrier phase is "left over" at the end of each code period to reset trig argument
    remcarrphase = 0.0;
	    
    %--------------------------------------------------------------------------
    % CARRIER TRACKING LOOP (FLL-PLL integration)
    PDIcarr = 0.001;
	pllbw = 25;
	fllbw = 250;
	k1 = PDIcarr*((pllbw/0.53)^2) + sqrt(2)*(pllbw/0.53);
	k2 = sqrt(2)*(pllbw/0.53);
	k3 = PDIcarr*(fllbw/0.25);

    % CODE TRACKING LOOP
    PDIcode = 0.001;		% integration time
	zeta = 1/sqrt(2);		% Butterworth
    dllbw = 0.5;			%DLL noise bandwidth
    wn = dllbw*8*zeta/(4*zeta^2+1);		% Natural Frequency;
    k = 1;					% Gain of the overall loop
    % loop filter coefficients
    t1 = k/(wn*wn);
    t2 = 2*zeta/wn;
    % filter values initialization
    t2_div_t1 = t2/t1;
    olderr_code = 0;
    oldnco_code = 0;
	olderr_carr = 0;
    oldnco_carr = 0;   
    delta_div_t1 = PDIcode/t1;
    
    e_i = zeros(1, msec);
    e_q = zeros(1, msec);
    p_i = zeros(1, msec);
    p_q = zeros(1, msec);
    l_i = zeros(1, msec);
    l_q = zeros(1, msec);
    C = zeros(1, msec);
    vet_freq = zeros(1, msec);
    vec_nco_carr = zeros(1, msec);
    vec_err_carr = zeros(1, msec);
    vec_freq_carr = zeros(1, msec);
    vec_nco_code = zeros(1, msec);
    vec_err_code = zeros(1, msec);
    vec_freq_code = zeros(1, msec);
    absoluteSample = zeros(1, msec);
    vec_freq_mean = zeros(1, msec);
	I1 = 0.001;
	Q1 = 0.001;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % while less than the number of specified code periods, continue processing
    while (loopcnt <=  msec)
        %update the phasestep based on code freq (variable) and sampling frequency (fixed)
        codephasestep=codefreq/fs;
        vet_freq(loopcnt) = codefreq;

        %find the size of a "block" or code period in whole samples
        %note - not starting from zero, but rather where you left off from last code period
        %this could be done using the CA code generator with sampling freq, code freq, and offset inputs
        %but this will be faster
        blksize=ceil((code_length-remcodephase)/codephasestep);
        %read in the appropriate number of samples to process this interation
		[glodata,scount] = fread(fid,blksize*floor(sampletype),'schar');
		data1=glodata(1:2:end);    
		data2=glodata(2:2:end);
        switch (sampletype)
			case 2
				glodata=data1 + 1i*data2;
			case 2.5
				glodata=data1;
        end
        currentSample = currentSample + blksize*floor(sampletype);
% fprintf(' %d: %d ~ %d ; ',loopcnt,blksize,currentSample);
        %[rawdata,scount] = fread(fid,blksize,'schar');
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
        %aa = find(rawdata > 8);
        %samples_noise = round(sqrt(5)*randn(1,length(rawdata)));
        %rawdata(aa)= samples_noise(aa);

        %bb = find(rawdata < -7);
        %samples_noise = round(sqrt(5)*randn(1,length(rawdata)));
        %rawdata(bb)= samples_noise(bb);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %if did not read in enough samples, then could be out of data - better exit
        if (scount ~= blksize*floor(sampletype))
            fprintf('Not able to read the specified number of samples, exiting...');
            e_i=0.0;e_q=0.0;p_i=0.0;p_q=0.0;l_i=0.0;l_q=0.0;
            fclose(fid);
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
        remcodephase = (tcode(blksize) + codephasestep) - code_length;

        %generate the carrier frequency to mix the signal to baseband
        time=(0:blksize) ./ fs;
        %get the argument to sin/cos functions
        trigarg = ((carrfreq * 2.0 * pi) .* time) + remcarrphase;
        %compute the "leftover" on sample after the last to start there next time
        remcarrphase=rem(trigarg(blksize+1),(2 * pi));
        %finally compute the signal to mix the collected data to bandband
        sigCarr = exp(1i .* trigarg(1:blksize));

        %generate the six standard accumulated values
        %first mix to baseband
        tempdatacos = real(sigCarr.*(glodata.'));
        tempdatasin = imag(sigCarr.*(glodata.'));

        %now get early, late, and prompt values for each
        e_i(loopcnt) = sum(earlycode .* tempdatasin);
        e_q(loopcnt) = sum(earlycode .* tempdatacos);
        p_i(loopcnt) = sum(promptcode .* tempdatasin);
        p_q(loopcnt) = sum(promptcode .* tempdatacos);
        l_i(loopcnt) = sum(latecode .* tempdatasin);
        l_q(loopcnt) = sum(latecode .* tempdatacos);

		% Find combined PLL/FLL error and update carrier NCO (FLL-assisted PLL) ------
        I2 = I1;
		Q2 = Q1;
        I1 = p_i(loopcnt);
		Q1 = p_q(loopcnt);
        cross = I1*Q2 - I2*Q1;
        dot = abs(I1*I2 + Q1*Q2);

        % FLL-assisted PLL discriminators and feed-back
        %implement carrier loop discriminator (phase detector)
        err_carr = atan(p_q(loopcnt)/p_i(loopcnt))/(2*pi);    % discriminator as arc-tan
		%Implement carrier loop discriminator (frequency detector)
		err_freq = atan2(cross, dot)/pi;
        vec_err_carr(loopcnt) = err_carr;
		%implement carrier loop filter and generate NCO command
		nco_carr = oldnco_carr + k1*err_carr - k2*olderr_carr - k3*err_freq;
        %modify carrier freq based on NCO command
        carrfreq =  mycarrfreq + nco_carr;
        % N.B! freq = d(phase)/dt;
        vec_nco_carr(loopcnt) = nco_carr;
        vec_freq_carr(loopcnt) = carrfreq;
        olderr_carr = err_carr;
        oldnco_carr = nco_carr;
		
        % DLL discriminator and feed-back
        err_code = (sqrt( e_i(loopcnt)^2 + e_q(loopcnt)^2 ) - sqrt( l_i(loopcnt)^2 + l_q(loopcnt)^2 ))/(sqrt( e_i(loopcnt)^2 + e_q(loopcnt)^2 ) + sqrt( l_i(loopcnt)^2 + l_q(loopcnt)^2 ));
        %errort_code = ((e_i(loopcnt) - l_i(loopcnt))*sign(p_i(loopcnt)))/(p_i(loopcnt));
        % N.B: Normalized discriminator, the gain of the overall loop is set to 1 (k=1)
        vec_err_code(loopcnt) = err_code;
        %implement code loop filter and generate NCO command
        nco_code = oldnco_code + (t2_div_t1 * (err_code - olderr_code)) + err_code*delta_div_t1;
        %modify code freq based on NCO command
        codefreq = code_rate - nco_code + (carrfreq - fif)/((frf0 + FCH*fch_step)/code_rate);
        vec_nco_code(loopcnt) = nco_code;
        vec_freq_code(loopcnt) = codefreq;
		olderr_code = err_code;
        oldnco_code = nco_code;
		
        %since it can be time consuming, do a periodic display of current stage
        %also a good point to display debugging information
%         if (rem(loopcnt,10)==0)
%             fprintf('\n er%f er%f nco%f',err_carr,err_freq,nco_carr);
%             fprintf(' %f %f %f %f',vec_freq_code(loopcnt),vec_freq_carr(loopcnt),remcodephase,remcarrphase);
%         end
        
		%absoluteSample(loopcnt) = ftell(fid); %small file
		absoluteSample(loopcnt) = currentSample/floor(sampletype)-remcodephase*(fs/1000)/code_length;
% fprintf(' %f :   %f \n',remcodephase,absoluteSample(loopcnt));
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
        title(['Channel [',num2str(FCH),'] - Prompt Inphase Correlation'])
        subplot(122),plot(p_q)
        grid
        xlabel('milliseconds')
        ylabel('quadrature amplitude')
        title(['Channel [',num2str(FCH),'] - Prompt Quadrature Correlation'])

        for i=1:length(vec_freq_carr)
            if (i<=500),
                vect_local = vec_freq_carr(1:i);
                vec_freq_mean(i) = sum(vect_local)/length(vect_local);
            else
                vect_local = vec_freq_carr(i-499:i);
                vec_freq_mean(i) = sum(vect_local)/length(vect_local);
            end
        end
        
%         figure
%         plot (vec_freq_carr - mycarrfreq,'r'),grid on;
%         hold on
%         plot (vec_freq_mean - mycarrfreq,'g.-');
%         xlabel('milliseconds');
%         ylabel('Frequency [Hz]')
%         title(['Channel [',num2str(FCH),'] - Increment of the local carrier frequency with respect to the nominal one'])
        
        %    figure
        %    plot (vec_errort,'k'),grid on;
        %disp(['Time required to process ',int2str(msec),' msec of data is:',num2str(toc),'sec']);
        
        mean(vet_freq);
        
        figure
        plot(p_i .^2 + p_q .^ 2, 'g.-')
        hold on
        grid
        plot(e_i .^2 + e_q .^ 2, 'bx-');
        plot(l_i .^2 + l_q .^ 2, 'r+-');
        xlabel('milliseconds')
        ylabel('amplitude')
        title(['Channel [',num2str(FCH),'] - Correlation Results'])
        legend('prompt','early','late')
        
%         figure
%         plot(vet_freq,'k.-'), grid on;
%         xlabel('milliseconds')
%         ylabel('Frequency [Hz]')
%         title(['Channel [',num2str(FCH),'] - Local Code Frequency'])
%         mean_freq = mean(vet_freq);
    end
%     mat_file = ['../Data_files/GPS_TrackPLL_PRN',num2str(prn)];
%     save (mat_file);
% else
%     mat_file = ['../Data_files/GPS_TrackPLL_PRN',num2str(prn)];
%     load(mat_file);
end


