% Acquisizione FFT in Time domain

function [doppler, code, st] = signal_acquisition_FFT(file,prn)
global f_sampling;
global nominalfreq;
global isNewRun;
global sampletype;
global code_rate;
global code_length;

if (isNewRun==1)

    acq_metric = 1.5;

    seek_sec = 30;
    
    %% Extract data
% if (sampletype == 1)
%     disp('file: 8 bit real samples S0,S1,S2,...\n');
% else
%     disp('file: 8 bit complex samples I0,Q0,I1,Q1,I2,Q2...');
% end
n_code_extr = 11; % extract 11 times of codelength
fseek(file,ceil(sampletype*f_sampling*seek_sec),-1); % skip samples
% [gpsdata,scount] = fread(file,n_code_extr*sampletype*f_sampling*code_length/code_rate,'schar'); % datatype = signed char, little endian
[gpsdata,scount] = fread(file,n_code_extr*sampletype*f_sampling*code_length/code_rate,'int8'); % datatype = int16, little endian
if (sampletype==2)
    data1=gpsdata(1:2:end);    
    data2=gpsdata(2:2:end);    
    gpsdata=data1 + 1i*data2;
end 
%%
%     fseek(file,ceil(f_sampling*seek_sec),-1);
%     [gpsdata,scount] = fread(file,40000*6,'schar'); %...solitamente leggere schar, da NordNav
    %gpsdata = gpsdata-128;  % ATT! minus 127 required when the data set has been acquired with the NavSAS front end
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %                        BLANKING
    % index = find(abs(gpsdata)>2);
    % gpsdata(index)=0;
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ts = 1/f_sampling;
    code_rate = 1.023*1e6;              % Nominal GPS C/A code rate.
    num_samples = f_sampling/code_rate; % Number of samples per chip.


    K=1; % Downsampling

    Fs=f_sampling/K;
    Rc=1.023e6;
    CodeLen=1023;
    fif=rem(nominalfreq,Fs/K);

    gpsdata=gpsdata(1:K:end);

    N=floor(Fs*CodeLen/Rc)+1; %lay 1ms samples - du 1 chu ky cua CA code
    Dopplerstep=250;

    C=[];

    S=prn;

    Loc = (double(gnssCACode(S,"NavIC L5-SPS")-0.5))';

    Loc = [Loc Loc(1)];  % extend the local code of 1 chip. The floor operation used with high sampling freq lead to sampling the local code at chip 1024.

    idx=1;

    for FD=-8000:Dopplerstep:8000;

        corr=zeros(1,N)+j*zeros(1,N);

        % Ricampiono

        k=0:N-1;
        %Lay mau local CA (1023 chips, 1 chips chua nhieu mau vi tan so lay
        %mau lon hon chip rate
        SigLOC=Loc(floor(k*Rc/Fs)+1);

        % FFT codice locale e complesso coniugato

        SigLOCFFT=conj(fft(SigLOC,N));

        argx=2*pi*(fif+FD)/Fs;
        carrI=cos(argx*k);
        carrQ=sin(argx*k);

        for M=0:5
            SigIN=gpsdata(N*M+1:N*M+N)';

            % Demodulo

%             I=SigIN.*carrI-SigIN.*carrQ;
%             Q=SigIN.*carrQ+SigIN.*carrQ;

            SigINIQ=SigIN.*(carrI+j*carrQ);

            corr=corr+abs(ifft(fft(SigINIQ,N).*(SigLOCFFT)));

        end

        C(idx,:)=corr;
        idx=idx+1;

    end

    % --- Find the main peak in the correlation floor and the corresponding frequency bin index
    [bb ind_mixf] = max(max(C'));
    [bb ind_mixc] = max(max(C));

    if (ind_mixc < ceil(num_samples/K)),
        vect_search_peak2 = [zeros(1,2*ceil(num_samples/K)), C(ind_mixf,(2*ceil(num_samples/K)):end)];
    elseif (ind_mixc < ceil(num_samples/K))
        vect_search_peak2 = [C(ind_mixf,1:(end-2*ceil(num_samples/K)):end), zeros(1,2*ceil(num_samples/K))];
    else
        vect_search_peak2 = [C(ind_mixf,1:(ind_mixc-ceil(num_samples/K))),zeros(1,2*ceil(num_samples/K)-1),C(ind_mixf,(ind_mixc+ceil(num_samples/K)):end)];
    end


    % --- Find the second highest peak in the correlation floor
    second_peak = max(vect_search_peak2);


    % --- compare the acquisition metric to a predefined threshold

    if ((bb/second_peak) > acq_metric)
        figure, surf(C), shading interp, title(prn);
        fprintf('%d ...acquired satellite PRN %i\n ',(bb/second_peak),prn)
        code = ceil(f_sampling*seek_sec)+((ind_mixc-1)*K);% - 20;
        % 20 is a sistematic error, which simulates the number of offset samples, which makes the acquisition estimation not correct, since it took time (not negligible).
        % The effect on the carrier frequency can be neglected, no problems for that!

        doppler = (nominalfreq-8e3) + (ind_mixf-1)*Dopplerstep;
        st = 1;
    else
        fprintf('...no satellite ')
        code = 0;
        doppler = 0;
        st = 0;
    end
    mat_file = ['../Data_files/GPS_Acq_PRN',num2str(prn)];
    save(mat_file);
else
    mat_file = ['../Data_files/GPS_Acq_PRN',num2str(prn)];
    load(mat_file);
end
