
% Control PLL.
% This program control the PLL status and estimates the C/No ratio. 


function [PLL_status] = control_PLL_v1(act_IP,act_QP,loop,thres);

global FIFO_IP;
global FIFO_QP;
global f_sampling;
global samplesPDI;
global CNo_Emanuela;
global CNo_bluebook;
global CNo_SNV;
global CNo_MM;
global CNo_Bea;

%%%%%EF Beq = 2*(f_sampling/samplesPDI);  % real signal component + complex noise estimation
Beq = (f_sampling/samplesPDI);  % real signal component + complex noise estimation

% Update the FIFO structure used to estimate the C/No

FIFO_IP = [act_IP, FIFO_IP(1:end-1)];
FIFO_QP = [act_QP, FIFO_QP(1:end-1)];

PLL_status = 0;
index = 0;

if (rem(loop*10,1000)==0)
    
    index = index + 1;
    %*** Estimate the C/No ratio - Emanuela's Method ***
    %disp('Here estimate the C/No and check the lock of the tracking loops')
    
    Pn = 2*var(FIFO_QP);    
    Ptot = var(FIFO_IP + i*FIFO_QP);
    Ps = Ptot-Pn;
    if Ps > 0
        SNR_RSCN = Ps/Pn;  
    else
        SNR_RSCN = 1;     % in case of Pn<0, it means that a phase rotation occured.  The estimator is not able to provide a correct value of SNR, 
                          % which is set equal to 1 in order to have 0 in the C/No.
    end
    
    CNo = 10*log10(SNR_RSCN * Beq);
    CNo_Emanuela = [CNo_Emanuela, CNo]; 
    % ***************************************************
    %*** Estimate the C/No ratio - Blue Book Method ***
    %disp('Here estimate the C/No and check the lock of the tracking loops')
    MSamplesPerBit = 2;
    
    K = length(FIFO_IP)/MSamplesPerBit; % represents the number of bits used in the estimation
    HTotSampleNumber = K*MSamplesPerBit;
    
    I_sam = FIFO_IP(1:HTotSampleNumber);
    Q_sam = FIFO_QP(1:HTotSampleNumber);
    
    CMatrix = reshape(I_sam,MSamplesPerBit,K) + j*reshape(Q_sam,MSamplesPerBit,K);
    
    WBP = sum(abs(CMatrix).^2);
    
    NBP = (sum(real(CMatrix))).^2 + (sum(imag(CMatrix))).^2;
    
    NP = NBP./WBP;
    
    mNP = mean(NP);
    
    CNo_PRM = (mNP-1)/(MSamplesPerBit-mNP)/10e-3;
    
    if CNo_PRM < 0
        CNo_PRM = NaN;
    end
    
    CNo = 10*log10(CNo_PRM);
    CNo_bluebook = [CNo_bluebook, CNo]; 
    % **************************************************  
    %*** Estimate the C/No ratio - SNV Method ***
    
    %%%%%EF Ps = mean(abs(FIFO_IP))^2 + mean(abs(FIFO_QP))^2;
    Ps = mean(abs(FIFO_IP))^2;
    Ptot = mean(abs(FIFO_IP+i*FIFO_QP).^2)
    Pn = Ptot-Ps;
    
    if Pn > 0
        SNR_SNV = Ps/Pn;
    else
        SNR_SNV = NaN;
    end
    
    CNo = 10*log10(SNR_SNV * Beq);
    CNo_SNV = [CNo_SNV, CNo];
    
    % ***************************************************
    %*** Estimate the C/No ratio - Moments Method ***
    
    m2 = mean(abs(FIFO_IP+i*FIFO_QP).^2);
    m4 = mean(abs(FIFO_IP+i*FIFO_QP).^4);    
    Ps = sqrt(2*m2^2 - m4);
    Pn = m2-Ps;    
    if Pn > 0
        SNR_MM = Ps/Pn;
    else
        SNR_MM = NaN;
    end
    
    CNo = 10*log10(SNR_MM * Beq)
    CNo_MM = [CNo_MM, CNo];
    
    % **************************************************  
    %*** Estimate the C/No ratio - Beaulieu ***
    
    %%%%%EF  
    
    Nsym = length(FIFO_IP)-1;
    X_re = FIFO_IP(1:end-1);
    Y_im = FIFO_IP(2:end);
    func2 = sum((abs(X_re) - abs(Y_im)).^2 ./ (X_re.^2 + Y_im.^2));
    if func2 > 0
        SNR_Bea = Nsym/func2/2;
    else
        SNR_Bea = NaN;
    end
    CNo = 10*log10(SNR_Bea * Beq)
    CNo_Bea = [CNo_Bea, CNo];
    % ************************************************** 

end

return