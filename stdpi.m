S4=[]

for k=60000:1000:365000
    NBP=[]
    WBP=[]
    for i=(k-59999):20:k
        NBP= [NBP sum(p_i(i:i+19))^2+sum(p_q(i:i+19))^2];
        WBP= [WBP sum(p_i(i:i+19).^2)+sum(p_q(i:i+19).^2)];
    end
    SIraw=NBP-WBP;
    SItrend=mean(SIraw);
    SI=SIraw./SItrend;

    S4= [ S4 sqrt((mean(SI.^2)-mean(SI)^2)/mean(SI)^2)];
%     S4= [ S4 1];
end
%         fid = fopen('s4.txt','a');
%         fprintf(fid,'%d %f \n',prn,S4);
%         fclose(fid);
%         
     