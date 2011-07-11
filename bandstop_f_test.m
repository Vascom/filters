%f = [0 f0 f1 f2 f3 1];
N = 256;                %filter length
Clock = 23.75;           %MHz
n2 = N/2+1-1;
n3_start = 1;
n3_stop = 128;
max_1 = 63;
Amp = [1 1 0 0 1 1];    %filter amplitude-frequency characteristic
Wei = [1 200 1];      %weight [0.4 1 0.4]
nn = zeros(95,n2);

%for SB = 0.04:0.04:0.08 %stop band width 0.01:0.01:0.16
    SB = 1/Clock;
    disp(SB*23.75)
    for k = 1:128 %1:83
        %disp(k);
        f0 = k*1/n2 + 1/n2; %k/100
        f1 = f0 + 0.01;
        f2 = f1 + SB;
        f3 = f2 + 0.01;
        f4 = f3 + 0.01;
        Freq = [0 f0 f1 f2 f3 1];
        
        b = firls(N,Freq,Amp,Wei);
        %freqz(b)
        %pause(2);
        n = round(b(1:n2)*8192/max(abs(b(1:n2+1))));
        %hist(abs(n(1:n2-1)),64);
        %pause(1);
        nn(k,:) = n;
        nn2(k,:) = b;
        for p = n3_start:(n3_stop-1) 
            if(abs(n(p)>max_1)) 
                %display(n(p));display(f0);display(SB); 
            end
            %if(abs(n(n2)~=512)) 
            %    disp(n(n2));
            %end
        end
    end
%end