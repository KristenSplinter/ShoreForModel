% Filter function
% [Low,High]=fftFilter(y,fs,cutoff,flag);
% y= original time-series
% fs= sampling frequency
% cutoff = cut-off frequency
% flag==1 for a plot, else no plot of the results
% Low = low-passed series
% High = High-passed series
function [Low,High]=fftFilter(y,fs,cutoff,flag)

n=length(y);
dt=1/fs;
bw=1/(n*dt);
time=[1:n]./fs; % -- Mark's original, time is in days from 0 to length
%time=[1:n]./fs+732159; % -- Mel's for 2004-2010 timeseries
%time=[1:n]./fs+721837; % -- Mels' for 1976-2010 timeseries
f=[0:n-1]*bw;

% Low-pass filter
jj=find(f<=cutoff);
m=max(jj);
P=fft(y);
Plow=P;
Plow(m+1:n-m+1)=0;    %Ks change from Plow(m:n-m+1)=0; 
Low=ifft(Plow);
Low=real(Low);

%  High-pass filter
Phigh=P;
Phigh(1:m)=0;     % Ks change from Phigh(1:m-1)=0;
Phigh(n-m+2:n)=0;  % Ks change from Phigh(n-m+2:n)=0;
High=ifft(Phigh);
High=real(High);

% % Add high and low components back together for plotting
% recons = High+Low;

if flag==1
% output data    
Fcutoff=m*bw;
disp('Filter Cut-off Frequency (1/dt)= '); Fcutoff
disp('Filter Cut-off Period (dt)    = '); 1/Fcutoff    
% plotting
figure
subplot(2,1,1)
plot(time,y)
ylabel('y')
title('Original Series')
datetick % ------------------------------------- Mel added
% Plot filtered series
subplot(2,1,2)
plot(time,Low)
datetick % ------------------------------------- Mel added
% Plot high frequency component
hold on
plot(time,High,'r')
title('Low-pass series (blue), High-pass series (red)')
% % Plot both together again
% subplot(3,1,3)
% plot(time,recons,'k')
% title('Reconstructed Time Series (low + high)')
% datetick % ------------------------------------- Mel added
%%%%%
end