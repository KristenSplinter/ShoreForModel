function [A, cutoff]= calibrate_model(H,T,w,dnum,data,m,dx)

global iCal iWaveCal ihardWire omegaFilterCutOff k omegaRecordMeanCal
global datesCal calStats betaCal tideCal tideSwitchCal contShoreCal
global zShoreCal dnumCal Alin mCal
mm=length(iCal);


% Initialise model data array & omegaFiltered array
    x=zeros(size(dnum));    %iCal
    omegaFiltered=zeros(size(dnum))';   %iCal

% Model/Data timestep - seconds - NB must be a regular time-step in data!
    dt=round((dnum(2)-dnum(1))*24*60*60);
        fs=(60*60)/dt;
        
% Linear model (for comparison)
    t=[0:mCal-1].*dt;t=t';
    C=data;
    %B2=[ones(mm,1),t(iWaveCal)];
    B2=[ones(mCal,1),t];
    Alin=B2(iWaveCal,:)\C(iCal);
    xb=B2*Alin;

    % Bench mark values using linear trend
    cr=corrcoef(data(iCal),xb(iWaveCal));
    rmstrend=sqrt(mean((data(iCal)-xb(iWaveCal)).^2));
    nmstrend = sum((data(iCal)-xb(iWaveCal)).^2)./sum((data(iCal)-mean(data(iCal))).^2);
    dataMean=mean(data);
    rmsErrorMean=sqrt(mean((data(iCal)-dataMean).^2));
    [AIClin]=akaikeIC(data(iCal),xb(iWaveCal),2);

    % compute raw omega series

    omega=H./(w.*T); 
    % Equilibrium omega value for the whole record (info only)
    omegaRecordMean=mean(omega);%
    omegaRecordMeanCal=omegaRecordMean;
    
% if (strcmp(model,'DLT10')==1)
%     F =  omega.^2; 
%     omegaFilterCutOff=NaN; % i.e. no filtering
% elseif (strcmp(model,'DMT11')==1)
    F = calcP(H,T); %(1025*9.81.^2/(64*pi))*H.^2 .* T; % Offshore Wave Power
%    [omegaStorm, high]=fftFilter(omega,fs,2*24,0); % for DMT11 Model

    %storm response term
%     omegaStorm=WS85Filter(omega, 30, 5, dt);
%     %now need to cutoff all the data that is nan because filter is looking
%     %backwards.
%     iCalFull=iCal;  %keep original counter
%     st = 30*24+1;
%     H = H(st+1:end);
%     T = T(st+1:end);
%     tideCal = tideCal(st+1:end);
%     F = F(st+1:end);
    % define omegaFiltered low pass cutt-off in hours
    
    if ihardWire==0
        omegaFilterCutOff=omegaFilterCutOff*24;            
    else
        maxYears=ceil(mCal/(365*24)) %KS change 04/05/2012
        minHours=7*24;%24*dt/(60*60);
        %omegaFilterCutOff=[minHours:24:24*69,24*70:7*24:24*100,124*24:30*24:maxYears*365*24];
        omegaFilterCutOff=[minHours:24:24*69,24*70:7*24:maxYears*365*24, maxYears*365*24]; %ensures record mean is included
        
        %omegaFilterCutOff=[minHours:23,24:24:24*69,24*70:7*24:24*100,124*24:30*24:maxYears*365*24];
        %omegaFilterCutOff=[minHours:24:maxYears*365*24];
    end
% end
% if strcmp(tideSwitchCal,'y')==1
%    RHigh = calcR2(betaCal,H,T)+tideCal;
%    zeroF = RHigh<contShoreCal;
%     F(zeroF)=0;
% end

% Define detrended shoreline data for optimisation
C=detrend(data);

% Begin loop for trying different lags in omega calculations
for kk=1:length(omegaFilterCutOff)
% Low-pass filter omega time-series with increasing cut-off period to simulate
% slowly changing equilibrium conditions
    cutoff=1/omegaFilterCutOff(kk);

%    if (strcmp(model,'DMT11')==1)
        [omegaFiltered, high]=fftFilter(omega,fs,cutoff,0); % for DMT11 Model
%    else
%        omegaFiltered=sum(omega.^2)./sum(omega); % for DLT10 model for k=1 (equivelent to weighting on H of k=0.5 in DMT11)
%    end

% Find erosion & accretion times
    ie= omegaFiltered-omega <=0;
    ia= omegaFiltered-omega  >0;

% First pass least squares operating on detrended wave and shoreline data
    X=F.^k.*(omegaFiltered-omega);% NB Using omega without trend initially
    Xorig=X; % Record orininal detrended X then claculate detrended value with mean preserved
    meanX=mean(X);
    X=detrend(X)+meanX;
    X2=X;
    X2(ie)=0;
    X(ia)=0;
    X=cumtrapz(X)*dt;
    X2=cumtrapz(X2)*dt;
    %B=[ones(mm,1),X(iCal),X2(iCal)];
    B=[ones(mCal,1),X,X2];
    initialA=B(iWaveCal,:)\C(iCal);  % least squares using detrended shoreline series
        
% Compute model data (x)
    % Re-make B-array for least squares with trends back in the wave data
    X=Xorig;
    X2=X;
    X2(ie)=0;
    X(ia)=0;
    X=cumtrapz(X)*dt;
    X2=cumtrapz(X2)*dt;
    %B=[ones(mm,1),X(iCal),X2(iCal)];
    B=[ones(mCal,1),X,X2];
    x=B*initialA; % Model with cross-shore trend only
    clear B
    %px = polyfit(t(iCal),x,1);waveTrendArray(kk)=px(1); waveOffsetArray(kk)=px(2);% Compute and record wave trend and offset
    px = polyfit(t,x,1);waveTrendArray(kk)=px(1); waveOffsetArray(kk)=px(2);% Compute and record wave trend and offset
  
% Second regression here to evaluate non cross-shore trend in the data
    B=[ones(mCal,1),t,x];
    regressCoeff=B(iWaveCal,:)\data(iCal);
    clear B
        
% Adjust coefficients and compute final modeled shoreline
    A(1)=(initialA(1).*regressCoeff(3))+regressCoeff(1); % Aggregate offsets
    A(2)=regressCoeff(2); % This is the isolated non-crosshore component of trend
    A(3)=initialA(2).*regressCoeff(3); % C-
    A(4)=initialA(3).*regressCoeff(3); % C+
    A=A';
    %B=[ones(mm,1),t(iCal),X(iCal),X2(iCal)]; % Full model including all four parameters - offset, trend and C+/C-
    B=[ones(mCal,1),t,X,X2]; % Full model including all four parameters - offset, trend and C+/C-
   x=B*A;  % Full model including cross-shore plus other trends (e.g. longshore transport)
    %tmp1 = regress(B(:,3).*A(3),[ones(length(iCal),1),t(iCal)]);
    %tmp2 = regress(B(:,4).*A(4),[ones(length(iCal),1),t(iCal)]);

    clear B
    

% Compute correlation coefficient, rms error & BSS
    c=corrcoef(data(iCal),x(iWaveCal));cc(kk)=c(1,2);
    clc; disp(['Cut-off period for omega (days) ' num2str(omegaFilterCutOff(kk)./24) ',' 'Correlation coefficient = ' num2str(cc(kk))]);
    rmsError(kk)=sqrt(mean((data(iCal)-x(iWaveCal)).^2));
    nmsError(kk) = sum((data(iCal)-x).^2)./sum((data(iCal)-mean(data(iCal))).^2);

    format SHORT ENG
    [BSS(kk)]=brierss(data(iCal),x(iWaveCal),xb(iWaveCal),dx);
    AIC(kk) = akaikeIC(data(iCal),x(iWaveCal),5);
    Aarray(:,kk)=A; clear A x
    

end % end major loop
% find the maximum correlation index and recompute the model for best
% values
    jj=max(find(cc==max(cc)));
    
    cutoff=1/omegaFilterCutOff(jj);
    [omegaFiltered, high]=fftFilter(omega,fs,cutoff,0);
    ie=find(omegaFiltered-omega <=0);
    ia=find(omegaFiltered-omega  >0); 

    % least squares
    X=F.^k.*(omegaFiltered-omega);%
    X2=X;
    X2(ie)=0;
    X(ia)=0;
    X=cumtrapz(X)*dt;
    X2=cumtrapz(X2)*dt;
    A=Aarray(:,jj);    % best least squares value
        
% Compute model data (x) 
    B=[ones(length(t),1),t,X,X2];
    x=B*A;
        
% Compute wave trend line
    xbWave=t*waveTrendArray(jj)+waveOffsetArray(jj);
    
% Compute correlation coefficient, rms error & BSS
    ccFinal=cc(jj);
    rmsErrorFinal=rmsError(jj);
    BSSFinal=BSS(jj);
    nmsErrorFinal = nmsError(jj);
    AICFinal = AIC(jj);

    % Plot the best results
    figure
    clf
    subplot(2,1,1)
    plot(dnum,omega,'Color',[0 0 0],'LineWidth',1)
    datetick('x','yyyy')
    ylabel('\Omega')
    grid on
    title('Calibration Results')
    subplot(2,1,2)
    plot(datesCal(iCal),data(iCal),'Color',[0.5 0.5 0.5],'LineWidth',2)
    hold on
    %plot(datesCal(iCal),x,'k','LineWidth',2)
    plot(dnumCal,x,'k','LineWidth',2)
    %plot(dnum,xb,'k')
    datetick('x','yyyy')
    xlabel('Time (years)')
    ylabel('\Delta shoreline position (m)')
    grid on
    legend('Data','Model')

if (ihardWire==1); 
figure 
subplot(2,2,1)% BSS Plot
    semilogx(omegaFilterCutOff./24,BSS,'k.--')
    hold on
    xlabel('Filter cut-off  (days)')
    ylabel('Brier Skills Score')
    grid on
    legend('BSS')
    plot([.01 2000], [0.3 0.3],'k--')
    plot([.01 2000], [0.6 0.6],'k--')
    plot([.01 2000], [0.8 0.8],'k--')

subplot(2,2,2) % Correlation Coefficient Plot
    semilogx(omegaFilterCutOff./24,cc,'k','LineWidth',2)
    hold on
    xlabel('Filter cut-off  (days)')
    ylabel('Correlation')
    grid on
    plot([.01 2000], [cr(1,2) cr(1,2)],'k--','LineWidth',2)
    legend('Correlation coefficient','Linear trend correlation coefficient')

subplot(2,2,3) % RMS error plot
    semilogx(omegaFilterCutOff./24,rmsError,'k-.')
    hold on
    xlabel('Filter cut-off (days)')
    ylabel('RMSE (m)')
    grid on
    plot([.01 2000], [rmstrend rmstrend],'k--')
    legend('RMSE','Linear trend RMSE')

% Plot Filtered omega series
subplot(2,2,4)
    plot(dnum,omegaFiltered,'k')
    hold on
    plot([min(dnum) max(dnum)],[omegaRecordMean omegaRecordMean],'k--')
    datetick('x','yyyy')
    xlabel('Time (years)')
    ylabel('Filtered Omega Series')
    grid on
    legend('Filtered Omega','Record Mean Value')
end    
calStats.datesCal=datesCal;
calStats.A=A;
calStats.waveTrend = waveTrendArray(jj);
calStats.xbWave = xbWave;
calStats.cutoff=cutoff;
calStats.cc=cc(jj);
calStats.rmse=rmsError(jj);
calStats.nmse=nmsError(jj)
calStats.BSS=BSS(jj);
calStats.AIC = AIC(jj)
calStats.omegaRecordMean=omegaRecordMean;
calStats.varOmega = var(omega);
calStats.varShore = var(data(iCal));
calStats.dx=dx;
calStats.cr=cr(1,2);
calStats.rmstrend=rmstrend;
calStats.nmstrend=nmstrend;
calStats.AIClin=AIClin;
calStats.omegaFilter=omegaFiltered;
calStats.omega=omega;
calStats.F=F;
calStats.k=k;
calStats.dt=dt
calStats.linModelA=Alin;
calStats.contourShore=contShoreCal;
