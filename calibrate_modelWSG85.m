function [A, cutoff]= calibrate_modelWSG85(H,T,w,dnum,data,m,dx)

global iCal iWaveCal ihardWire omegaFilterCutOff k omegaRecordMeanCal
global datesCal calStats betaCal tideCal tideSwitchCal contShoreCal
global zShoreCal dnumCal Alin mCal D phi
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
    [AIClin]=akaikeIC(data(iCal),xb(iWaveCal),2);
    rmstrend=sqrt(mean((data(iCal)-xb(iWaveCal)).^2));
    nmstrend = sum((data(iCal)-xb(iWaveCal)).^2)./sum((data(iCal)-mean(data(iCal))).^2);
    dataMean=mean(data);
    rmsErrorMean=sqrt(mean((data(iCal)-dataMean).^2));

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
        %maxYears=ceil(mCal/(365*24)) %KS change 04/05/2012
        %minHours=7*24;%24*dt/(60*60);
        %omegaFilterCutOff=[minHours:24:24*69,24*70:7*24:24*100,124*24:30*24:maxYears*365*24];
        %omegaFilterCutOff=[minHours:24:24*69,24*70:7*24:maxYears*365*24, maxYears*365*24]; %ensures record mean is included
        
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
phitest = [5:5:50 60:10:100 120:20:380];
Dtest = 2.*phitest;
combo = [Dtest(:) phitest(:)];

% Begin loop for trying different lags in omega calculations
for kk=1:length(Dtest)
% Low-pass filter omega time-series with increasing cut-off period to simulate
% slowly changing equilibrium conditions
    %cutoff=1/omegaFilterCutOff(kk);

%    if (strcmp(model,'DMT11')==1)
       % [omegaFiltered, ~]=fftFilter(omega,fs,cutoff,0); % for DMT11 Model
%    else
%        omegaFiltered=sum(omega.^2)./sum(omega); % for DLT10 model for k=1 (equivelent to weighting on H of k=0.5 in DMT11)
%    end
    omegaFiltered = nan(size(omega));
    omegaFiltered(iWaveCal(1)-24*3600./dt.*Dtest(kk)+1:end)=WS85FilterConv(omega(iWaveCal(1)-24*3600./dt.*Dtest(kk)+1:end), Dtest(kk), phitest(kk), dt);
    omegaFiltered(iWaveCal(1)-24*3600./dt.*Dtest(kk)+1:iWaveCal(1)-1)=NaN;
    % Find erosion & accretion times
    ie= omegaFiltered-omega <=0;
    ia= omegaFiltered-omega  >0;

% First pass least squares operating on detrended wave and shoreline data
    X=F.^k.*(omegaFiltered-omega);% NB Using omega without trend initially
    % Kristen add to try and do the filterering WSG85
    keep = find(isnan(X)==0);
    X = X(keep);
    ie = ie(keep); ia = ia(keep);
    mCal2 = length(keep);
    iCal2 = iCal(find(iWaveCal>=keep(1),1,'first'):end);
    iWaveCal2 = iWaveCal(find(iWaveCal>=keep(1),1,'first'):end)-keep(1)+1;
    
    t2 = t(keep)-t(keep(1));
    Xorig=X; % Record orininal detrended X then claculate detrended value with mean preserved
    meanX=mean(X);
    X=detrend(X)+meanX;
    X2=X;
    X2(ie)=0;
    X(ia)=0;
    X=cumtrapz(X)*dt;
    X2=cumtrapz(X2)*dt;
    %B=[ones(mm,1),X(iCal),X2(iCal)];
    B=[ones(mCal2,1),X,X2];
    initialA=B(iWaveCal2,:)\C(iCal2);  % least squares using detrended shoreline series
        
% Compute model data (x)
    % Re-make B-array for least squares with trends back in the wave data
    X=Xorig;
    X2=X;
    X2(ie)=0;
    X(ia)=0;
    X=cumtrapz(X)*dt;
    X2=cumtrapz(X2)*dt;
    %B=[ones(mm,1),X(iCal),X2(iCal)];
    B=[ones(mCal2,1),X,X2];
    x=B*initialA; % Model with cross-shore trend only
    clear B
    %px = polyfit(t(iCal),x,1);waveTrendArray(kk)=px(1); waveOffsetArray(kk)=px(2);% Compute and record wave trend and offset
    px = polyfit(t2,x,1);waveTrendArray(kk)=px(1); waveOffsetArray(kk)=px(2);% Compute and record wave trend and offset
  
% Second regression here to evaluate non cross-shore trend in the data
    B=[ones(mCal2,1),t2,x];
    regressCoeff=B(iWaveCal2,:)\data(iCal2);
    clear B
        
% Adjust coefficients and compute final modeled shoreline
    A(1)=(initialA(1).*regressCoeff(3))+regressCoeff(1); % Aggregate offsets
    A(2)=regressCoeff(2); % This is the isolated non-crosshore component of trend
    A(3)=initialA(2).*regressCoeff(3); % C-
    A(4)=initialA(3).*regressCoeff(3); % C+
    A=A';
    %B=[ones(mm,1),t(iCal),X(iCal),X2(iCal)]; % Full model including all four parameters - offset, trend and C+/C-
    B=[ones(mCal2,1),t2,X,X2]; % Full model including all four parameters - offset, trend and C+/C-
   x=B*A;  % Full model including cross-shore plus other trends (e.g. longshore transport)
    %tmp1 = regress(B(:,3).*A(3),[ones(length(iCal),1),t(iCal)]);
    %tmp2 = regress(B(:,4).*A(4),[ones(length(iCal),1),t(iCal)]);

    clear B
    

% Compute correlation coefficient, rms error & BSS
    c=corrcoef(data(iCal2),x(iWaveCal2));cc(kk)=c(1,2);
    clc; disp(['Cut-off period for omega (days) D = ' num2str(Dtest(kk)) ', phi = ' num2str(phitest(kk)), ' Correlation coefficient = ' num2str(cc(kk))]);
    rmsError(kk)=sqrt(mean((data(iCal2)-x(iWaveCal2)).^2));
    format SHORT ENG
    [BSS(kk)]=brierss(data(iCal2),x(iWaveCal2),xb(iWaveCal2+keep(1)-1),dx);
    [AIC(kk)]=akaikeIC(data(iCal2),x(iWaveCal2),5);
    nmsError(kk) = sum((data(iCal2)-x(iWaveCal2)).^2)./sum((data(iCal2)-mean(data(iCal2))).^2);
    Aarray(:,kk)=A; clear A x
    

end % end major loop
% find the maximum correlation index and recompute the model for best
% values
    jj=max(find(cc==max(cc)));
    cutoff = [Dtest(jj) phitest(jj)]
    D = Dtest(jj);
    phi = phitest(jj);
    %cutoff=1/omegaFilterCutOff(jj);
    %[omegaFiltered, ~]=fftFilter(omega,fs,cutoff,0);
    %omegaFiltered=WS85Filter(omega, D, phi, dt);
    omegaFiltered = nan(size(omega));
    omegaFiltered(iWaveCal(1)-24*3600./dt.*D+1:end)=WS85FilterConv(omega(iWaveCal(1)-24*3600./dt.*D+1:end), D, phi, dt);
    omegaFiltered(iWaveCal(1)-24*3600./dt.*Dtest(kk)+1:iWaveCal(1)-1)=NaN;

    ie= omegaFiltered-omega <=0;
    ia= omegaFiltered-omega  >0;

    % least squares
    X=F.^k.*(omegaFiltered-omega);%
    %KS
    keep = find(isnan(X)==0);
    X = X(keep);
    ie = ie(keep); ia = ia(keep);
    mCal2 = length(keep);
    iCal2 = iCal(find(iWaveCal>=keep(1),1,'first'):end);
    iWaveCal2 = iWaveCal(find(iWaveCal>=keep(1),1,'first'):end)-keep(1)+1;
    dnumCal2 = dnumCal(keep);
    t2 = t(keep)-t(keep(1));
    %end KS
    X2=X;
    X2(ie)=0;
    X(ia)=0;
    X=cumtrapz(X)*dt;
    X2=cumtrapz(X2)*dt;
    A=Aarray(:,jj);    % best least squares value
        
% Compute model data (x) 
    B=[ones(length(t2),1),t2,X,X2];
    x=B*A;
        
% Compute wave trend line
    xbWave=t2*waveTrendArray(jj)+waveOffsetArray(jj);
    
% Compute correlation coefficient, rms error & BSS
    ccFinal=cc(jj);
    rmsErrorFinal=rmsError(jj);
    BSSFinal=BSS(jj);
    AkaikeFinal = AIC(jj);
    nmsErrorFinal = nmsError(jj);

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
    plot(datesCal(iCal2),data(iCal2),'Color',[0.5 0.5 0.5],'LineWidth',2)
    hold on
    %plot(datesCal(iCal),x,'k','LineWidth',2)
    plot(dnumCal2(iWaveCal2),x(iWaveCal2),'k','LineWidth',2)
    %plot(dnum,xb,'k')
    datetick('x','yyyy')
    xlabel('Time (years)')
    ylabel('\Delta shoreline position (m)')
    grid on
    legend('Data','Model')

if (ihardWire==1); 
figure 
subplot(2,2,1)% BSS Plot
    plot(phitest, BSS,'k.--')
    hold on
    xlabel('phi  (days)')
    ylabel('Brier Skills Score')
    grid on
    legend('BSS')
    plot([0 phitest(end)], [0.3 0.3],'k--')
    plot([0 phitest(end)], [0.6 0.6],'k--')
    plot([0 phitest(end)], [0.8 0.8],'k--')
    
    subplot(2,2,2) % Akaike Info criteria Plot
    plot(phitest,AIC-AIClin,'k','LineWidth',2)
    hold on
    xlabel('phi (days)')
    ylabel('Relative AIC')
    grid on
    plot([0 phitest(end)], [-1 -1],'k--','LineWidth',2)
    %legend('relative AIC')


subplot(2,2,3) % Correlation Coefficient Plot
    plot(phitest,cc.^2,'k','LineWidth',2)
    hold on
    xlabel('phi (days)')
    ylabel('Correlation squared')
    grid on
    plot([0 phitest(end)], [cr(1,2).^2 cr(1,2).^2],'k--','LineWidth',2)
    legend('R^2','Linear trend R^2')

subplot(2,2,4) % NMS error plot
    plot(phitest, nmsError,'k-.')
    hold on
    xlabel('phi (days)')
    ylabel('NMSE (m)')
    grid on
    plot([0 phitest(end)], [nmstrend nmstrend],'k--')
    legend('NMSE','Linear trend NMSE')

% % Plot Filtered omega series
% subplot(2,2,4)
%     plot(dnum,omegaFiltered,'k')
%     hold on
%     plot([min(dnum) max(dnum)],[omegaRecordMean omegaRecordMean],'k--')
%     datetick('x','yyyy')
%     xlabel('Time (years)')
%     ylabel('Filtered Omega Series')
%     grid on
%     legend('Filtered Omega','Record Mean Value')
end    

calStats.datesCal=datesCal;
calStats.Nyr = (datesCal(end)-datesCal(1))/365;
calStats.A=A;
calStats.waveTrend = waveTrendArray(jj);
calStats.xbWave = xbWave;
calStats.cutoff=cutoff;
calStats.phi = phi;
calStats.D = D;
calStats.cc=cc(jj);
calStats.rmse=rmsError(jj);
calStats.nmse = nmsError(jj);
calStats.BSS=BSS(jj);
calStats.AIC = AIC(jj);
calStats.omegaRecordMean=omegaRecordMean;
calStats.varOmega = var(omega);
calStats.varShore = var(data(iCal));
calStats.dx=dx;
calStats.cr=cr(1,2);
calStats.rmstrend=rmstrend;
calStats.nmstrend=nmstrend;
calStats.AIClin = AIClin;
calStats.omegaFilter=omegaFiltered;
calStats.omega=omega;
calStats.F=F;
calStats.k=k;
calStats.dt=dt;
calStats.linModelA=Alin;
calStats.contourShore=contShoreCal;
