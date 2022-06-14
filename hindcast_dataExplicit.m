function hindcast_dataExplicit(cutoff,A,H,T,w,dnum,data)

global iShore iWave mm m  k  dates_survey iWave omegaRecordMean hindStats dx
global omegaFiltered  contShore hindDir

% Initialise model data array & omegaFiltered array
    x=zeros(size(iShore));

% Model/Data timestep - seconds - NB must be a regular time-step in data!
    dt=round((dnum(2)-dnum(1))*24*60*60);
    % Linear model (for comparison)
    t=[0:m-1].*dt;t=t';
    C=data;
    B2=[ones(mm,1),t(iWave)];
    A2=B2\C(iShore)
    xb=B2*A2;

% Bench mark values using linear trend
    cr=corrcoef(data(iShore),xb);
    rmstrend=sqrt(mean((data(iShore)-xb).^2));
    dataMean=mean(data);
    rmsErrorMean=sqrt(mean((data(iShore)-dataMean).^2));

    % compute raw omega series
    omega=H./(w*T); 
    % Equilibrium omega value for the whole record (info only)
    omegaRecordMean=mean(omega);%

% if (strcmp(model,'DLT10')==1)
%     F =  omega.^2; 
%     omegaFiltered=sum(omega.^2)./sum(omega); % for DLT10 model for k=1 (equivelent to weighting on H of k=0.5 in DMT11)
%  
% elseif (strcmp(model,'DMT11')==1)
    F = calcP(H,T); %(1025*9.81.^2/(64*pi))*H.^2 .* T; % 
    fs=(60*60)/dt;
    [omegaFiltered, high]=fftFilter(omega,fs,cutoff,0);
% end
% if strcmp(tideSwitch,'y')==1
%    RHigh = calcR2(beta,H,T)+tide
%    zeroF = RHigh<contShore;
%    F(zeroF)=0;
% end
    ie= omegaFiltered-omega <=0;
    ia= omegaFiltered-omega  >0; 

% least squares
%     if iConstC ~= 1
        X=F.^k.*(omegaFiltered-omega);% 
        X2=X;
        X2(ie)=0;
        X(ia)=0;
        X=cumtrapz(X)*dt;
        X2=cumtrapz(X2)*dt;
%     else
%         X=F.^k.*(omegaFiltered-omega)*dt;%
%         X=cumtrapz(X)*dt;
%         X2=t;
%     end
        B=[ones(mm,1),t(iWave),X(iWave),X2(iWave)];
%        x=B*A; clear B
if strcmp(hindDir,'forward')
    xs = data(1)+B(:,3:4)*A(3:4);   
    x = xs + B(:,2)*A(2);
elseif strcmp(hindDir,'backward')
    Bb = [ones(mm,1),t(iWave)-t(end),X(iWave)-X(end),X2(iWave)-X2(end)];
    xs = data(end)+Bb(:,3:4)*A(3:4);   
    x = xs + Bb(:,2)*A(2);
end

    c=corrcoef(data(iShore),x); cc=c(1,2);
    rmsError=sqrt(mean((data(iShore)-x).^2));
    BSS=brierss(data(iShore),x,xb,dx);
   cs=corrcoef(data(iShore),xs); ccs=cs(1,2);
    rmsErrors=sqrt(mean((data(iShore)-xs).^2));
    BSSs=brierss(data(iShore),xs,xb,dx);
    
    % Plot the best results
    figure
    subplot(2,1,1)
    plot(dnum,omega,'Color',[0 0 0],'LineWidth',1)
    datetick('x','yyyy')
    ylabel('\Omega')
    grid on
    title('Hindcast Results')
    subplot(2,1,2)
    plot(dates_survey(iShore),data(iShore),'Color',[0.5 0.5 0.5],'LineWidth',2)
    hold on
    plot(dates_survey(iShore),xs,'k--','LineWidth',2)
    plot(dates_survey(iShore),x,'k','LineWidth',2)
    plot(dates_survey(iShore),xb,'r','LineWidth',2)
    datetick('x','yyyy')
    xlabel('Time (years)')
    ylabel('Shoreline position (m)')
    grid on
    legend('Data','Model no trend','Model with Trend','linear Model')

hindStats.dates=dates_survey;
hindStats.ccFull=cc;
hindStats.rmseFull=rmsError;
hindStats.BSSFull=BSS;
hindStats.ccShort=ccs;
hindStats.rmseShort=rmsErrors;
hindStats.BSSShort=BSSs;

hindStats.omegaRecordMean=omegaRecordMean;
hindStats.dx=dx;
hindStats.cr=cr(1,2);
hindStats.rmstrend=rmstrend;
hindStats.omegaFilter=omegaFiltered;
hindStats.F=F;
hindStats.k=k;
hindStats.dt=dt;
hindStats.omega=omega;
hindStats.x=x;
hindStats.bTrend=A2;
hindStats.contourShore=contShore;
hindStats.linModel = xb;

