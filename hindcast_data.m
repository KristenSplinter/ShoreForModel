function hindcast_data(cutoff,A,H,T,w,dnum,data)

global iShore iWave mm m model k iConstC iLinTerm dates_survey iWave omegaRecordMean hindStats dx
global omegaFiltered beta tide tideSwitch contShore Alin

D = cutoff(1);
phi = cutoff(2);
% Initialise model data array & omegaFiltered array
    x=zeros(size(iShore));

% Model/Data timestep - seconds - NB must be a regular time-step in data!
    dt=round((dnum(2)-dnum(1))*24*60*60);
    % Linear model (for comparison)
    t=[0:m-1].*dt;t=t';
    C=data;
    B2=[ones(m,1),t];
    %A2=B2\C(iShore)
    % loop to solve for best fit offset in hindcast data set
    Atest = [min(data):1:max(data)]
    for ii=1:length(Atest)
    xbt(ii)=0;
    end
    xb=B2*Alin;
    xbr = data-xb(iWave);
    Atest = [ones(size(xbr),1)]\xbr
    Alin(1)=Alin(1)+Atest; clear xb
   xb=B2*Alin;
% Bench mark values using linear trend
    cr=corrcoef(data(iShore),xb(iWave));
    [AIClin]=akaikeIC(data(iShore),xb(iWave),2);
    rmstrend=sqrt(mean((data(iShore)-xb(iWave)).^2));
    nmstrend = sum((data(iShore)-xb(iWave)).^2)./sum((data(iShore)-mean(data(iShore))).^2);
    dataMean=mean(data);
    rmsErrorMean=sqrt(mean((data(iShore)-dataMean).^2));

    % compute raw omega series
    omega=H(:)./(w*T(:)); 
    % Equilibrium omega value for the whole record (info only)
    omegaRecordMean=mean(omega);%

% if (strcmp(model,'DLT10')==1)
%     F =  omega.^2; 
%     omegaFiltered=sum(omega.^2)./sum(omega); % for DLT10 model for k=1 (equivelent to weighting on H of k=0.5 in DMT11)
%  
% elseif (strcmp(model,'DMT11')==1)
    F = calcP(H,T); %(1025*9.81.^2/(64*pi))*H.^2 .* T; % 
    fs=(60*60)/dt;
    %omegaFiltered=WS85FilterConv(omega, D, phi, dt);
    omegaFiltered = nan(size(omega));
    omegaFiltered(iWave(1)-24*3600./dt.*D+1:end)=WS85FilterConv(omega(iWave(1)-24*3600./dt.*D+1:end), D, phi, dt);
    omegaFiltered(iWave(1)-24*3600./dt.*D+1:iWave(1)-1)=NaN;

    %[omegaFiltered, ~]=fftFilter(omega,fs,cutoff,0);
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

    keep = find(isnan(X)==0);
    X = X(keep);
    ie = ie(keep); ia = ia(keep);
    mm2 = length(keep);
    iShore2 = iShore(find(iWave>=keep(1),1,'first'):end);
    iWave2 = iWave(find(iWave>=keep(1),1,'first'):end)-keep(1)+1;
    t2 = t(keep)-t(keep(1));

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
%        B=[ones(mm,1),t(iWave),X(iWave),X2(iWave),ones(mm,1)];
        B=[ones(mm2,1),t2,X,X2];    %only calculates points where you have data

        x=B(iWave2,:)*A; 
        xr = data(iShore2)-x;   %adjust offset term if needed
        Atest = [ones(size(xr),1)]\xr;
        A(1)=A(1)+Atest; clear xr x
        x=B(iWave2,:)*A; 

       
        clear B
        xb = xb(keep);
    c=corrcoef(data(iShore2),x); cc=c(1,2);
    rmsError=sqrt(mean((data(iShore2)-x).^2));
    nmsError = sum((data(iShore2)-x).^2)./sum((data(iShore2)-mean(data(iShore2))).^2);
    BSS=brierss(data(iShore2),x,xb(iWave2),dx);
    AIC = akaikeIC(data(iShore2),x,5);
  
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
    plot(dates_survey(iShore2),x,'k','LineWidth',2)
    plot(dates_survey(iShore),xb(iWave2),'r','LineWidth',2)
    datetick('x','yyyy')
    xlabel('Time (years)')
    ylabel('\Delta shoreline position (m)')
    grid on
    legend('Data','Model','linear Model')

hindStats.dates=dates_survey;
hindStats.varShore = var(data(iShore));
hindStats.varOmega = var(omega);
hindStats.Nyr = (dates_survey(iShore(end))-dates_survey(iShore(1)))/365;
hindStats.Alin = Alin;
hindStats.A = A
hindStats.cc=cc;
hindStats.rmse=rmsError;
hindStats.nmse=nmsError;
hindStats.BSS=BSS;
hindStats.AIC=AIC;
hindStats.omegaRecordMean=omegaRecordMean;
hindStats.dx=dx;
hindStats.cr=cr(1,2);
hindStats.rmstrend=rmstrend;
hindStats.nmstrend=nmstrend;
hindStats.AIClin=AIClin;
hindStats.omegaFilter=omegaFiltered;
hindStats.F=F;
hindStats.k=k;
hindStats.dt=dt;
hindStats.omega=omega;
hindStats.x=x;
hindStats.contourShore=contShore;
hindStats.linModel = xb;


