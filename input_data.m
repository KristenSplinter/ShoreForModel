function input_data

global H T dnum data dates_survey w dx m mm dtShoreline iShore omega
global iCal iWave k ihardWire iWaveCal iConstC
global HCal TCal wCal dnumCal dataCal mCal dxCal mmCal datesCal
global h hCal A cutoff site omegaFilterCutOff contShore contShoreCal
global dataFile HFile filePathName iswitchCal filePathNameCal dataFileCal s1 s2 HFileCal
global tide tideCal tideSwitch tideSwitchCal beta betaCal zShore zShoreCal calStats hindStats hindDir

try
    load([pwd '\userData.mat'])
catch
    genUserData
    load('userData.mat')
end

%% HINDCAST DATA SET
    load(filePathName);
    data=eval(char(dataFile));  
    keep = find(isnan(data)==0);
    data=data(keep);
    dates_survey=eval(char(datesFile));
    [y,mo,d,hr,mi,s]=datevec(dates_survey(keep));
    dates_survey=datenum([y,mo,d,round(hr+mi./60),zeros(size(y)),zeros(size(y))]);
    H=eval(char(HFile));
    T=eval(char(TFile));
    dnum=eval(char(datesFileWave));
    clear y mo d hr mi s
    [y,mo,d,hr,mi,s]=datevec(dnum(:));
    dnum=datenum([y,mo,d,round(hr+mi./60),zeros(size(y)),zeros(size(y))]);
   
%      if strcmp(tideSwitch,'y')==1  
%         tide=eval(char(TideFile))+tideOffset;
%         datesTide=eval(char(datesFileTide));
%      end
     zShore=eval(char(zFile));
    ch=std(diff(dnum));
    if ch>1E-6
        dnum2 = [dnum(1):1/24:dnum(end)]';
        H2 = interp1(dnum,H,dnum2);
        T2 = interp1(dnum,T,dnum2);
%         if strcmp(tideSwitchCal,'y')==1  
%             tide2 = interp1(datesTide,tide,dnum2);
%             clear tide datesTide
%             tide = tide2; clear tide2
%          end
        clear H T dnum 
        H = H2; clear H2
        T = T2;, clear T2
        dnum = dnum2; clear dnum2
    end
  
%check data is equally spaced and only covers the time frame covered in
  
     w = fallvelocity(d50/1000,15);
     hindStats.d50=d50;
    disp(['fall velocity used (assumes tempT = 15 C) = ', num2str(w), 'm/s'])

% Length of time-series
    m=length(dnum);
    [~,iShore,iWave]=intersect(dates_survey,dnum);
    dtShoreline=24*mean(diff(dates_survey)); % Average shoreline sampling in hours
    mm=length(iShore);

%%CALIBRATION DATA SET
   if iswitchCal==0
       iCal=iShore;
       iWaveCal=iWave;
       HCal=H;
       TCal=T;
       dataCal=data;
       dnumCal=dnum;
       datesCal=dates_survey;
       wCal=w;
       calStats.d50=d50;
       dxCal=dx;
       mCal=length(dnumCal);
      mmCal=length(iCal);
      filePathNameCal=filePathName;
      HFileCal=HFile;
      dataFileCal=dataFile;
      zShoreCal=zShore;
      contShoreCal=contShore;
%       if strcmp(tideSwitch,'y')==1
%       TideFileCal=TideFile;
%       datesFileTideCal=datesFileTide;
%       tideOffsetCal=tideOffset;
%       betaCal=beta;
%       tideCal=tide;
%       end
%        tideSwitchCal=tideSwitch;

   elseif iswitchCal==1
       t1 = datenum(s1);
       t2 = datenum(s2);
       icalW = find(dnum>=t1 & dnum <= t2);
       dnumCal=dnum(icalW);
       HCal=H(icalW);
       TCal=T(icalW);
       icalS = find(dates_survey>=t1 & dates_survey <= t2);
       datesCal = dates_survey(icalS);
       dataCal=data(icalS);
       [~,iCal,iWaveCal]=intersect(dates_survey(icalS),dnum(icalW));
        wCal=w;
        calStats.d50=d50;
        dxCal=dx;
        mCal=length(dnumCal);
        mmCal=length(iCal);
       filePathNameCal=filePathName;
      HFileCal=HFile;
        dataFileCal=dataFile;
      zShoreCal=zShore(icalS);
      contShoreCal=contShore;
%       if strcmp(tideSwitch,'y')==1
%       tideCal=tide(icalW);
%       TideFileCal=TideFile;
%       datesFileTideCal=datesFileTide;
%       tideOffsetCal=tideOffset;
%       betaCal=beta;
%       end
%         tideSwitchCal=tideSwitch;
   else
       Htmp=H; Ttmp=T; datatmp=data; dnumtmp=dnum; tidetmp=tide; dates_surveytmp=dates_survey;
       % load file
       clear H T data dnum dates_survey keep
       load(filePathNameCal);
    dataCal=eval(char(dataFileCal));  
    keep = find(isnan(dataCal)==0);
    dataCal=dataCal(keep);
    datesCal=eval(char(datesFileCal));
    clear y mo d hr mi s
    [y,mo,d,hr,mi,s]=datevec(datesCal(keep));
    datesCal=datenum([y,mo,d,round(hr+mi/60),zeros(size(y)),zeros(size(y))]);
    HCal=eval(char(HFileCal));
    dnumCal=eval(char(datesFileWaveCal));
    zShoreCal=eval(char(zFileCal));
    clear y mo d hr mi s
   [y,mo,d,hr,mi,s]=datevec(dnumCal);
    dnumCal=datenum([y,mo,d,round(hr+mi./60),zeros(size(y)),zeros(size(y))]);

    TCal = eval(char(TFileCal));
%     if strcmp(tideSwitchCal,'y')==1
%     tideCal=eval(char(TideFileCal))+tideOffsetCal;
%     datesTideCal=eval(char(datesFileTideCal));
%     end
%check data is equally spaced and only covers the time frame covered in
  
    %shoreline data
    ch=std(diff(dnumCal));
    if ch>1E-10
        dnum2 = [dnumCal(1):1/24:dnumCal(end)]';
        H2 = interp1(dnumCal,HCal,dnum2);
        T2 = interp1(dnumCal,TCal,dnum2);
%      if strcmp(tideSwitchCal,'y')==1       
%         tide2 = interp1(datesTideCal,tideCal,dnum2);
%         clear tideCal
%         tideCal=tide2; clear tide2
%      end
        clear HCal  TCal dnumCal
        HCal = H2; clear H2
        TCal = T2;, clear T2
        dnumCal = dnum2; clear dnum2
    end

    wCal = fallvelocity(d50Cal/1000,15);
    disp(['fall velocity used (assumes tempT = 15 C) = ', num2str(wCal), 'm/s'])
    calStats.d50=d50Cal;
% Length of time-series
    mCal=length(dnumCal);
    [tBoth,iCal,iWaveCal]=intersect(datesCal,dnumCal);
    mmCal=length(iCal);
         %make sure full data sets are correct
        clear H T dates_survey dnum data
        H = Htmp;
        T = Ttmp;
        dates_survey=dates_surveytmp;
        dnum = dnumtmp;
        data=datatmp;
%        tide=tidetmp;
   end
