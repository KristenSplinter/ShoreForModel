function output_results
global site ihardWire iLinTerm dnum dnumCal iShore data datesCal dates_survey iCal
global dataFile HFile filePathName iswitchCal filePathNameCal dataFileCal s1 HFileCal 
global  calStats hindStats m  hindDir vShoreFor



% Model/Data timestep - seconds - NB must be a regular time-step in data!
    dt=round((dnum(2)-dnum(1))*24*60*60);
    dtCal=round((dnumCal(2)-dnumCal(1))*24*60*60);
    dtShoreline=mean(diff(dates_survey));
    dtShoreCal = mean(diff(datesCal));
        

% Output model results to screen:
        clc
        disp('++++++++++++++++++++PROGRAM RESULTS+++++++++++++++++++++++++++++++++++++++++++++')
        disp(['Site modelled : ' site])
        disp(['Model used : ShoreFor ' vShoreFor])
        
        disp('++++++++++++++++++++CALIBRATION RESULTS+++++++++++++++++++++++++++++++++++++++++++++')
        disp(['File analysed : ' filePathNameCal])
        disp(['Shoreline data analysed : ' char(dataFileCal)])
        disp(['Shoreline contour elevation (m) : ' num2str(calStats.contourShore)])
        disp(['Median grain size, d50 (mm) : ' num2str(calStats.d50)])
        disp(['Wave height variable selected : ' char(HFileCal)])
        disp(['Shoreline data sampling interval (days) = ' num2str(dtShoreCal)])
        disp(['Wave data sampling interval (hours) = ' num2str(dtCal./(60*60))])        
 %       if (strcmp(model,'DMT11')==1)
       disp(['Omega low pass filter cutoff in days = D = ' num2str(calStats.D), ' phi = ', num2str(calStats.phi)])
 %       end
 %       disp(['Include tidal exposure index (y/n) : ' char(tideSwitchCal)])
        disp(['Omega record mean = ' num2str(calStats.omegaRecordMean)])
        disp(['Linear trend term in linear model (m/yr) = ' num2str(calStats.linModelA(2).*(60*60*24*365))])
%        if (iConstC ==0 )
            disp(['Optimal coefficients (a, b, c+, c-) = '])
            disp(['Shoreline offset (m) a = ' num2str(calStats.A(1))])
            disp(['Shoreline trend (m/yr) b = ' num2str(calStats.A(2).*(60*60*24*365))])
            disp(['Shoreline trend due to trend in waves (m/yr) bwave = ' num2str(calStats.waveTrend.*(60*60*24*365))])          
            disp(['Shoreline accretion response rate ((m/s).(W/m).^0.5) = c+ ' num2str(calStats.A(4))])
            disp(['Shoreline erosion response rate ((m/s).(W/m).^0.5) = c- ' num2str(calStats.A(3))])
            disp(['Average response rate ((m/s).(W/m).^0.5) = ' num2str(mean(calStats.A(3:4)))]) 
            disp(['Accretion:Erosion response rate ratio = ' num2str(calStats.A(4)./calStats.A(3))]) 
%        else
%             disp(['Optimal coefficients (a, b, c) = '])
%             disp(['Shoreline offset (m) = ' num2str(calStats.A(1))])
%             disp(['Shoreline trend (m/yr) = ' num2str(calStats.A(2).*(60*60*24*365))])
%             disp(['Shoreline response rate (m/s for DLT10 or (m/s).(W/m).^0.5 for DMT11) = ' num2str(calStats.A(3))])
%             disp(['Shoreline trend remaining (m/yr) = ' num2str(calStats.A(4).*(60*60*24*365))])
%         end
        disp(['RMS Error for ShoreFor(m) = ' num2str(calStats.rmse)])
        disp(['RMS Error for linear model (m) = ' num2str(calStats.rmstrend)])
        disp(['NMS Error for ShoreFor  = ' num2str(calStats.nmse)])
        disp(['NMS Error for linear model  = ' num2str(calStats.nmstrend)])
        disp(['Correlation coefficient = ' num2str(calStats.cc)])
        disp(['Linear trend corr coeff = ' num2str(calStats.cr)])
       disp(['Akaike Information Criteria ShoreFor = ' num2str(calStats.AIC)]);
        disp(['Akaike Information Criteria  linear model = ' num2str(calStats.AIClin)]);
        disp(['Brier Skills Score = ' num2str(calStats.BSS)])
        disp(['Error specified for shoreline position (m) = ' num2str(calStats.dx)]) 
        disp(['Number of survey data points used in calibration= ' num2str(length(iCal))])
        disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
         disp('++++++++++++++++++++HINDCAST RESULTS+++++++++++++++++++++++++++++++++++++++++++++')
        disp(['File analysed : ' filePathName])
        disp(['Shoreline data analysed : ' char(dataFile)])
        disp(['Shoreline contour elevation (m) : ' num2str(hindStats.contourShore)])
%        disp(['Include tidal exposure index (y/n) : ' char(tideSwitch)])
%        disp(['Hindcast direction (forward/backward) : ' char(hindDir)])
        disp(['Median grain size, d50 (mm) : ' num2str(hindStats.d50)])
        disp(['Wave height variable selected : ' char(HFile)])
        disp(['Shoreline data sampling interval (days) = ' num2str(dtShoreline)])
        disp(['Wave data sampling interval (hours) = ' num2str(dt./(60*60))])
        disp(['Error specified for shoreline position (m) = ' num2str(hindStats.dx)]) 
        disp(['Number of survey data points used in hindcast= ' num2str(length(iShore))])
        disp(['Omega record mean = ' num2str(hindStats.omegaRecordMean)])
        disp(['Optimal coefficients (a, b, c+, c-) = '])
        disp(['Shoreline offset (m) a = ' num2str(hindStats.A(1))])
        disp(['Shoreline trend (m/yr) b = ' num2str(hindStats.A(2).*(60*60*24*365))])
        disp(['Shoreline accretion response rate ((m/s).(W/m).^0.5) = c+ ' num2str(hindStats.A(4))])
        disp(['Shoreline erosion response rate ((m/s).(W/m).^0.5) = c- ' num2str(hindStats.A(3))])
        disp(['Average response rate ((m/s).(W/m).^0.5) = ' num2str(mean(hindStats.A(3:4)))]) 
        disp(['Accretion:Erosion response rate ratio = ' num2str(hindStats.A(4)./hindStats.A(3))]) 

       disp('++++++LINEAR MODEL+++++++++++++++++++++++++++++++++++++++++++++')
%        disp(['Linear trend term in linear model (m/yr) = ' num2str(hindStats.bTrend(2).*(60*60*24*365))])
        disp(['RMS Error for linear model (m) = ' num2str(hindStats.rmstrend)])
        disp(['NMS Error for linear model (m) = ' num2str(hindStats.nmstrend)])

        disp(['Linear trend corr coeff = ' num2str(hindStats.cr)])
        disp(['Akaike Information Criteria  linear model = ' num2str(hindStats.AIClin)]);

        disp('+++++++ShoreFor +++++++++++++++++++++++++++++++')
        disp(['RMS Error for new model (m) = ' num2str(hindStats.rmse)])
        disp(['NMS Error for new model (m) = ' num2str(hindStats.nmse)])
        disp(['Correlation coefficient = ' num2str(hindStats.cc)])
        disp(['Brier Skills Score = ' num2str(hindStats.BSS)])
       disp(['Akaike Information Criteria ShoreFor = ' num2str(hindStats.AIC)]);

%        disp('+++++++FULL MODEL (with trend term) +++++++++++++++++++++++++++++++')
%         disp(['RMS Error for new model (m) = ' num2str(hindStats.rmseFull)])
%         disp(['Correlation coefficient = ' num2str(hindStats.ccFull)])
%         disp(['Brier Skills Score = ' num2str(hindStats.BSSFull)])
%        disp('+++++++FULL MODEL (without trend term) +++++++++++++++++++++++++++++++')
%         disp(['RMS Error for new model (m) = ' num2str(hindStats.rmseShort)])
%         disp(['Correlation coefficient = ' num2str(hindStats.ccShort)])
%         disp(['Brier Skills Score = ' num2str(hindStats.BSSShort)])
%         disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
% 
        
 % Save output first go is quick and dirty - i.e. save whole workspace!
 %uisave 
 fid = fopen(['Results_' site '.txt'],'w')
  disp('Saving output to file')
% Output model results to file:
        fprintf(fid, '%s\r\n', '++++++++++++++++++++PROGRAM RESULTS+++++++++++++++++++++++++++++++++++++++++++++');
        fprintf(fid, '%s\r\n', ['Site modelled : ' site]);
        fprintf(fid, '%s\r\n', ['Model used : ShoreFor ' vShoreFor]);
        
        fprintf(fid, '%s\r\n', '++++++++++++++++++++CALIBRATION RESULTS+++++++++++++++++++++++++++++++++++++++++++++');
        fprintf(fid, '%s\r\n', ['File analysed : ' filePathNameCal]);
        fprintf(fid, '%s\r\n', ['Shoreline data analysed : ' char(dataFileCal)]);
        fprintf(fid, '%s\r\n', ['Shoreline contour elevation (m) : ' num2str(calStats.contourShore)]);
        fprintf(fid, '%s\r\n', ['Median grain size, d50 (mm) : ' num2str(calStats.d50)]);
        fprintf(fid, '%s\r\n', ['Wave height variable selected : ' char(HFileCal)]);
        fprintf(fid, '%s\r\n', ['Shoreline data sampling interval (days) = ' num2str(dtShoreCal)]);
        fprintf(fid, '%s\r\n', ['Wave data sampling interval (hours) = ' num2str(dtCal./(60*60))]);        
 %       if (strcmp(model,'DMT11')==1)
            fprintf(fid, '%s\r\n', ['Omega low pass filter cutoff in days = D = ' num2str(calStats.D), ' phi = ', num2str(calStats.phi)])
 %       end
 %       fprintf(fid, '%s\r\n', ['Include tidal exposure index (y/n) : ' char(tideSwitchCal)])
        fprintf(fid, '%s\r\n', ['Omega record mean = ' num2str(calStats.omegaRecordMean)]);
        fprintf(fid, '%s\r\n', ['Linear trend term in linear model (m/yr) = ' num2str(calStats.linModelA(2).*(60*60*24*365))]);
%        if (iConstC ==0 )
            fprintf(fid, '%s\r\n', ['Optimal coefficients (a, b, c+, c-) = ']);
            fprintf(fid, '%s\r\n', ['Shoreline offset (m) = ' num2str(calStats.A(1))]);
            fprintf(fid, '%s\r\n', ['Shoreline trend (m/yr) = ' num2str(calStats.A(2).*(60*60*24*365))]);
            fprintf(fid, '%s\r\n', ['Shoreline accretion response rate ((m/s)/(W/m).^0.5) = ' num2str(calStats.A(4))]);
            fprintf(fid, '%s\r\n', ['Shoreline erosion response rate ((m/s)/(W/m).^0.5) = ' num2str(calStats.A(3))]);
            fprintf(fid, '%s\r\n', ['Average response rate ((m/s)/(W/m).^0.5) = ' num2str(mean(calStats.A(3:4)))]) ;
            fprintf(fid, '%s\r\n', ['Accretion:Erosion response rate ratio = ' num2str(calStats.A(4)./calStats.A(3))]) ;
%        else
%             fprintf(fid, '%s\r\n', ['Optimal coefficients (a, b, c) = '])
%             fprintf(fid, '%s\r\n', ['Shoreline offset (m) = ' num2str(calStats.A(1))])
%             fprintf(fid, '%s\r\n', ['Shoreline trend (m/yr) = ' num2str(calStats.A(2).*(60*60*24*365))])
%             fprintf(fid, '%s\r\n', ['Shoreline response rate (m/s for DLT10 or (m/s).(W/m).^0.5 for DMT11) = ' num2str(calStats.A(3))])
%             fprintf(fid, '%s\r\n', ['Shoreline trend remaining (m/yr) = ' num2str(calStats.A(4).*(60*60*24*365))])
%         end
        fprintf(fid, '%s\r\n', ['RMS Error for ShoreFor (m) = ' num2str(calStats.rmse)]);
        fprintf(fid, '%s\r\n', ['RMS Error for linear model (m) = ' num2str(calStats.rmstrend)]);
        fprintf(fid, '%s\r\n', ['NMS Error for ShoreFor  = ' num2str(calStats.nmse)]);
        fprintf(fid, '%s\r\n', ['NMS Error for linear model  = ' num2str(calStats.nmstrend)]);
        fprintf(fid, '%s\r\n', ['Correlation coefficient = ' num2str(calStats.cc)]);
        fprintf(fid, '%s\r\n', ['Linear trend corr coeff = ' num2str(calStats.cr)]);
        fprintf(fid, '%s\r\n', ['Akaike Information Criteria ShoreFor = ' num2str(calStats.AIC)]);
        fprintf(fid, '%s\r\n', ['Akaike Information Criteria  linear model = ' num2str(calStats.AIClin)]);
        fprintf(fid, '%s\r\n', ['Brier Skills Score = ' num2str(calStats.BSS)]);
        fprintf(fid, '%s\r\n', ['Error specified for shoreline position (m) = ' num2str(calStats.dx)]) ;
        fprintf(fid, '%s\r\n', ['Number of survey data points used in calibration= ' num2str(length(iCal))]);
        fprintf(fid, '%s\r\n', '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
         fprintf(fid, '%s\r\n', '++++++++++++++++++++HINDCAST RESULTS+++++++++++++++++++++++++++++++++++++++++++++');
        fprintf(fid, '%s\r\n', ['File analysed : ' filePathName]);
        fprintf(fid, '%s\r\n', ['Shoreline data analysed : ' char(dataFile)]);
        fprintf(fid, '%s\r\n', ['Shoreline contour elevation (m) : ' num2str(hindStats.contourShore)]);
%        fprintf(fid, '%s\r\n', ['Include tidal exposure index (y/n) : ' char(tideSwitch)])
%        fprintf(fid, '%s\r\n', ['Hindcast direction (forward/backward) : ' char(hindDir)])
        fprintf(fid, '%s\r\n', ['Median grain size, d50 (mm) : ' num2str(hindStats.d50)]);
        fprintf(fid, '%s\r\n', ['Wave height variable selected : ' char(HFile)]);
        fprintf(fid, '%s\r\n', ['Shoreline data sampling interval (days) = ' num2str(dtShoreline)]);
        fprintf(fid, '%s\r\n', ['Wave data sampling interval (hours) = ' num2str(dt./(60*60))]);
        fprintf(fid, '%s\r\n', ['Error specified for shoreline position (m) = ' num2str(hindStats.dx)]) ;
        fprintf(fid, '%s\r\n', ['Number of survey data points used in hindcast= ' num2str(length(iShore))]);
        fprintf(fid, '%s\r\n', ['Omega record mean = ' num2str(hindStats.omegaRecordMean)]);
        fprintf(fid, '%s\r\n', ['Shoreline offset (m) a = ' num2str(hindStats.A(1))])
        fprintf(fid, '%s\r\n', ['Shoreline trend (m/yr) b = ' num2str(hindStats.A(2).*(60*60*24*365))])
        fprintf(fid, '%s\r\n', ['Shoreline accretion response rate ((m/s).(W/m).^0.5) = c+ ' num2str(hindStats.A(4))])
        fprintf(fid, '%s\r\n', ['Shoreline erosion response rate ((m/s).(W/m).^0.5) = c- ' num2str(hindStats.A(3))])
        fprintf(fid, '%s\r\n', ['Average response rate ((m/s).(W/m).^0.5) = ' num2str(mean(hindStats.A(3:4)))]) 
        fprintf(fid, '%s\r\n', ['Accretion:Erosion response rate ratio = ' num2str(hindStats.A(4)./hindStats.A(3))]) 

        fprintf(fid, '%s\r\n', '++++++LINEAR MODEL+++++++++++++++++++++++++++++++++++++++++++++');
%        fprintf(fid, '%s\r\n', ['Linear trend term in linear model (m/yr) = ' num2str(hindStats.bTrend(2).*(60*60*24*365))]);
        fprintf(fid, '%s\r\n', ['RMS Error for linear model (m) = ' num2str(hindStats.rmstrend)]);
        fprintf(fid, '%s\r\n', ['NMS Error for linear model (m) = ' num2str(hindStats.nmstrend)]);
       fprintf(fid, '%s\r\n', ['Akaike Information Criteria  linear model = ' num2str(hindStats.AIClin)]);
        fprintf(fid, '%s\r\n', ['Linear trend corr coeff = ' num2str(hindStats.cr)]);
        fprintf(fid, '%s\r\n', '+++++++ShoreFor +++++++++++++++++++++++++++++++');
        fprintf(fid, '%s\r\n', ['RMS Error for ShoreFor(m) = ' num2str(hindStats.rmse)]);
        fprintf(fid, '%s\r\n', ['NMS Error for ShoreFor  = ' num2str(hindStats.nmse)]);
        fprintf(fid, '%s\r\n', ['Correlation coefficient = ' num2str(hindStats.cc)]);
        fprintf(fid, '%s\r\n', ['Brier Skills Score = ' num2str(hindStats.BSS)]);
        fprintf(fid, '%s\r\n', ['Akaike Information Criteria   = ' num2str(hindStats.AIC)]);

        fclose(fid)
           tmp = [dtShoreCal length(iCal) calStats.Nyr calStats.varShore calStats.varOmega calStats.omegaRecordMean calStats.phi hindStats.A(1) calStats.A(2).*(60*60*24*365) calStats.waveTrend.*(60*60*24*365)...
       calStats.A(4) calStats.A(3) calStats.A(4)./calStats.A(3) hindStats.cc hindStats.BSS hindStats.AIC-hindStats.AIClin hindStats.rmse hindStats.nmse hindStats.Nyr hindStats.varShore hindStats.varOmega];
   fid = fopen('OutputStats.txt','a');
    fprintf(fid, '%6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.0f\t %6.2f\t %6.2f\t %6.2f\t %6.4e\t %6.4e\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.4f\t %6.2f\t %6.2f\t %6.2f\n',tmp)
   fclose(fid)
