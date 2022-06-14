function genUserData
clc
site=input('What site are you testing? ','s')
disp('Which data set do you wish to use to TEST (hindcast) the model?')
% load file
    [FileName,PathName,FilterIndex] = uigetfile('*.mat','Please choose a hindcast data file (*.mat)');
    filePathName=[PathName,FileName];
%hindDir = input('What direction do you want to hindcast in? (forward or backward)','s')
    % sort out datafiles and wave height data for selection
    fileContents=who('-file', filePathName);        
    jj=1;
    kk=1;
    zz=1;
    yy=1;
    mm=1;
    kk2=1;
    for ii=1:length(fileContents)
        character=char(fileContents(ii));
        % look for data files
        if length(character)>=4 & character(1:4)=='data';
            dataFileTmp(jj)=fileContents(ii);
            jj=jj+1;
        end
        % Look for wave data
        if  character(1:1)=='H';
            HFileTmp(kk)=fileContents(ii);
            kk=kk+1;
        end
         if  character(1:1)=='T';
            TFileTmp(kk2)=fileContents(ii);
            kk2=kk2+1;
        end
        % Look for shore elevation data
        if  character(1:1)=='z';
            zFileTmp(mm)=fileContents(ii);
            mm=mm+1;
        end
        % Look for date data
        if  length(character)>=5 & character(1:5)=='dates';
            datesFileTmp(zz)=fileContents(ii);
            zz=zz+1;
        end
        % Look for tide data
%         if  length(character)>=4 & lower(character(1:4))=='tide';
%             TideFileTmp(yy)=fileContents(ii);
%             yy=yy+1;
%         end
    end
     % display shoreline data files and select data to be analysed  
    nf=length(dataFileTmp);
    if nf==1
        fswitch=1;
        dataFile=dataFileTmp(fswitch);
    else
        disp('Enter data file number ')
        for ii=1:nf
            disp(['file No.' num2str(ii)   dataFileTmp(ii)])    
        end
        fswitch=input('');
        dataFile=dataFileTmp(fswitch);
    end
    nz=length(zFileTmp);
    if nz==1
        fswitch=1;
        zFile=zFileTmp(fswitch);
    else
        disp('Enter data file number ')
        for ii=1:nz
            disp(['file No.' num2str(ii)   zFileTmp(ii)])    
        end
        fswitch=input('');
        zFile=zFileTmp(fswitch);
    end
    contShore = input('What elevation is this contour wrt tidal datum (m)? ')
    nT=length(datesFileTmp);
    if nT==1
        fswitch=1;
       datesFile=datesFileTmp(fswitch);      
    else
        disp('Enter date file number for shoreline ')
        for ii=1:nT
            disp(['file No.' num2str(ii)   datesFileTmp(ii)])    
        end
        fswitch=input('');
        datesFile=datesFileTmp(fswitch);  
    end
    
 % display H available data and select data to be analysed 
    nH=length(HFileTmp);
    if nH==1
        Hswitch=1;
        HFile=HFileTmp(Hswitch);
    else
        disp('Enter data file number ')
        for ii=1:nH
            disp(['file No.' num2str(ii)   HFileTmp(ii)])    
        end
        Hswitch=input('');
        HFile=HFileTmp(Hswitch);
    end
        nTT=length(TFileTmp);
    if nTT==1
        Tswitch=1;
        TFile=TFileTmp(Tswitch);
    else
        disp('Enter data file number ')
        for ii=1:nTT
            disp(['file No.' num2str(ii)   TFileTmp(ii)])    
        end
        Tswitch=input('');
        TFile=TFileTmp(Tswitch);
    end
nT=length(datesFileTmp);
    if nT==1
        fswitch=1;
        datesFileWave=datesFileTmp(fswitch);              
    else
        disp('Enter date file number for waves ')
        for ii=1:nT
            disp(['file No.' num2str(ii)   datesFileTmp(ii)])    
        end
        fswitch=input('');
        datesFileWave=datesFileTmp(fswitch);
    end
%    tideSwitch=input('Do you wish to include tidal exposure index (This requires tide info) (y/n)? ','s');
%     if strcmp(tideSwitch,'y')==1
%          if nT==1
%             fswitch=1;
%             datesFileTide=datesFileTmp(fswitch);              
%          else
%             disp('Enter date file number for Tides ')
%             for ii=1:nT
%                 disp(['file No.' num2str(ii)   datesFileTmp(ii)])    
%             end
%             fswitch=input('');
%             datesFileTide=datesFileTmp(fswitch);
%          end 
%           nTi=length(TideFileTmp);
%     if nTi==1
%         Tiswitch=1;
%         TideFile=TideFileTmp(Tiswitch);
%     else
%         disp('Enter Tide file number ')
%         for ii=1:nTi
%             disp(['file No.' num2str(ii)   TideFileTmp(ii)])    
%         end
%         Tiswitch=input('');
%         TideFile=TideFileTmp(Tiswitch);
%     end
%    tideOffset=input('Enter vertical offset (m) in datum from MSL (ex. NGVD at Duck = -0.1926m): ');
%     beta = input('Enter mean beach slope (to use in runup calc) : ')
    %end
    %h = input('What water depth is the wave data measured in?')
   d50 = input('Please input d50 (mm) ');
   dx=input('Please input the error in shoreline measurements (m) ');

 %% CALIBRATION DATA SET
   disp('Which data set do you wish to use to CALIBRATE the model?')
   disp('--- Input 0 for FULL data set')
   disp('--- Input 1 for PARTIAL data set')
   disp('--- Input 2 for DIFFERENT data set')
   disp('--- Input 3 for use known calibration coefficients')
   iswitchCal= input('');
   if iswitchCal==1
       s1=input('Enter start date for calibration (ex. Jan 01, 2000): ', 's');
       s2=input('Enter end date for calibration (ex. Jan 01, 2000): ', 's');
   elseif iswitchCal==2
        [FileNameCal,PathNameCal,FilterIndexC] = uigetfile('*.mat','Please choose a Calibration shoreline data file (*.mat)');
        filePathNameCal=[PathNameCal,FileNameCal];
       % sort out datafiles and wave height data for selection
        fileContents=who('-file', filePathNameCal);        
        jj=1;
        kk=1;
        zz=1;
        yy=1;
        mm=1;
        nn=1;
        clear dataFileTmp HFileTmp datesFileTmp TideFileTmp zFileTmp TFileTmp
        for ii=1:length(fileContents)
            character=char(fileContents(ii));
            % look for data files
            if length(character)>=4 & character(1:4)=='data';
                dataFileTmp(jj)=fileContents(ii);
                jj=jj+1;
            end
            % Look for wave data
            if  character(1:1)=='H';
                HFileTmp(kk)=fileContents(ii);
                kk=kk+1;
            end
            % Look for wave data
            if  character(1:1)=='T';
                TFileTmp(nn)=fileContents(ii);
                nn=nn+1;
            end
            % Look for date data
            if  length(character)>=5 & character(1:5)=='dates';
                datesFileTmp(zz)=fileContents(ii);
                zz=zz+1;
            end
                    % Look for shore elevation data
        if  character(1:1)=='z';
            zFileTmp(mm)=fileContents(ii);
            mm=mm+1;
        end
%                    % Look for tide data
%         if  length(character)>=4 & character(1:4)=='Tide';
%             TideFileTmp(yy)=fileContents(ii);
%             yy=yy+1;
%         end

        end
         % display shoreline data files and select data to be analysed  
        nf=length(dataFileTmp);
        if nf==1
            fswitch=1;
            dataFileCal=dataFileTmp(fswitch);
        else
            disp('Enter data file number ')
            for ii=1:nf
                disp(['file No.' num2str(ii)   dataFileTmp(ii)])    
            end
            fswitch=input('');
            dataFileCal=dataFileTmp(fswitch);
        end
           nz=length(zFileTmp);
        if nz==1
        fswitch=1;
        zFileCal=zFileTmp(fswitch);
    else
        disp('Enter data file number ')
        for ii=1:nz
            disp(['file No.' num2str(ii)   zFileTmp(ii)])    
        end
        fswitch=input('');
        zFileCal=zFileTmp(fswitch);
    end
            contShoreCal = input('What elevation is this contour wrt tidal datum (m)? ')
        nT=length(datesFileTmp);
        if nT==1
            fswitch=1;
           datesFileCal=datesFileTmp(fswitch);      
        else
            disp('Enter date file number for shoreline ')
            for ii=1:nT
                disp(['file No.' num2str(ii)   datesFileTmp(ii)])    
            end
            fswitch=input('');
            datesFileCal=datesFileTmp(fswitch);  
        end

     % display H available data and select data to be analysed 
        nH=length(HFileTmp);
        if nH==1
            Hswitch=1;
            HFileCal=HFileTmp(Hswitch);
        else
            disp('Enter data file number ')
            for ii=1:nH
                disp(['file No.' num2str(ii)   HFileTmp(ii)])    
            end
            Hswitch=input('');
            HFileCal=HFileTmp(Hswitch);
        end
             % display H available data and select data to be analysed 
        nTw=length(TFileTmp);
        if nTw==1
            Tswitch=1;
            TFileCal=TFileTmp(Tswitch);
        else
            disp('Enter data file number ')
            for ii=1:nTw
                disp(['file No.' num2str(ii)   TFileTmp(ii)])    
            end
            Tswitch=input('');
            TFileCal=TFileTmp(Tswitch);
        end
    nT=length(datesFileTmp);
        if nT==1
            fswitch=1;
            datesFileWaveCal=datesFileTmp(fswitch);              
        else
            disp('Enter date file number for waves ')
            for ii=1:nT
                disp(['file No.' num2str(ii)   datesFileTmp(ii)])    
            end
            fswitch=input('');
            datesFileWaveCal=datesFileTmp(fswitch);
        end
         %hCal = input('What water depth is the wave data measured in?')
%              tideSwitchCal=input('Do you wish to include tidal exposure index (This requires tide info) (y/n)? ','s');
%     if strcmp(tideSwitchCal,'y')==1
%          if nT==1
%             fswitch=1;
%             datesFileTideCal=datesFileTmp(fswitch);              
%          else
%             disp('Enter date file number for Tides ')
%             for ii=1:nT
%                 disp(['file No.' num2str(ii)   datesFileTmp(ii)])    
%             end
%             fswitch=input('');
%             datesFileTideCal=datesFileTmp(fswitch);
%          end 
%           nTi=length(TideFileTmp);
%     if nTi==1
%         Tiswitch=1;
%         TideFileCal=TideFileTmp(Tiswitch);
%     else
%         disp('Enter Tide file number ')
%         for ii=1:nTi
%             disp(['file No.' num2str(ii)   TideFileTmp(ii)])    
%         end
%         Tiswitch=input('');
%         TideFileCal=TideFileTmp(Tiswitch);
%     end
%     tideOffsetCal=input('Enter vertical offset (m) in datum from MSL (ex. NGVD at Duck = -0.1926m): ');
%         betaCal = input('Enter mean beach slope (to use in runup calc: ')
% 
%     end

       d50Cal = input('Please input d50 (mm) ');
       dxCal=input('Please input the error in shoreline measurements (m) ')
       
   end
    
   %Model info
   %    model=input('Which shoreline model do you want (DLT10 or DMT11)? ','s')
   %    if strcmp(model,'DMT11')==1
       if iswitchCal~=3
           disp('Do you wish to hard-wire filter cut-off?')
            disp('--- Input 0 for hard-wire filter')
            disp('--- Input 1 for optimise filter')
            ihardWire=input('');
            if (ihardWire==0)
                omegaFilterCutOff=input('Input filter value in days ')
            end
       end
              %iLinTerm=input('Do you wish to include a best fit linear
              %trend term (b), y/n?','s')  %obsolete in this version
   %    end

   %    disp('Response Rate Parameter (c or c+/-)')
   %    disp('--- Input 0 for variable erosion & accretion response rate coefficients (c+/c-) ')
   %    disp('--- Input 1 for constant response rate coefficient (c) ')
   %    iConstC= input('');
        iConstC=0;
   if iswitchCal==3
           disp(['Enter Optimal coefficients (c+, c-) = '])
            %a=input('Enter Shoreline offset  (m): a = ')
%             if strcmp(iLinTerm,'y')==1
%                 b=input('Shoreline trend (m/yr): b = ')
%             end
%            if iConstC==0
                cp=input('Enter Shoreline accretion response rate (m/s for DLT10 or (m/s).(W/m).^0.5 for DMT11): c+ = ')
                cm=input('Enter Shoreline erosion response rate (m/s for DLT10 or (m/s).(W/m).^0.5 for DMT11): c- = ')
%            else
%                c = input('Enter Shoreline response rate (m/s for DLT10 or (m/s).(W/m).^0.5 for DMT11): c = ')
%            end
            %if (strcmp(iLinTerm,'y')==1 & iConstC==0)
            %    A = [a;b;cm;cp];
           % elseif(strcmp(iLinTerm,'y')==1 & iConstC==1)
            %    A = [a;b;c];
            %elseif(strcmp(iLinTerm,'n')==1 & iConstC==0)
            %    A = [a;cm;cp];
            %end
 %           if(iConstC==0)  %%I think this still needs work.
                A(3:4) = [cm; cp];
                A(1:2)=[];
 %           else
 %               A(3)=[c];
 %               A(1:2)=[];
 %           end
            cutoff=input('Enter optimal cutoff frequency (hours) ')
       end
 k=0.5  % exponent in equation (hard wired)

save userData.mat