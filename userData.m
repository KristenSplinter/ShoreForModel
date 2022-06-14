% File used to generate userData.mat that is read in from input_data

%Hindcast Data Set
site = ''               %site name
filePathNameShore = ''; %complete file name of hindcast data set (contains shore, time, waves)
dataFile = '';          %variable name (in filePathName) of shoreline data
datesFile = '';         %variable name for dates corresponding to shoreline data
HFile = '';             %variable name of wave height
TFile = '';             %variable name of wave period
datesFileWave = '';     %variable name of dates file for waves
d50 = [];               %median grain size, mm
dx = [];                %error in shoreline measurement, m. USed in Brier Skill score.
%tideSwtitch='';       %User switch (y/n) if you wish to include tide info (and R2) in shoreline exposure
%if tideSwitch=y
    beta=[];         %beach slope
%    tideFile='';     %tidal data file
%    datesFileTide='';   %tidal dates file
%    tideOffset=[];  %offset of datum relative to MSL
%hindDir='';         %'backward' or 'forward' direction in which you do the hindcast

%h = [];                 %depth of wave measurements

%Calibration options
iswitchCal=[];          %0 = full hindcast, 1 = partial hindcast, 2=different
%if iswitchCal=1
    s1 = '';                %Start date for partial hindcast (ex. Jan. 10, 2000)
    s2 = '';                %End date for partial hindcast (ex. Jan. 10, 2000)
%if iswitchCal=2
    filePathNameCal = '';   %complete file name of calibration data set (contains shore, time, waves)
    dataFileCal = '';          %variable name (in filePathName) of shoreline data
    datesFileCal = '';         %variable name for dates corresponding to shoreline data
    HFileCal = '';             %variable name of wave height
    TFileCal = '';             %variable name of wave period
    datesFileWaveCal = '';     %variable name of dates file for waves
    d50Cal= [];               %median grain size, mm
    dxCal=[];               %error in shoreline measurement, m. USed in Brier Skill score.
%    hCal=[];               %depth of wave measurements for cal data set
%if iswitchCal=3    %already know my calibration data
    cutoff = [];            %cutoff frequency (in hours)
    A = [];                   %Appropriate coef. matrix. based on your model 1 col. 
                            %should include a;b;cm;cp (DMT11) OR a;cm;cp
                            %(DMT11) OR a;b;c (DLT10 or DMT11)

%Model options
model = '';             %Model name (current options DLT10, DMT11)
ihardWire = [];         %0 for hard-wire filter and 1 for optimized
%if ihardWire=0
    omegaFilterCutoff=[];   %cut-off in days
iConstC=[];             %0 for variable erosion/accretion and 1 for constant
iLinTerm='';           %'n' for no linear trend (b), 'y' for include
k=0.5;                  % exponent in equation (hard wired)
tideSwitchCal='';       %User switch (y/n) if you wish to include tide info (and R2) in shoreline exposure
%if tideSwitchCal=y
    betaCal=[];         %beach slope
    tideFileCal='';     %tidal data file
    datesFileTideCal='';   %tidal dates file
    tideOffsetCal=[];  %offset of datum relative to MSL
save userData.mat