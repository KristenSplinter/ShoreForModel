% SHOREFOR v 1.2
% Main program that calls subroutines to calibrate, hindcast, forecast and
% output shoreline predictions for a given site. Main model based on DMT11.
% Mark Davidson, Melissa Mole, Kristen Splinter, Ian Turner
% Nov. 10, 2011

% Additional files needed to run this programme include:
% Data files: RichardsData.mat (Gold Coast) and MMTest2.mat (Narrabeen)
% fftfilter.m - FFT filter subroutine by Mark Davidson
% brierss.m - Brierr Skills score subroutine by Mark Davidson
% Modified 25/5/11 to only use data recorded at the morphological time-step
% for calibration and skills scores. (little difference made for weekly
% samples)

% clear everything
%     clear variables; close all;clc;
%     clearvars -global
    format SHORT ENG
% version number
    vShoreFor='v1.0';
% an input file (userData.mat) listing file names and variables can be 
% created instead of manullay inputting this info in input_data. A matlab 
% program (genUserData.m)is providded for this purpose. input_data
% checks if this file exists and if it doesn't, prompts user for info.

% list all variables needed to be passed among codes here
global H T dnum data dates_survey w 
global  A cutoff vShoreFor
global HCal TCal wCal dnumCal dataCal mCal dxCal 
global  mCal dates_survey iswitchCal
global  hindStats calStats iShore site

disp('importing data')
input_data
%check sizes
if length(data) ~=length(dates_survey); disp('Warning: Length of date data file does not match the length of data file, programme aborted'); break; end    
if length(H) ~=length(dnum); disp('Warning: Length of date data file does not match the length of wave file, programme aborted'); break; end    

if iswitchCal~=3
disp('model calibration')
[A, cutoff]=calibrate_modelWSG85(HCal,TCal,wCal,dnumCal,dataCal,mCal,dxCal)
end
disp('hindcasting data')
%hindcast_dataExplicit(cutoff,A,H,T,w,dnum,data-mean(data))
hindcast_data(cutoff,A,H,T,w,dnum,data)

disp('forecasting data')
%forecast
theta=zeros(size(H));
%forecast_dataSim(w, H, T,theta,dnum,data-mean(data),cutoff,A)
output_results
