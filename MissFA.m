%% clear workspace
clearvars;
close all;
clc;

%% enter in the data folder
%dataFolder = 'C:\Users\Mario\Dropbox\Data Nicole Gillett\';
dataFolder = '/Users/NicoleGillett/Documents/Neuroscience MSc/Code - receptive field/';
%% load session metadata
load([dataFolder 'Meta Data/metaData'])
fprintf('There are %d datasets\n',length(metaData.fileName));
%% choose session
idataset = 5; %example good one: 10
%%
addpath([dataFolder 'Eye Analysis'])
addpath([dataFolder 'PlottingToolbox'])
%% load raw data
load([dataFolder 'Raw Data/rawData_' metaData.fileName{idataset}])
%% load eye position and time frames

time_eye = rawData.time_eye;
nt = length(time_eye); %length of eye position vector

eyepos = rawData.eyepos; %horizontal (1) and vertical (2) eye position

timeFrames = rawData.timeFrames; %this can be arbitrary and outside the function
nf = length(timeFrames); %length of calcium frame vector 




%% interpolate nan values in eye position

nanVals = zeros(1,nt);
nanFrames = find(isnan(mean(eyepos))); %find nan values of eye position
realFrames = find(~isnan(mean(eyepos))); %find correct values of eye position
nanVals(nanFrames) = ones(1,length(nanFrames));

if sum(isnan(time_eye))==length(time_eye)
    error('Bug in eye time frames\n')
end
nanSync = interp1(time_eye, nanVals, timeFrames);

%% interpolate eyeposition using time of calcium frames
eyepos_sync = nan(2,nf);
for ix = 1:2
    eyepos_real = eyepos(ix,:);
    eyepos_real(nanFrames) = interp1(realFrames, eyepos(ix,realFrames), nanFrames, 'pchip');
    eyepos_sync(ix,:) = interp1(time_eye, eyepos_real, timeFrames);
end

%% compare misses and false alarms with manually counted saccades

M = dlmread('Manual.csv'); %load manually counted saccades

ST = 10.^(0:0.05:3); %set step threshold
ns = length(ST);
HIT = nan(1,ns);
MISS = nan(1,ns);
FA = nan(1,ns);

for ij = 1:ns

    StepThre = ST(ij);
    saccades = findSaccades(eyepos_sync,nanSync,StepThre); %automatically find saccades
 
    hit = intersect(M,saccades); %find hits between manual and automatic saccades
    HIT(ij) = length(hit);
    
    miss = setdiff(M,hit); %find misses from the difference between hits and manual saccades
    MISS(ij)= length(miss);
    
    s2 = saccades-1; %shift saccades
    h2 = intersect(miss,s2); %calculate new number of hits
    h2(ij) = length(h2);
    
    fa = setdiff(saccades,hit); %find false alarms from the difference between hits and automatic saccades
    FA(ij) = length(fa);
end
    %scale misses and false alarms
    MissM=MISS./length(M);
    FaM=FA./length(M);
    
    MissFa=MissM+FaM; %total number of misses and false alarms
    
    figure;
    hold on
    
    plot(ST,MissM,'o-')
    plot(ST,FaM,'o-')
    plot(ST,MissFa,'o-')
    xlabel('StepThre')
    ylabel('/manual')
    xlim([-inf inf])
    ylim([-inf 2])
    legend('MISS','FA','MISS+FA')
    set(gca,'XScale','Log')
    