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



%% check threshold
StepThre=17;
saccades = findSaccades(eyepos_sync,nanSync,StepThre);

max_eye = max(eyepos,[],2);
min_eye = min(eyepos,[],2);

figure;
hold on
    
plot(timeFrames,eyepos_sync(1,:),'c')
    for is = saccades
        plot((timeFrames(is-1)+timeFrames(is)).*[1 1]/2,[min_eye(1) max_eye(1)],'k:')
    end
xlim([-inf inf])
ylim([-inf inf])
xlabel('Time (s)')
ylabel('x pos. (Deg)')
title('5')