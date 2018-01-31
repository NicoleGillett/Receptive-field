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
idataset = 2; %example good one: 10
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

%% shift saccades

manual = dlmread('Manual2.csv'); %load manually counted saccades

StepThre = 17; %set step threshold
ns = length(StepThre);

saccades = findSaccades(eyepos_sync,nanSync,StepThre); %automatically find saccades
length_saccades = length(saccades);

hit = intersect(manual,saccades); %find hits between manually counted saccades and automatic detection
HIT = length(hit);
    
miss = setdiff(manual,hit); %saccades missed by the algorithm
MISS = length(miss);

rh = zeros(10,length_saccades); %'remove hits' matrix to store values left once hits have been removed
s = zeros(10,length_saccades); %the remaining values are shifted and stored in matrix s
s(1,:) = saccades;
h = zeros(10,HIT); %h stores the new number of hits following the shift
h(1,:) = hit;
m = zeros(10,MISS); %m stores the new number of misses following the shift
m(1,:) = miss;

for jj = 2:10 %set 10 shifts 
        
    last_str_index = jj-1;
    str_index = jj;
    
    last_s = s(last_str_index,:);
    last_h = h(last_str_index,:);
    current_rh = rh(str_index,:);
    last_m = m(last_str_index,:);
    current_s = s(str_index,:);
    current_h = h(str_index,:);
    
    remove_hits = setdiff(last_s,last_h); %remove hits detected by alrogithm on previous turn
    
    if mod(jj,2) == 0 %shifts -1 +2 -3 +4 etc. from starting position
        ii = -1;
    else
        ii = 1;
    end
    
    shift_saccades = remove_hits + (ii*(jj-1)); 
    remove_hits(1,length_saccades) = 0; %pad vector with zeros
    shift_saccades(1,length_saccades) = 0;
    rh(str_index,:) = remove_hits;
    s(str_index,:) = shift_saccades;
    new_hits = intersect(last_m,shift_saccades); %new hits by comparing last misses with shifted saccades 
    new_hits(1,HIT) = 0;
    h(str_index,:) = new_hits;
    new_misses = setdiff(last_m,new_hits); %new misses by comparing last misses with new hits
    new_misses(1,MISS) = 0;
    m(str_index,:) = new_misses; 
    
end
        
figure;
hold on

plot(timeFrames,eyepos_sync(1,:),'c') %plot time course of eye position
xlim([-inf inf])
ylim([-inf inf])
xlabel('Time (s)')
ylabel('x pos. (Deg)')

T=transpose(manual);

max_eye = max(eyepos,[],2);
min_eye = min(eyepos,[],2);

%plot manually detected saccades
for is = T
        plot((timeFrames(is-1)+timeFrames(is)).*[1 1]/2,[min_eye(1) max_eye(1)],'k')
end

%plot saccades detected by algorithm
for is = s(1,:)
        plot((timeFrames(is-1)+timeFrames(is)).*[1 1]/2,[min_eye(1) max_eye(1)],'r--','LineWidth',2)
end

%plot shifted saccades
for iz = 2:10
    for is = s(iz,(s(iz,:)~=0))
        plot((timeFrames(is-1)+timeFrames(is)).*[1 1]/2,[min_eye(1) max_eye(1)],'r:')
    end
end