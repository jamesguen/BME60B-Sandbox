clear all; clc; close all;

%% setup

% Ask user to choose the .mat file
disp('Select valid ECG .mat file.');
[fileName, filePath] = uigetfile('*.mat', 'Select the ECG .mat file');
if isequal(fileName, 0)
    disp('No file selected. Exiting...');
    return;
else
    clc
    disp(' Thank you! Now analyzing...');
    disp(' ');
end
load(fullfile(filePath, fileName));

% Allowed graph numbers
allowedGraphs = [1:5, 7:23, 25:30, 150:162];

% Prompt the user to enter a graph number
while true
    disp('Valid graph numbers span [1-5,7-23,25-30, 150-162]. You may choose a graph number');
    disp('within these boundaries. (Graphs not included are too unstable for analysis');
    disp(' ');
    Graph = input('Enter a graph number: ');
    
    % Check if the input is valid
    if ismember(Graph, allowedGraphs)
        fprintf('You entered a valid graph number: %d\n', Graph);
        disp(' ');
        disp('Now analyzing...');
        break; % Exit the loop for valid input
    else
        fprintf('Invalid input. Please try again.\n');
    end
end

Data=ECGData.Data(Graph,:);
fs=128; %max frequency is 64Hz
%Data=highpass(Data,1,fs); %filters frequencies to isolate the R spike

%% Filtering and R wave Detection

HPF=designfilt('highpassfir','StopbandFrequency',1, 'PassbandFrequency',2,SampleRate=fs); % Bandpass Filter Parameters
dataFilt=filtfilt(HPF,Data); % Data Filtering
time=0:1/fs:size(dataFilt,2)/fs-1/fs;
[yR,RpeakX]=findpeaks((dataFilt.^3),"MinPeakHeight",2*std(dataFilt.^3),'MinPeakDistance',10); % R wave detection
[nyR,nRpeakX]=findpeaks(-1*(dataFilt.^3),"MinPeakHeight",2*std(dataFilt.^3),'MinPeakDistance',10); % R wave detection
if length(nRpeakX)>length(RpeakX)
    RpeakX=nRpeakX;yR=nyR*-1;
    dataFilt=dataFilt*-1;
    negtativeEcg=1;
else
    negtativeEcg=0;
end
clc

%% Ploting Full ECG Data and Detected R-waves

% Prompt the user for input on whether to display the figure
userInput = input('Would you like to see the ECG plot with detected R waves? (y/n): ', 's');

% Brief description of the figure
disp(' ');
disp('Figure 1, titled "ECG Signal with Detected R Waves", depicts the filtered ECG signal');
disp('with detected R wave markers. The filtered ECG signal is shown in blue, and the R wave');
disp('markers are displayed in red.');
disp(' ');
if lower(userInput) == 'y'
    % Plotting ECG Data with R Wave Markers (Cleaned)
    figure; % Initialize a new figure
    hold on;

    % Plot the filtered ECG data with a thinner line
    plot(time, dataFilt, 'b', 'LineWidth', 0.8, 'DisplayName', 'Filtered ECG Signal');

    % Plot the R wave markers with a distinct color and smaller size
    scatter(time(RpeakX), dataFilt(RpeakX), 30, 'r', 'filled', 'DisplayName', 'R Wave');

    % Add a title and labels
    title('ECG Signal with Detected R Waves', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Time (s)', 'FontSize', 10);
    ylabel('Amplitude (mV)', 'FontSize', 10);

    % Adjust x-axis limits for better visualization
    xlim([0, max(time)]);

    % Improve readability with a grid
    grid on;

    % Add a legend in a non-intrusive location
    legend('Location', 'best', 'FontSize', 9);

    % Set a clean background color
    set(gca, 'Color', [0.95, 0.95, 0.95]); % Light gray background

    % Finalize the plot
    hold off;
else
    disp('You chose not to display the figure.');
end

%% Wave Extraction (Ploting) for Visualization

% Define the window around the R wave to extract
preRspike = round(fs * 0.440); % Duration (in samples) before the R wave
postRspike = round(fs * 0.500); % Duration (in samples) after the R wave

% Total number of detected R wave peaks
numSpikes = length(RpeakX);

% Initialize matrix to store individual cardiac cycles
% Each row corresponds to one cycle, spanning preRspike + postRspike samples
wave = zeros(numSpikes - 2, preRspike + postRspike + 1);

% Ask the user if they would like to see the figure
userInput = input('Would you like to see the figure showing extracted cardiac cycles and the mean waveform? (y/n): ', 's');

% Provide a brief description of the figure
disp(' ');
disp('Figure 2, titled "Extracted Cardiac Cycles and Mean Waveform", displays ');
disp('individual cardiac cycles (in black) extracted from the ECG signal, ');
disp('along with the average waveform (in red) that represents the mean cycle.');
disp(' ');

% Create a figure for visualizing individual cardiac cycles if the user chooses 'y'
if lower(userInput) == 'y'
    figure; 
    hold on;

    % Loop through each detected R wave (skipping the first and last for safety)
    for i = 2:numSpikes-1
        % Define the start and end indices for the current cardiac cycle
        idxStart = RpeakX(i) - preRspike;
        idxEnd = RpeakX(i) + postRspike;

        % Extract the current cycle and store it in the matrix
        wave(i, :) = dataFilt(idxStart:idxEnd);

        % Plot the extracted cycle in black for visualization
        plot(dataFilt(idxStart:idxEnd), 'k', 'DisplayName', sprintf('Cycle %d', i - 1));
    end

    % Compute the average wave across all cycles for visualization
    meanWave = mean(wave, 1);

    % Overlay the average wave in red for emphasis
    plot(meanWave, 'r', 'LineWidth', 2, 'DisplayName', 'Mean Cycle');

    % Add labels, title, and legend to the plot
    title('Extracted Cardiac Cycles and Mean Waveform', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Sample Points (Relative to R Wave)', 'FontSize', 12);
    ylabel('Amplitude (mV)', 'FontSize', 12);

    % Set x-axis limits to show only the relevant part of the wave
    xlim([1, preRspike + postRspike + 1]);

    % Improve visualization with grid and light background
    grid on;
    set(gca, 'Color', [0.95, 0.95, 0.95]);

    % Finalize the plot
    hold off;
end

%% Feature Extraction

windowSize = 4; % Number of samples in the sliding window
rangeP=25; % How far back is the P wave expected from Q wave
rangeQ=15; % How far back is the Q wave expected from R spike
rangeS=10; % How far forward is the S wave expected after R spike
rangeT=45; % How far forward is the T wave expected after S spike

% Variable Preallocation
prInterval=zeros(1,numSpikes-3);
rrInterval = zeros(1, numSpikes - 3);
qtInterval = zeros(1, numSpikes - 3);
stSegment = zeros(1, numSpikes - 3);
prSegment = zeros(1, numSpikes - 3);
pAmplitude = zeros(1, numSpikes - 3);
rAmplitude = zeros(1, numSpikes - 3);
tAmplitude = zeros(1, numSpikes - 3);
qAmplitude = zeros(1, numSpikes - 3);
sAmplitude = zeros(1, numSpikes - 3);


pWidth=zeros(1,numSpikes-3);
tWidth = zeros(1, numSpikes - 3);
qWidth = zeros(1, numSpikes - 3);
sWidth = zeros(1, numSpikes - 3);
qrsWidth = zeros(1, numSpikes - 3);

RRdist = diff(RpeakX); % Computes differences between consecutive R-peak indices
for i=2:numSpikes-1  
%for i=2:20 %spikes 2-10 for plotting
    idxStart=RpeakX(i)-preRspike; % start index of extracted wave
    idxEnd=RpeakX(i)+postRspike; % end index of extracted wave
    wave(i,:)=dataFilt(idxStart:idxEnd);% Extracts ith spike   
    firstDerivative = diff(wave(i, :)); % First derivative
    center=preRspike+1; % Index of R spike
    
    % Q wave Location
    [Qy,indQ]=min(wave(i,center-rangeQ:center)); % find min=Q
    indQ=indQ+preRspike-rangeQ; % Fix index in relation to 0
    Qwindow = indQ-windowSize-1:indQ-1; % Define Qwindow as a backward range from indQ
    Qsegment = wave(i, Qwindow);   % Extract corresponding segment
    QamplitudeDiffs = abs(Qsegment - wave(i, indQ)); % Difference from Q wave amplitude
    % Find the point with the minimum difference
    slidingMean = movmean(QamplitudeDiffs, 3); % Adjust window size if needed
    [~, minIndex] = min(diff(slidingMean)); % Find where differences minimize
    % Adjust minIndex to the original Qwindow
    qStartIndex = Qwindow(minIndex); % Map local index to global index
    
    % P wave Location
    [heightP,locP,widthP,promP]=findpeaks(wave(i,qStartIndex-rangeP:qStartIndex), ...
        'MinPeakDistance',length(wave(i,qStartIndex-rangeP:qStartIndex))-2); % findpeak P
    locP=locP+qStartIndex-rangeP-1; % Fix index in relation to 0 (Peak of P)
    startPInd=locP-round(widthP); % Index at which P wave starts

    % S wave Location
    [Sy,indS]=min(wave(i,center:rangeS+center)); % fin min=S
    indS=indS+center-1; % Fix index in relation to 0
    % Define Swindow as a forward range from indS
    Swindow = indS+1:indS + windowSize+1; % Ensure valid forward range thats not indS
    Ssegment = wave(i, Swindow);      % Extract corresponding segment
    % Calculate absolute differences relative to S wave amplitude
    SamplitudeDiffs = abs(Ssegment - wave(i, indS)); % Difference from S wave amplitude
    % Apply a moving average to smooth differences
    slidingMean = movmean(SamplitudeDiffs, 3); % Smooth with a window size of 3
    % Find the point with the minimum difference
    [~, minIndex] = min(diff(slidingMean)); % Find where differences minimize
    % Adjust minIndex to the original Swindow
    sEndIndex = Swindow(minIndex); % Map local index to global index
    
    % T wave Location
    [heightT,locT,widthT,promT]=findpeaks(wave(i,sEndIndex:sEndIndex+rangeT), ...
        'MinPeakDistance',length(wave(i,sEndIndex:sEndIndex+rangeT))-2); % findpeak T
    locT=locT+sEndIndex-1; % Fix index in relation to 0 (Peak of T)
   
    %% Storing Variables

    % Compute Intervals
    prInterval(i) = (qStartIndex - startPInd) / fs * 1000; % PR Interval in ms
    qtInterval(i) = (locT + round(widthT) - qStartIndex) / fs * 1000; % QT Interval in ms
    rrInterval=1000*(RRdist/fs);
    % Compute Segments
    stSegment(i) = (locT - round(widthT) - sEndIndex) / fs * 1000; % ST Segment in ms
    prSegment(i) = (qStartIndex - locP - round(widthP)) / fs * 1000; % PR Interval in ms

    % Compute Amplitudes
    pAmplitude(i) = promP; % P Wave Amplitude
    rAmplitude(i) = wave(i, center); % R Wave Amplitude
    tAmplitude(i) = promT; % T Wave Amplitude
    qAmplitude(i) = abs(Qy-wave(i,qStartIndex));
    sAmplitude(i) = abs(Sy-wave(i,sEndIndex));
    % Compute Wave Widths
    qWidth(i) = 2000*(indQ - qStartIndex) / fs; % Q wave width in ms
    sWidth(i) = 2000*(sEndIndex - indS) / fs; % S Wave Width in ms
    pWidth(i) = 2*widthP / fs * 1000; % P Wave Width in ms  
    qrsWidth(i) = (indS-indQ) / fs * 1000; % QRS Complex Width in ms

end

%% Diagnostics

% Summary Metrics
meanPR = mean(prInterval, 'omitnan');
meanRR = mean(rrInterval, 'omitnan');
meanQT = mean(qtInterval, 'omitnan');
meanSTs = mean(stSegment, 'omitnan');
meanPRs = mean(prSegment, 'omitnan');
meanQRS = mean(qrsWidth,  'omitnan');
meanpAmp = mean(pAmplitude, 'omitnan');
meanrAmp = mean(rAmplitude, 'omitnan');
meantAmp = mean(tAmplitude, 'omitnan');
meanqAmp = mean(qAmplitude, 'omitnan');
meansAmp = mean(sAmplitude, 'omitnan');
meanpWidth = mean(pWidth, 'omitnan');
meanqWidth = mean(qWidth, 'omitnan'); if meanqWidth==0;meanqWidth=NaN; end
meansWidth = mean(sWidth, 'omitnan'); if meansWidth==0;meansWidth=NaN; end
heartRateMean = mean(60*(RRdist/fs).^-1,'omitnan');

%CV calculation
muRR = mean(RRdist);
sigmaRR = std(RRdist);
CV = (sigmaRR / muRR) * 100;

disp(' ');
disp('Analyze completed! Now displaying results:');
disp(' ');
 
% Initialize the diagnostic table
Metrics = {
    'PR Interval (ms)';          % PR Interval in ms
    'RR Interval (ms)';          % RR Interval in ms
    'QT Interval (ms)';          % QT Interval in ms
    'QRS Duration (ms)';         % QRS Complex Duration in ms
    'ST Segment (ms)';           % ST Segment in ms
    'P Wave Amplitude (mV)';     % P Wave Amplitude in mV
    'R Wave Amplitude (mV)';     % R Wave Amplitude in mV
    'T Wave Amplitude (mV)';     % T Wave Amplitude in mV
    'Q Wave Depth (mV)';              %Q Wave Depth
    'S Wave Depth (mV)';              %S Wave Depth
    'Heart Rate (bpm)';          % Heart Rate in bpm
    'P Wave Duration (ms)';      % P Wave Width in ms
    'Q Wave Duration (ms)';      % Q Wave Width in ms
    'S Wave Duration (ms)';      % S Wave Width in ms
    'Coefficent of Variation (%)'
};

Recorded = [
    meanPR;          % Mean PR Interval (ms)
    meanRR;          % Mean RR Interval (ms)
    meanQT;          % Mean QT Interval (ms)
    meanQRS;         % Mean QRS Duration (ms)
    meanSTs;         % Mean ST Segment (ms)
    meanpAmp;        % Mean P Wave Amplitude (mV)
    meanrAmp;        % Mean R Wave Amplitude (mV)
    meantAmp;        % Mean T Wave Amplitude (mV)
    meanqAmp;        % Mean Q Wave Amplitude (mV)
    meansAmp;        % Mean S Wave Amplitude (mV)
    heartRateMean;   % Mean Heart Rate (bpm)
    meanpWidth;      % Mean P Wave Width (ms)
    meanqWidth;      % Mean Q Wave Width (ms)
    meansWidth;      % Mean S Wave Width (ms)
    CV               %Coeficient of Variation (%)
];

% Assign healthy ranges for each metric
HealthyRange = [
    120, 200;   % PR Interval: 120-200 ms
    600, 1200;  % RR Interval: 600-1200 ms (50-100 bpm)
    0, 400;      % QT Interval: <400 ms
    80, 120;    % QRS Duration: 80-120 ms
    0, 120;     % ST Segment: 0-120 ms
    0.1, 0.3;   % P Wave Amplitude: 0.1-0.3 mV
    0.5, 1.5;   % R Wave Amplitude: 0.5-1.5 mV
    0, 0.5;   % T Wave Amplitude: <0.5 mV
    0, 0.2;   % Q Wave Depth: <0.2 mV
    0, 1.5;   % S Wave Depth: <0.5 mV
    60, 100;    % Heart Rate: 60-100 bpm
    0, 120;    % P Wave Width: 60-110 ms
    0, 40;     % Q Wave Width: <40 ms
    30, 50      % S Wave Width: 30-50 ms
    0,10        % Standard CV range
];

% Sensitivity Threshold for Significant Deviation
sensitivityThreshold = 1; % Multiply by 1 for more sensitivity

% Calculate the range of each metric
range = HealthyRange(:, 2) - HealthyRange(:, 1);
% Default Diagnosis (Normal by default)
Diagnosis = repmat({'Normal'}, size(Metrics));

% Combine abnormalities into diagnostic categories
overallDiagnosis = '';
if heartRateMean > 100 && CV > 10
    overallDiagnosis = 'Tachycardia with Arrhythmia';
elseif ismember('STEMI', Diagnosis)
    overallDiagnosis = 'Acute Myocardial Infarction';
elseif ismember('Bundle Branch Block', Diagnosis)
    overallDiagnosis = 'Bundle Branch Block';
elseif ismember('1st Degree AV Block', Diagnosis)
    overallDiagnosis = '1st Degree AV Block';
elseif meanPR > 280
    overallDiagnosis = '2nd Degree AV Block or Complete AV Block';
elseif meanPR < 120
    overallDiagnosis = 'Wolff-Parkinson-White Syndrome';
elseif meanQRS > 120
    overallDiagnosis = 'Bundle Branch Block';
elseif meanSTs > 120
    overallDiagnosis = 'Possible Myocardial Ischemia or Injury';
elseif heartRateMean < 60
    overallDiagnosis = 'Bradycardia';
elseif heartRateMean > 100
    overallDiagnosis = 'Tachycardia';
elseif CV > 15
    overallDiagnosis = 'Highly Irregular Heart Rate (Arrhythmia)';
elseif CV > 10
    overallDiagnosis = 'Borderline Irregular Heart Rate';
else
    overallDiagnosis = 'Normal';
end

% Calculate deviations and normalize by range
lowerDeviation = ((HealthyRange(:, 1) - Recorded) ./ range) .* (HealthyRange(:, 1) > Recorded); % Below lower bound
upperDeviation = ((Recorded - HealthyRange(:, 2)) ./ range) .* (HealthyRange(:, 2) < Recorded); % Above upper bound
deviations = max(lowerDeviation, upperDeviation); % Take the larger normalized deviation

% Update Most Abnormal Metric
[maxDeviation, idx] = max(deviations);
mostAbnormalMetric = Metrics{idx};
mostAbnormalValue = Recorded(idx);
mostAbnormalRange = HealthyRange(idx,:);
mostAbnormalDiagnosis = Diagnosis{idx};

% Compare each metric against the healthy range and assign adaptive diagnoses
for i = 1:length(Metrics)
    if Recorded(i) < HealthyRange(i, 1)
        if deviations(i) <= sensitivityThreshold
                Diagnosis{i} = 'Abnormal Benign (Low)';
        else
                Diagnosis{i} = 'Abnormal Pathological (Low)';
        end
    elseif Recorded(i) > HealthyRange(i, 2)
        if deviations(i) <= sensitivityThreshold
                Diagnosis{i} = 'Abnormal Benign (High)';
        else
                Diagnosis{i} = 'Abnormal Pathological (High)';
        end
     end
end

% Create a diagnostic table
DiagnosticTable = table(Metrics, HealthyRange, Recorded, Diagnosis);

% Display the table
disp(DiagnosticTable);
% Display Overall Diagnosis
disp(['Most Abnormal Metric: ', mostAbnormalMetric]);
disp(['Classification: ', mostAbnormalDiagnosis]);
disp(['Overall Diagnosis: ', overallDiagnosis]);
disp(' ');

%% Explanation of Metrics

% Prompt the user to decide if they want to see the explanation
userInput = input('Display further details regarding metrics? (y/n): ', 's');

% Create a figure for visualizing individual cardiac cycles if the user chooses 'y'
if lower(userInput) == 'y'

% Display abnormality interpretation
disp(' ');
disp(' ');
disp(' ');
disp('Abnormalities in any of the following metrics may indicate the following:');
disp(' ');
disp('- **PR Interval**: Abnormalities may indicate heart block (delayed electrical signal).');
disp('- **RR Interval**: Abnormalities may indicate arrhythmias (irregular heart rhythms).');
disp('- **QT Interval**: Prolonged or shortened QT interval can indicate an increased risk of arrhythmias or sudden cardiac death.');
disp('- **QRS Duration**: Prolonged QRS duration can suggest conduction problems, including bundle branch blocks.');
disp('- **ST Segment**: Elevation or depression of the ST segment can be a sign of myocardial ischemia or heart attack.');
disp(['- **P Wave A' ...
    'mplitude**: Abnormal amplitude may indicate atrial enlargement or electrical conduction problems.']);
disp('- **R Wave Amplitude**: Low or high amplitude may suggest ventricular issues or enlargement.');
disp('- **T Wave Amplitude**: Abnormal T wave amplitude can indicate issues with ventricular recovery or ischemia.');
disp('- **Heart Rate**: Tachycardia (high heart rate) or bradycardia (low heart rate) may indicate various heart conditions.');
disp('- **P Wave Width**: Prolonged P wave width could be a sign of atrial enlargement or conduction delays.');
disp('- **T Wave Width**: Abnormalities may indicate issues with ventricular recovery.');
disp('- **Q Wave Width**: Wide or deep Q waves can indicate a previous heart attack or heart damage.');
disp('- **S Wave Width**: Prolonged S wave width may indicate ventricular conduction problems.');
disp(' ');

end

%% Repeat process per user request
% Ask the user if they want to analyze another file or end the session
userChoice = input('Would you like to analyze another ECG file/Graph Number? (y/n): ', 's');

% Check the user's input and proceed accordingly
if lower(userChoice) == 'y'
    % Clear variables and figures to start over
    clear; clc; close all;
    disp('Restarting the analysis...');

    % Call the script or function to start analyzing again (this assumes you are working in a script)
    run('Sandbox_Final.m'); % Replace with whatever the code's name is
else
    disp('Ending session. Goodbye!');
end
