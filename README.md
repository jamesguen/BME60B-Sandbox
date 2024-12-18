**README: ECG Signal Analysis MATLAB Script**
Overview
This MATLAB script is designed to analyze ECG data from PhysioNet databases (https://drive.google.com/drive/folders/1r4cgBXT63XsNcLTBYUu0NfStxiLZeT7A?usp=sharing) or any compatible .mat file. It is optimized for ECG signals sampled at 128 Hz but can handle other sample rates with appropriate adjustments to the script.
Analysis Criteria
Detect R-peaks using signal processing techniques.
Extract individual cardiac waveforms.
Calculate diagnostic features: wave amplitudes, durations, and time intervals.
Generate a comprehensive table of metrics with classifications of abnormalities.
Evaluate deviations from healthy ranges and classify findings as Benign or Pathological.

**Key Features**
Interactive Workflow:
User prompts for file input, graph selection, and optional visualizations.
Visualization:
Plots the filtered ECG signal with detected R-peaks.
Displays extracted cardiac cycles and overlays the mean waveform.
Efficiency:
Optimized for matrix operations for fast processing of large datasets.
Interpretability:
Generates a clear diagnostic table of recorded metrics, healthy ranges, and classifications.
Optionally explains the clinical significance of each metric for user interpretation.
Flexibility:
Allows users to analyze multiple ECG files or graphs without restarting the script.


**Signal Processing Pipeline**

Filtering and R-Wave Detection

Uses a high-pass FIR filter to remove low-frequency noise.
Detects R-peaks with the findpeaks() function.
Accounts for both positive and negative polarity ECG signals.
Wave Component Identification

P wave: Detected using findpeaks() within a restricted backward range.
Q wave: Located by finding the minimum value in a restricted range before the R-peak.
S wave: Located by finding the minimum value in a restricted range after the R-peak.
T wave: Detected using findpeaks() in a restricted forward range.
Feature Extraction

Calculates key diagnostic features using matrix operations and sliding windows for efficiency.
Extracts the following metrics:
Intervals: PR, QT, RR, ST Segments
Amplitudes: P, Q, R, S, T waves
Durations: QRS width, P wave width, S wave width
Abnormality Evaluation

Evaluates deviations of recorded metrics from predefined healthy ranges.
Normalizes deviations to account for the range of each metric.
Classifies findings as:
Normal
Abnormal Benign (High/Low)
Abnormal Pathological (High/Low)
Diagnostics

Summarizes the results in a diagnostic table.
Highlights the most abnormal metric.
Provides an overall diagnosis based on thresholds, such as:
Tachycardia with Arrhythmia
Bradycardia
Bundle Branch Block
Acute Myocardial Infarction
