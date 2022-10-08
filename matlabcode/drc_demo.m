% A demo of how to use simple DRC
%
% Alastair Turl / a.c.turl@cs.bham.ac.uk / 2016-03-20

clear; clc;

% simple_drc(filename, ratio, percentile, time)
%   filename = audio file to apply DRC to
%   ratio = compression ratio (1: no compression, Inf: hard clip)
%   percentile = (0-100: determine the threshold)
%   time = attack & release time in milliseconds

[audio, Fs] = simple_drc('tabla.aif',4,75,0);

% parameters
play = 0; % 1 to play the resulting DRC processed audio
plotHist = 1; % 1 to plot a histogram of the processed audio

if play == 1
    soundsc(audio,Fs);
end

if plotHist == 1
    figure(1);
    subplot(1,2,1)
    hist(audio(:,1),50);
    title(['DRC audio kurtosis (L): ' ...
        num2str(kurtosis(audio(:,1)))]);
    subplot(1,2,2)
    hist(audio(:,2),50);
    title(['DRC audio kurtosis (R): ' ...
        num2str(kurtosis(audio(:,2)))]);
end