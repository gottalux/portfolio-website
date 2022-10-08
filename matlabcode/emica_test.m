% EM-ICA implementation from Max Welling's 2001 paper
%
% Alastair Turl, 2016-02-22

%% GENERATE SOURCES AND MIXTURES

clear; clc;
numSamples = 88200; % number of samples
t = 0:numSamples-1; % time
Fs = 44100; % sampling rate

f1 = 453; f2 = 323; f3 = 563; % frequencies (Hz)

sources = []; % source signal matrix
sources(1,:) = sin(2*pi*f1*(t/Fs)); % sine wave (subgauss)
sources(2,:) = sign(sin(2*pi*f2*(t/Fs))); % square wave (subgauss)
%sources(3,:) = sawtooth(2*pi*f3*(t/Fs)); % sawtooth wave (subgauss)

numSources = size(sources,1); % number of sources
Mi = 2; % number of gaussian components per source
numComp = numSources*Mi; % total number of gaussian components

maxIter = 7; % number of EM iterations to perform

%mixingMatrix = 0.5*rand(numSources); % random mixing matrix
%mixingMatrix = [0.15 0.2 0.3; 0.2 0.1 0.23; 0.25 0.3 0.12];
mixingMatrix = [0.15 0.8; 0.2 0.32];

% noise parameter (beta)
beta = 0.02;

% mixtures: As + n
mixtures = (mixingMatrix*sources) + normrnd(0,beta,numSources,numSamples); 
numMixtures = size(mixtures,1); % number of mixtures
origMixtures = mixtures;

%% INITIALISE PARAMETERS

A = rand(numSources); % estimated mixing matrix
betasq = 0.05; % isotropic noise parameter

logPimc = normrnd(1,0.05,numComp,1); % gaussian component mixing coefficients
logMu = normrnd(1,0.05,numComp,1); % gaussian component means
logSigmasq = normrnd(1,0.05,numComp,1); % gaussian component variances

% all statistics use gammasq, b and alpha (Welling paper, eq. 3.9):
logGammasq = zeros(numComp,1);
logB = zeros(numComp,numSamples);
logAlpha = zeros(numComp,numSamples);

%% CENTRING AND WHITENING

mixturemu = mean(mixtures,2); % calculate mean vector mu
% subtract mu from each term
mixtures = bsxfun(@minus, mixtures, mixturemu);
[E,D] = eig(cov(mixtures')); % eigenvalue decomposition
whiten = sqrt(D)\E'; % whitening matrix
mixtures =  whiten * mixtures; % whiten mixtures

%% E-STEP: CALCULATE THE STATISTICS

for iteration = 1:maxIter

u = A\mixtures; % estimated sources (at least from whitened data)
logbe2si2 = log(betasq)+logSigmasq;
logGammasq = logbe2si2-log(exp(logSigmasq)+betasq-exp(logbe2si2));

for zi = 1:numComp
    temp = ((1-betasq)/betasq)*u(ceil(zi/Mi),:)+(exp(logMu(zi)-logSigmasq(zi)));
    logB(zi,:) = logGammasq(zi)+log(temp);
end
    
logAlpha = zeros(numComp,numSamples);
for n = 1:numSamples
    for zi = 1:numComp
        var = (betasq/(1-betasq))+exp(logSigmasq(zi));
        var = real(var); % ensure var is real for gaussian below
        logAlpha(zi,n) = logPimc(zi)+loggausspdf(u(ceil(zi/Mi),n),exp(logMu(zi)),var);
    end
end

% now ensure, for each sample n, alphas for each source sum to 1
logAlphaSums = zeros(numSources,numSamples);
for si = 1:numSources
    logAlphaSums(si,:) = logsumexp(logAlpha([(si*Mi)-(Mi-1):si*Mi],:),1);
    for n = 1:numSamples
        logAlpha([(si*Mi)-(Mi-1):si*Mi],n) = logAlpha([(si*Mi)-(Mi-1):si*Mi],n)-logAlphaSums(si,n);
    end
end

%%

% first order statistics <si>_zi^n and <si>_n
logE_s = logB;

logE_sSum = zeros(numSources,numSamples);
for si = 1:numSources
    logE_sSum(si,:) = logsumexp(logAlpha([(si*Mi)-(Mi-1):si*Mi],:)+logB([(si*Mi)-(Mi-1):si*Mi],:),1);
end

% second order statistics <si^2>_zi^n and <si^2>_n
logE_s2 = zeros(numComp,numSamples);
for zi = 1:numComp
    logE_s2(zi,:) = log(exp(logGammasq(zi))+exp(2*logB(zi,:)));
end

logE_s2Sum = zeros(numSources,numSamples);
for si = 1:numSources
    logE_s2Sum(si,:) = logsumexp(logAlpha([(si*Mi)-(Mi-1):si*Mi],:)+logE_s2([(si*Mi)-(Mi-1):si*Mi],:),1);
end

%% M-STEP: UPDATE PARAMETER ESTIMATES

% update mixing coefficient (Welling paper, eq. 4.8)
logPimc = logsumexp(logAlpha,2)-log(numSamples);

% update gaussian component means (Welling paper, eq. 4.8)
logMu = logsumexp(logAlpha+logE_s,2)-logsumexp(logAlpha,2);

% update gaussian component variances (Welling paper, eq. 4.8)
logSigmasq = real(log(exp(logsumexp(logAlpha+logE_s2,2)-logsumexp(logAlpha,2))-exp(logMu)));

% update rotation matrix R (Welling paper, eq. 4.3)
V = real(exp(log(mixtures*transpose(exp(logE_sSum)))-log(numSamples)));
R = orth(V); % orthonormal basis for V
A = sqrt(1-betasq)*R; % scaled version of R

% update beta (Welling paper, eq. 4.4)
% a1 = trace(R*exp(logE_sSum)*transpose(mixtures))/(numSamples*numSources);
% a2 = sum(sum(exp(logE_s2Sum)))/(numSamples*numSources);
a1 = exp(log(trace(R*exp(logE_sSum)*transpose(mixtures)))-log(numSamples)-log(numSources));
a2 = exp(logsumexp(logsumexp(logE_s2Sum))-log(numSamples)-log(numSources));
a1 = real(a1); a2 = real(a2);
cubicRoots = roots([1 -a1 a2 -a1]); % calculate cubic roots
w = real(cubicRoots(find(abs(imag(cubicRoots))<10e-9))); % select real root
betasq = sqrt(1-w^2); % calculate new beta
betasq = beta^2;

end