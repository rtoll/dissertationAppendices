clear all
close all
clc

% number of time points in the signal, here 60 s * 100 Hz
n = 60 * 100;

% rayleigh parameter
b = .5;

% define coherence to enforce a desired "true" correlation of 0.5.
% Correlation increases monotonically with increasing coherence
coherence = .8;

% number of iterations
nIterations = 100;

% number of resamples per iteration
nResamples = 1000;

% initialize figure
fig = figure;
fig.Position = get(0,'ScreenSize');
fig.Color = 'w';
ax = gca;
ax.Color = 'w';
ax.XColor = 'k';
ax.YColor = 'k';
ax.FontName = 'Calibri';
ax.FontWeight = 'bold';
ax.FontSize = 16;
xlabel('Correlation')
ylabel('Probability')
ax.LineWidth = 2;
hold on

for iIteration=1:nIterations

   % initialize the plain (RX) and orthogonalized (RO) correlation arrays
   % containing the correlation value for each resample
   RX = NaN(nResamples,1);
   RO = NaN(nResamples,1);

   for iResamples=1:nResamples
      
      % define signal 1 (s1)
      rayleighRandomAmplitude = raylrnd(b, [1,n]);
      randomPhase = exp(1i * 2 * pi * rand(1,n));
      
      s1 = rayleighRandomAmplitude .* randomPhase;
      s1 = s1 / norm(s1);
      
      % define random signal (sr)
      uniformRandomAmplitude = rand(1,n);
      randomPhase = exp(1i * 2 * pi * rand(1,n));
      
      sr = uniformRandomAmplitude .* randomPhase;
      sr = sr / norm(sr);
      
      % define signal 2 (s2)
      s2 = coherence * s1 + sqrt(1 - coherence^2) * sr;
      s2 = s2 / norm(s2);
      
      % preserve the relationship in amplitude between s1 and s2 but
      % randomize phase
      randomPhase = exp(1i * 2 * pi * rand(1,n));
      s1 = abs(s1) .* randomPhase;
      
      randomPhase = exp(1i * 2 * pi * rand(1,n));
      s2 = abs(s2) .* randomPhase;
      
      % calculate the plain power envelopes
      PEplain1 = s1 .* conj(s1);
      PEplain2 = s2 .* conj(s2);
            
      % calculate the power envelope of s1 orthogonalized with respect to
      % s2
      PEortho1 = (abs(imag(s1 .* (conj(s2) ./ abs(s2)))) .^ 2);
      
      % calculate the power envelope of s2 orthogonalized with respect to
      % s1
      PEortho2 = (abs(imag(s2 .* (conj(s1) ./ abs(s1)))) .^ 2);
      
      % calculate the correlation of the logarithm (to render more normal)
      % of the plain power envelopes (RX)
      RX(iResamples,1) = corr(log(PEplain1'),log(PEplain2'));
      
      % calculate the correlation of the logarithm (to render more normal)
      % of the orthogonalized power envelopes (RO) and take the mean of
      % both values to account for directionality
      RO(iResamples,1) = (corr(log(PEplain1'),log(PEortho2')) + corr(log(PEplain2'),log(PEortho1'))) / 2;
      
   end

   % calculate the probability densities of the correlation values in RX and
   % RO
   [probPlain,probPlainX] = ksdensity(RX);
   [probOrtho,probOrthoX] = ksdensity(RO);

   probPlain = probPlain / sum(probPlain);
   probOrtho = probOrtho / sum(probOrtho);

   % plot the probability densities for this iteration
   plot(probPlainX,probPlain,'b','LineWidth',2)
   plot(probOrthoX,probOrtho,'r','LineWidth',2)

end

% apply legend to figure
legend({'plain','orthogonalized'},'Location','best')