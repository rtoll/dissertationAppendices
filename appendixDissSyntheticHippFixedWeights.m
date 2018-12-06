clear all
close all
clc

% construct the imaging kernel (K) of 5 fixed weights projecting to 2 electrodes
K = [1,1,0,1.25,.75; 0,1,1,.75,1.25];

% number of time points in the source space signal, here 60 s * 100 Hz
n = 60 * 100;

% rayleigh parameter
b = .5;

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

      % initialize source space (SS)
      SS = complex(NaN(size(K,2),n));

      for v=1:size(K,2) % v = vertex of SS
         
         rayleighRandomAmplitude = raylrnd(b, [1,n]);
         randomPhase = exp(1i * 2 * pi * rand(1,n));

         % generate a vertex of SS by multiplying amplitude by phase
         SS(v,:) = rayleighRandomAmplitude .* randomPhase;
         % scale the vertex of SS by its norm
         SS(v,:) = SS(v,:) / norm(SS(v,:));
         
      end

      % calculate the signal detected by the EEG electrodes
      EEG = K * SS;

      % calculate the plain power envelopes of the electrodes
      PEplain1 = EEG(1,:) .* conj(EEG(1,:));
      PEplain2 = EEG(2,:) .* conj(EEG(2,:));

      % calculate the power envelope of electrode 2 orthogonalized with respect to
      % electrode 1
      PEortho2 = (abs(imag(EEG(2,:) .* (conj(EEG(1,:)) ./ abs(EEG(1,:))))) .^ 2);

      % calculate the correlation of the logarithm (to render more normal)
      % of the plain power envelopes (RX)
      RX(iResamples,1) = corr(log(PEplain1'),log(PEplain2'));

      % calculate the correlation of the logarithm (to render more normal)
      % of the orthogonalized power envelopes (RO)
      RO(iResamples,1) = corr(log(PEplain1'),log(PEortho2'));

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