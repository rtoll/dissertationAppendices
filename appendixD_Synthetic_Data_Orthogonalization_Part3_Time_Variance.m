clear all
close all
clc

% construct the imaging kernel (K) of 5 fixed weights projecting to 2 electrodes
K = [1,1,0,1.25,.75; 0,1,1,.75,1.25];

% rayleigh parameter
b = .5;

% number of iterations
nIterations = 100;

% number of resamples per iteration
nResamples = 1000;

% set the number of time points to vary based on the iteration by defining
% an array of number of time points (N) to draw from, ranging from 1 second
% to 600 seconds (10 minutes), each multiplied by 100 Hz
N = floor(linspace(1 * 100,600 * 100,nIterations));

% initialize figures
fig1 = figure;
fig1.Position = get(0,'ScreenSize');
fig1.Color = 'w';
ax1 = gca;
ax1.Color = 'w';
ax1.XColor = 'k';
ax1.YColor = 'k';
ax1.FontName = 'Calibri';
ax1.FontWeight = 'bold';
ax1.FontSize = 16;
xlabel('Correlation')
ylabel('Probability')
ax1.LineWidth = 2;

fig2 = figure;
fig2.Position = get(0,'ScreenSize');
fig2.Color = 'w';
ax2 = gca;
ax2.Color = 'w';
ax2.XColor = 'k';
ax2.YColor = 'k';
ax2.FontName = 'Calibri';
ax2.FontWeight = 'bold';
ax2.FontSize = 16;
xlabel('Time (s)')
ylabel('Standard Deviation of Correlation Coefficients')
ax2.LineWidth = 2;

for iIteration=1:nIterations

   % set the number of time points to vary based on the iteration
   n = N(iIteration);
   
   % set the color of each line to vary based on the iteration
   plainLineColor = [0,1 - (iIteration / nIterations),1];
   orthoLineColor = [1,0,1 - (iIteration / nIterations)];

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

   set(0, 'currentfigure', fig1);
   set(fig1, 'currentaxes', ax1);
   hold on
   
   % plot the probability densities for this iteration
   plot(probPlainX,probPlain,'Color',plainLineColor,'LineWidth',2)
   plot(probOrthoX,probOrtho,'Color',orthoLineColor,'LineWidth',2)
   
   hold off
   
   set(0, 'currentfigure', fig2);
   set(fig2, 'currentaxes', ax2);
   hold on
   % alternate plot order for better visibility
   if mod(iIteration,2) == 0
      plot(n/100,std(RX),'Marker','o','MarkerSize',10,'MarkerEdgeColor',plainLineColor,'MarkerFaceColor',plainLineColor)
      plot(n/100,std(RO),'Marker','o','MarkerSize',10,'MarkerEdgeColor',orthoLineColor,'MarkerFaceColor',orthoLineColor)
   else
      plot(n/100,std(RO),'Marker','o','MarkerSize',10,'MarkerEdgeColor',orthoLineColor,'MarkerFaceColor',orthoLineColor)
      plot(n/100,std(RX),'Marker','o','MarkerSize',10,'MarkerEdgeColor',plainLineColor,'MarkerFaceColor',plainLineColor)
   end
   hold off

end