clear all
close all
clc

% number of time points in the source space signal, here 180 s * 100 Hz
n = 180 * 100;

% rayleigh parameter
b = .5;

% number of resamples per coherence increment
nResamples = 1000;

% number of coherence increments
nCoherenceIncrements = 2000;

% define coherence values, negative values are used to calculate negative
% correlations, the actual coherence value's magnitude is used
coherenceValues = linspace(-1,1,nCoherenceIncrements);

% initialize the plain (RX) and orthogonalized (RO) correlation arrays
% containing the correlation value for each resample at each coherence
% increment
RX = NaN(nResamples,nCoherenceIncrements);
RO = NaN(nResamples,nCoherenceIncrements);

for iCoherenceIncrement=1:nCoherenceIncrements
   
   coherence = coherenceValues(iCoherenceIncrement);
disp(coherence)   
   for iResamples=1:nResamples

      % define random signal (sr)
      uniformRandomAmplitude = rand(1,n);
      randomPhase = exp(1i * 2 * pi * rand(1,n));

      sr = uniformRandomAmplitude .* randomPhase;
      sr = sr / norm(sr);

      % define signal 1 (s1)
      rayleighRandomAmplitude = raylrnd(b, [1,n]);
      randomPhase = exp(1i * 2 * pi * rand(1,n));
      
      s1 = rayleighRandomAmplitude .* randomPhase;
      s1 = s1 / norm(s1);
      
      if coherence < 0
         
         % define signal 1 (s1) to map from the inverse of the cdf for
         % negative correlations
         rayleighRandomAmplitudeMirror = raylinv(1 - raylcdf(rayleighRandomAmplitude,b),b);
         
         s1Mirror = rayleighRandomAmplitudeMirror .* randomPhase;
         s1Mirror = s1Mirror / norm(s1Mirror);
         
         % define signal 2 (s2)
         s2 = abs(coherence) * s1Mirror + sqrt(1 - coherence^2) * sr;
         s2 = s2 / norm(s2);
         
      else
         
         % define signal 2 (s2)
         s2 = abs(coherence) * s1 + sqrt(1 - coherence^2) * sr;
         s2 = s2 / norm(s2);

      end
      
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
      RX(iResamples,iCoherenceIncrement) = corr(log(PEplain1'),log(PEplain2'));

      % calculate the correlation of the logarithm (to render more normal)
      % of the orthogonalized power envelopes (RO) and take the mean of
      % both values to account for directionality
      RO(iResamples,iCoherenceIncrement) = (corr(log(PEplain1'),log(PEortho2')) + corr(log(PEplain2'),log(PEortho1'))) / 2;
            
   end
   
end

% calculate the best robust fit trendline to the data
robustTrendline = robustfit(RX(:),RO(:));
robustTrendlineIntercept = robustTrendline(1);
robustTrendlineSlope = robustTrendline(2);

% plot results
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
xlabel('Plain Power Envelope Correlation (\rho)')
ylabel('Orthogonalized Power Envelope Correlation (\rho)')
ax.LineWidth = 2;
axis([-1 1 -1 1])
axis('square')
hold on

% apply gridlines on the major axis value 0
gridline = linspace(-1,1,nCoherenceIncrements);
plot(zeros(size(gridline)),gridline,'k:','LineWidth',2)
plot(gridline,zeros(size(gridline)),'k:','LineWidth',2)

% plot correlation pairs
plot(RX(:),RO(:),'b.')

% plot robust trendline
robustTrendlineX = gridline;
robustTrendlineY = robustTrendlineX .* robustTrendlineSlope + robustTrendlineIntercept;
plot(robustTrendlineX,robustTrendlineY,'r','LineWidth',3)

% display trendline slope
t1 = text;
t1.String = sprintf('slope = %0.6f',robustTrendlineSlope);
t1.Position = [.3 -.5];
t1.FontName = 'Calibri';
t1.FontWeight = 'bold';
t1.FontSize = 16;