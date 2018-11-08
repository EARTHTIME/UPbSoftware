% NBS982_DTalpha_DAT.m
% Created on 2014-04-29, TUE, 19:30, by K.M. Samperton, Department of Geosciences, Princeton University, Princeton NJ USA
% Updated and adapted for Sector 54 .DAT files on 2017-02-27 by C.B. Keller, Berkeley Geochronology Center, Berkeley CA USA

% GOALS:  1) Reduce raw ID-TIMS NBS982 intensity data and calculate Pb
%            deadtime and fractionation values for Daly PM measurements.
%         2) For a single data file, automate 2-sigma outlier screening 
%            of raw data and fractionation/deadtime calculations.
%         3) Apply cycle-by-cycle fraction corrections to each file.
%         3) Visualize outputs of reducing a single NBS982 run.

% NOTES:  1) Variable nomenclature is as follows: 208Pb/204Pb = r84.
%         2) No corrections for BaPO4 or Tl isobaric interferences are 
%            applied. This may need to be addressed in future versions,
%            but is likely neglibile for large NBS982 loads.
%         3) This script is structured to extract information from a single
%            data matrix, called 'data', usually imported from an Excel 
%            file. The matrix should have the form of date+name (e.g., 
%            "2014-05-13 NBS982 20140105_1"), where the first 10 characters
%            are the date as yyyy-mm-dd and the following characters are 
%            the sample name preceded by a space.
%         4) 'data' should be in the form of a nx4 matrix, where columns 
%            1-4 are raw, non-deadtime corrected Pb intensities in order of
%            increasing mass ([i204 i206 i207 i208]), and rows are single 
%            data cycles. Importantly, this raw data should be corrected 
%            for beam interpolation but NOT corrected for fractionation 
%            and deadtime.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PART A. DATA IMPORT AND OUTLIER EXCLUSION

% STEP 1.  Time script
    tic;

% Step 2. Import the raw data from Sector 54 .dat files
  
% % Check we're in the right folder
%     cd ~/Desktop/'Calvin NBS982'/

% Enumerate files to import
    files = dir('*.DAT');

    i = 1;

    system(['grep -e ''^\( .[0-9]\.[0-9]*E.[0-9]*\)\{5\}'' ' files(i).name ' > temp.csv']);
    data = load('temp.csv');
    % Matrix format:  [mV 206, %beamgrowth/min., 206/204, 206/207, 206/208]


% Step 3. Set NBS982 isotopic composition from Condon et al. (2014): 
  r64s = 36.7569;   % known 206Pb/204Pb
  r76s = 0.466967;  % known 207Pb/206Pb
  r86s = 1.000249;  % known 208Pb/206Pb
  r84s = r64s*r86s; % known 208Pb/204Pb
  r87s = r86s/r76s; % known 208Pb/207Pb

% Step 4. Calculate ratio vectors from the raw intensity data:
  i206_raw = data(:,1)*6.424*10^(18-11-3); % Intensity of 206Pb (convert from mV)
  r64_raw = data(:,3); % measured 206Pb/204Pb (raw)
  r76_raw = 1./data(:,4); % measured 207Pb/206Pb (raw)clear
  r86_raw = 1./data(:,5); % measured 208Pb/206Pb (raw)
  r84_raw = r86_raw.*r64_raw; % measured 208Pb/204Pb (raw)
  r87_raw = r86_raw.*data(:,4); % measured 208Pb/207Pb (raw)
  
  
  % Step 5. Find the indices of the raw data within a 2-sigma envelope of 
% the mean 208Pb/206Pb. Extract these data and sort the screened data in 
% descending order to produce an outlier-screened, ordered data matrix:
  I = find(r86_raw > (mean(r86_raw) - 2*std(r86_raw)) & ...
           r86_raw < (mean(r86_raw) + 2*std(r86_raw)));
  data_screened = data(I,:);
  DATA = flipud(sortrows(data_screened));
  Cycles_TOTAL   = size(data,1) %#ok<NOPTS>
  Cycles_QC_PASS = size(DATA,1) %#ok<NOPTS>
  Cycles_QC_FAIL = Cycles_TOTAL - Cycles_QC_PASS %#ok<NOPTS>
    
% Step 6. Calculate ratio vectors from the screened intensity data:
  i206_raw_s = DATA(:,1)*6.424*10^(18-11-3); % Intensity of 206Pb (convert from mV)
  r64_raw_s = DATA(:,3); % measured 206Pb/204Pb (screened)
  r76_raw_s = 1./DATA(:,4); % measured 207Pb/206Pb (screened)
  r86_raw_s = 1./DATA(:,5); % measured 208Pb/206Pb (screened)
  r84_raw_s = r86_raw_s.*r64_raw_s; % measured 208Pb/204Pb (screened)
  r87_raw_s = r86_raw_s.*DATA(:,4); % measured 208Pb/207Pb (screened)
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PART B. Pb MASS-DEPENDENT FRACTIONATION CORRECTION

% Step 1. Calculate Pb mass-dependent fractionation from the difference 
% of the known 208Pb/206Pb and measured 208Pb/206Pb values (units: %/amu).
% Apply a cycle-by-cycle (i.e., "real time") fractionation correction to
% all data. Compile the fractionation values:
  data_alpha  = zeros(size(DATA,1),4);
  PbalphaCOMP = zeros(size(DATA,1),1);
  for k = 1:size(DATA,1)
      Pbalpha = 100*((r86s - r86_raw_s(k))/r86s)/2;
      r64_alpha = r64_raw_s(k)/((-Pbalpha/100)*2 + 1); 
      r76_alpha = r76_raw_s(k)/((-Pbalpha/100)*1 + 1);
      r86_alpha = r86_raw_s(k)/((-Pbalpha/100)*2 + 1);
      r84_alpha = r84_raw_s(k)/((-Pbalpha/100)*4 + 1);
      r87_alpha = r87_raw_s(k)/((-Pbalpha/100)*1 + 1);  
      data_alpha(k,1) = i206_raw_s(k)./r64_alpha;
      data_alpha(k,2) = i206_raw_s(k);
      data_alpha(k,3) = i206_raw_s(k).*r76_alpha;
      data_alpha(k,4) = i206_raw_s(k).*r86_alpha;
      PbalphaCOMP(k) = Pbalpha; % #ok<SAGROW>
  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PART C. Pb DEADTIME CORRECTION

% Step 1. Select a deadtime range (units: nanoseconds) over which to 
% evaluate the raw intensity data:
  DT_range = 0:0.01:100;
  
% Step 2. Correct the raw intensity data over the specified deadtime range:
  DT_residuals = zeros(numel(DT_range),1);
  p_fitvalues  = zeros(numel(DT_range),1);
  h_fitvalues  = zeros(numel(DT_range),1);
  r76_refarray = r76s*ones(size(data_alpha,1),1);
  for i = 1:numel(DT_range)
      i206_alpha_DT_test = data_alpha(:,2)./(1-(DT_range(i).*1e-9)...
          .*data_alpha(:,2));
      i207_alpha_DT_test = data_alpha(:,3)./(1-(DT_range(i).*1e-9)...
          .*data_alpha(:,3));
      r76_alpha_DT_test  = i207_alpha_DT_test./i206_alpha_DT_test;
      DT_residuals(i,1) = sum((abs((r76s*ones(size(data_alpha,1),1)) ...
          - r76_alpha_DT_test)).^2);
      [h,p,ci,stats] = ttest2(r76_refarray,r76_alpha_DT_test);
      h_fitvalues(i,1) = h;
      p_fitvalues(i,1) = p;
  end
  
  % Step 2b. Calculate confidence interval bounds from T-Test results:
  [Tval,loc] = max(p_fitvalues);
  bestTfitDT = DT_range(loc);
  lowerlimit = DT_range(find(h_fitvalues==0,1));
  upperlimit = DT_range(find(h_fitvalues==0,1,'last'));
  
% Step 3. Determine the best-fitting deadtime value based on residuals:
  [val,row] = min(abs(DT_residuals));
  Pbdeadtime = DT_range(row) %#ok<NOPTS>
  
% Step 4. Create a fractionation- and deadtime-corrected intensity matrix 
% and ratio vectors:
  i204_alpha_DT = data_alpha(:,1)./(1-(Pbdeadtime.*1e-9).*data_alpha(:,1));
  i206_alpha_DT = data_alpha(:,2)./(1-(Pbdeadtime.*1e-9).*data_alpha(:,2));
  i207_alpha_DT = data_alpha(:,3)./(1-(Pbdeadtime.*1e-9).*data_alpha(:,3));
  i208_alpha_DT = data_alpha(:,4)./(1-(Pbdeadtime.*1e-9).*data_alpha(:,4));
  data_alpha_DT = [i204_alpha_DT i206_alpha_DT ...
      i207_alpha_DT i208_alpha_DT];
  
  r64_alpha_DT = data_alpha_DT(:,2)./data_alpha_DT(:,1); 
  r76_alpha_DT = data_alpha_DT(:,3)./data_alpha_DT(:,2); 
  r86_alpha_DT = data_alpha_DT(:,4)./data_alpha_DT(:,2);
  r84_alpha_DT = data_alpha_DT(:,4)./data_alpha_DT(:,1);
  r87_alpha_DT = data_alpha_DT(:,4)./data_alpha_DT(:,3);
  
% Step 5. Print Pb fractionation mean and standard deviation:
  Pbalpha_mean  = mean(PbalphaCOMP) %#ok<NOPTS>
  Pbalpha_stdev = std(PbalphaCOMP) %#ok<NOPTS>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PART D. DATA VISUALIZATION AND PLOTTING

% Figure 1, Panel 1. 206Pb/204Pb vs. intensity:
  figure
  subplot(2,3,1)
  plot(data_alpha_DT(:,4), r64_raw_s, 'ok', ...
      'MarkerFaceColor', 'r', 'MarkerSize', 6); hold on
  plot(data_alpha_DT(:,4), r64_alpha_DT, 'ok', ...
      'MarkerFaceColor', 'b', 'MarkerSize', 6); hold on
  plot([0 max(data_alpha_DT(:,4))], [r64s, r64s], 'm-')
  xlabel('^{208}Pb intensity (cps)', 'fontsize', 14)
  ylabel('206/204', 'fontsize', 14)
  set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')
  xlim([0 max(data_alpha_DT(:,4))])
  set(gca, 'fontsize', 14)
  grid on
  title('^{206}Pb/^{204}Pb vs. intensity')

% Figure 1, Panel 2. 207Pb/206Pb vs. intensity:
  subplot(2,3,2)
  plot(data_alpha_DT(:,4), r76_raw_s, 'ok', ...
      'MarkerFaceColor', 'r', 'MarkerSize', 6); hold on
  plot(data_alpha_DT(:,4), r76_alpha_DT, 'ok', ...
      'MarkerFaceColor', 'b', 'MarkerSize', 6); hold on
  plot([0 max(data_alpha_DT(:,2))], [r76s, r76s], 'm-')
  xlabel('^{208}Pb intensity (cps)', 'fontsize', 14)
  ylabel('207/206', 'fontsize', 14)
  set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')
  xlim([0 max(data_alpha_DT(:,4))])
  set(gca, 'fontsize', 14)
  grid on
  title('^{207}Pb/^{206}Pb vs. intensity')

% Figure 1, Panel 3. 208Pb/206Pb vs. intensity:
  subplot(2,3,3)
  plot(data_alpha_DT(:,4), r86_raw_s, 'ok', ...
      'MarkerFaceColor', 'r', 'MarkerSize', 6); hold on
  plot(data_alpha_DT(:,4), r86_alpha_DT, 'ok', ...
      'MarkerFaceColor', 'b', 'MarkerSize', 6); hold on
  plot([0 max(data_alpha_DT(:,2))], [r86s, r86s], 'm-')
  xlabel('^{208}Pb intensity (cps)', 'fontsize', 14)
  ylabel('208/206', 'fontsize', 14)
  xlim([0 max(data_alpha_DT(:,4))])
  set(gca, 'fontsize', 14)
  title('^{208}Pb/^{206}Pb vs. intensity')
  grid on

% Figure 1, Panel 4. 208Pb/204Pb vs. intensity:
  subplot(2,3,4)
  plot(data_alpha_DT(:,4), r84_raw_s, 'ok', ...
      'MarkerFaceColor', 'r', 'MarkerSize', 6); hold on
  plot(data_alpha_DT(:,4), r84_alpha_DT, 'ok', ...
      'MarkerFaceColor', 'b', 'MarkerSize', 6); hold on
  plot([0 max(data_alpha_DT(:,2))], [r84s, r84s], 'm-')
  xlabel('^{208}Pb intensity (cps)', 'fontsize', 14)
  ylabel('208/204', 'fontsize', 14)
  set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')
  xlim([0 max(data_alpha_DT(:,4))])
  set(gca, 'fontsize', 14)
  title('^{208}Pb/^{204}Pb vs. intensity')
  grid on
  
  legend('Raw (screened)', '\alpha + DT corrected', ...
    'True NBS982 IC', 'Location', 'SouthWest')

% Figure 1, Panel 5. 208Pb/207Pb vs. intensity:
  subplot(2,3,5)
  plot(data_alpha_DT(:,4), r87_raw_s, 'ok', ...
      'MarkerFaceColor', 'r', 'MarkerSize', 6); hold on
  plot(data_alpha_DT(:,4), r87_alpha_DT, 'ok', ...
      'MarkerFaceColor', 'b', 'MarkerSize', 6); hold on
  plot([0 max(data_alpha_DT(:,2))], [r87s, r87s], 'm-')
  xlabel('^{208}Pb intensity (cps)', 'fontsize', 14)
  ylabel('208/207', 'fontsize', 14)
  set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')
  xlim([0 max(data_alpha_DT(:,4))])
  set(gca, 'fontsize', 14)
  title('^{208}Pb/^{207}Pb vs. intensity')
  grid on

% Figure 1, Panel 6. Weighting factor vs. intensity:
  subplot(2,3,6)
  plot(data_alpha_DT(:,4), PbalphaCOMP, '^k', ...
      'MarkerFaceColor', 'g', 'MarkerSize', 6); hold on
  plot([0 max(data_alpha_DT(:,4))], [Pbalpha_mean Pbalpha_mean], ...
      'LineWidth', 2, 'Color', 'b'); hold on
  plot([0 max(data_alpha_DT(:,4))], [Pbalpha_mean+Pbalpha_stdev ...
      Pbalpha_mean+Pbalpha_stdev], '--', 'LineWidth', 2, ...
      'Color', 'b'); hold on
  plot([0 max(data_alpha_DT(:,4))], [Pbalpha_mean-Pbalpha_stdev ...
      Pbalpha_mean-Pbalpha_stdev], '--', 'LineWidth', 2, ...
      'Color', 'b'); hold on
  xlabel('^{208}Pb intensity (cps)', 'fontsize', 14)
  ylabel('Pb \alpha (%/amu)', 'fontsize', 14)
  xlim([0 max(data_alpha_DT(:,4))])
  set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')
  set(gca, 'fontsize', 14)
  title('Pb fractionation vs. intensity (mean ± 1sig)')
  grid on
  set(gcf, 'color', 'w');
    
  toc
  
  