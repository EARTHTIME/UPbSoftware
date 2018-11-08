% U500_DTalpha_SingleAnalysis_DAT.m
% Created on 2014-04-29, TUE, 19:30, by K.M. Samperton, Department of Geosciences, Princeton University, Princeton NJ USA
% Updated and adapted for Sector 54 .DAT files on 2017-02-27 by C.B. Keller, Berkeley Geochronology Center, Berkeley CA USA

% GOALS:  1) Reduce raw ID-TIMS U500 intensity data and calculate U
%            deadtime and fractionation values for Daly PM measurements.
%         2) Automate 2-sigma outlier exclusion of raw data, and 
%            fractionation and deadtime calculation of screened data.
%         3) Visualize fractionation and deadtime effects on U500 data.

% NOTES:  1) Variable nomenclature is as follows: 238U/235U = r85.
%         2) No oxide correction is applied as the interference at mass 270 
%            is negligible (18O-18O-234U on 16O-16O-238U).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PART 1. U MASS-DEPENDENT FRACTIONATION CORRECTION
    tic; % Time script

% STEP 1. Import the raw data from Sector 54 .dat files

% % Check we're in the right folder
%     cd ~/Desktop/'Calvin U500'/

% Enumerate files to import
    files = dir('*.DAT');

    i = 2;

    system(['grep -e ''^\( .[0-9]\.[0-9]*E.[0-9]*\)\{4\}'' ' files(i).name ' > temp.csv']);
    data = load('temp.csv');
    % Matrix format:  [mV 270, %beamgrowth/min., 267/270, 266/270]

% STEP 2. Set the isotopic composition of U500 from Condon et al. (2010): 
    r85s = 0.999781; % known 238U/235U
    r84s = 95.9364;  % known 238U/234U

% STEP 3. Calculate ratio vectors from the raw intensity data:
    r85_raw = 1./data(:,3); % measured 238U/235U (raw)
    r84_raw = 1./data(:,4); % measured 238U/234U (raw)
    r54_raw = data(:,3)./data(:,4); % measured 235U/234U (raw)
  
% STEP 4. Calculate the U mass-dependent fractionation from the difference 
% of the known 238U/235U and measured 238U/235U values (units: %/amu):
    Ualpha = 100*((r85s - mean(r85_raw))/r85s)/3

% STEP 5. Correct the raw ratios and intensities for mass fractionation:
    r85_alpha = r85_raw/((-Ualpha/100)*3 + 1);
    r84_alpha = r84_raw/((-Ualpha/100)*4 + 1);
    r54_alpha = r54_raw/((-Ualpha/100)*1 + 1);

    i238_alpha = data(:,1)*6.424*10^(18-11-3); % Intensity of 238U (convert from mV)
    i235_alpha = i238_alpha./r85_alpha;
    i234_alpha = i238_alpha./r84_alpha;
    data_alpha = [i234_alpha i235_alpha i238_alpha];
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PART 2. U DEADTIME CORRECTION

% STEP 6. Select a deadtime range (units: nanoseconds) over which to 
% evaluate the raw intensity data:
    DT_range = 0:0.01:100;

% STEP 7. Correct the raw intensity data over the specified deadtime range:
    DT_residuals = zeros(numel(DT_range),1);
    for i = 1:numel(DT_range)
      i234_alpha_DT_test = i234_alpha./(1-(DT_range(i).*1e-9).*i234_alpha);
      i238_alpha_DT_test = i238_alpha./(1-(DT_range(i).*1e-9).*i238_alpha);
      r84_alpha_DT_test  = i238_alpha_DT_test./i234_alpha_DT_test;
      DT_residuals(i,1)  = abs(r84s - mean(r84_alpha_DT_test));
    end
  
% STEP 8. Find the best-fitting deadtime value based on residuals:
    [val,row] = min(abs(DT_residuals));
    Udeadtime = DT_range(row)
  
% STEP 9. Create a fractionation- and deadtime-corrected intensity matrix 
% and ratio vectors:
    i234_alpha_DT = i234_alpha./(1-(Udeadtime.*1e-9).*i234_alpha);
    i235_alpha_DT = i235_alpha./(1-(Udeadtime.*1e-9).*i235_alpha);
    i238_alpha_DT = i238_alpha./(1-(Udeadtime.*1e-9).*i238_alpha);
    data_alpha_DT = [i234_alpha_DT i235_alpha_DT i238_alpha_DT];

    r85_alpha_DT = data_alpha_DT(:,3)./data_alpha_DT(:,2); 
    r84_alpha_DT = data_alpha_DT(:,3)./data_alpha_DT(:,1); 
    r54_alpha_DT = data_alpha_DT(:,2)./data_alpha_DT(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PART 3. DATA VISUALIZATION

% STEP 10. Plot the data:
% Panel 1. 238U/234U vs. intensity:
    figure
    subplot(1,3,1)
    plot(data_alpha_DT(:,3), r84_raw, 'ok', ...
      'MarkerFaceColor', 'r', 'MarkerSize', 6); hold on
    plot(data_alpha_DT(:,3), r84_alpha_DT, 'ok', ...
      'MarkerFaceColor', 'b', 'MarkerSize', 6); hold on
    plot([0 .7*10^6], [r84s, r84s], 'm-')
    xlabel('^{238}U intensity (cps)', 'fontsize', 14)
    ylabel('238/234', 'fontsize', 14)
    set(gca, 'fontsize', 14)
    grid on
    title('^{238}U/^{234}U vs. intensity')
    legend('Raw data', 'Corrected data', 'U500 value', ...
      'Location', 'SouthWest')

% Panel 2. 235U/234U vs. intensity:
    subplot(1,3,2)
    plot(data_alpha_DT(:,3), r54_raw, 'ok', ...
      'MarkerFaceColor', 'r', 'MarkerSize', 6); hold on
    plot(data_alpha_DT(:,3), r54_alpha_DT, 'ok', ...
      'MarkerFaceColor', 'b', 'MarkerSize', 6); hold on
    plot([0 .7*10^6], [r84s/r85s, r84s/r85s], 'm-')
    xlabel('^{238}U intensity (cps)', 'fontsize', 14)
    ylabel('235/234', 'fontsize', 14)
    set(gca, 'fontsize', 14)
    grid on
    title('^{235}U/^{234}U vs. intensity')

% Panel 3. 238U/235U vs. intensity:
    subplot(1,3,3)
    plot(data_alpha_DT(:,3), r85_raw, 'ok', ...
      'MarkerFaceColor', 'r', 'MarkerSize', 6); hold on
    plot(data_alpha_DT(:,3), r85_alpha_DT, 'ok', ...
      'MarkerFaceColor', 'b', 'MarkerSize', 6); hold on
    plot([0 .7*10^6], [r85s, r85s], 'm-')
    xlabel('^{238}U intensity (cps)', 'fontsize', 14)
    ylabel('238/235', 'fontsize', 14)
    set(gca, 'fontsize', 14)
    grid on
    title('^{238}U/^{235}U vs. intensity')
    set(gcf, 'color', 'w');

    toc;