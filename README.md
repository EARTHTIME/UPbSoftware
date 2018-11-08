# Repository of software for ID-TIMS U-Pb geochronology

This space is meant to be a community-editable compilation of useful codes and software for ID-TIMS U-Pb geochronology, from reducing/interpreting/compiling mass spectrometer data to interpreting ages in geological context.  To add to or edit this document, email \<noahmc at ku dot edu\> or open a pull request.

## Table of Contents:

### <u>Mass Spectrometry:</u>

#### Method files
| Instrument | File | Author(s) | Description |
|----------| :---------:| :---------: | ----------|
| Phoenix62 |  ... | ... | ... |

#### Deadtime and standards
| Software | Source | Author(s) | Notes |
|----------| :---------:| :---------: | ----------|
| PU code for Pb standards | url | Samperton et al. |
| KU code for Pb standards | [GitHub](https://github.com/noahmclean/TIMSLAB/tree/master/PbStandardsAnalysis) | McLean |
| BSU deadtime spreadsheet | url | McLean |
| BGC code for Pb alpha, DT | [NBS982_DTalpha_DAT.m](src/NBS982_DTalpha_DAT.m) | Samperton, Keller | Adapted for Sector 54 by Keller<br> from Samperton src |
| BGC code for U alpha, DT | [U500_DTalpha_DAT.m](src/U500_DTalpha_DAT.m) | Samperton, Keller | Adapted for Sector 54 by Keller<br> from Samperton src  |


#### Other
| Software | Source | Author(s) | Notes |
|----------| :---------:| :---------: | ----------|
| Collector Efficiencies | url | Princeton? | |
| Average Pb blank IC | url | McLean | |
| How long to measure baselines? | [GitHub](https://www.noahmclean.org/baseline-times/) | McLean | |


### <u>Data Reduction/Uncertainty Propagation:</u>

| Software | Source | Author(s) | Notes |
|----------| :---------:| :---------: | ----------|
| Tripoli | [CIRDLES website](http://cirdles.org/projects/tripoli/) | Jim Bowring, Noah McLean | |
| Schmitz and Schoene 'UPbR' spreadsheet | [BSU website](https://earth.boisestate.edu/isotope/labshare/data-reduction-software/) | Mark Schmitz, Blair Schoene | |
| ET_Redux | [Github](https://github.com/CIRDLES/ET_Redux/releases) | Bowring(s), McLean | |



### <u>Age Interpretation:</u>

| Software | Source | Author(s) |
|----------| :---------:| ----------|
| Keller et al., (2018) GPL 8, 31-35<br> Estimate eruption age from mineral age spectrum | [GitHub Repository](https://github.com/brenhinkeller/BayeZirChron.c), <br>[Demo Jupyter Notebook](https://mybinder.org/v2/gh/brenhinkeller/BayeZirChron.c/master?filepath=julia%2Fdemo.ipynb) | Keller, Schoene, Samperton |
| Bayesian stratigraphic age model | url | Robin Trayler
| Chron.jl: A Bayesian Framework for Integrated <br>Eruption Age and Age-Depth Modelling. <br> [doi:10.17605/osf.io/TQX3F](https://doi.org/10.17605/osf.io/TQX3F) | [GitHub Repository](https://github.com/brenhinkeller/Chron.jl), <br>[Demo Jupyter Notebook](https://mybinder.org/v2/gh/brenhinkeller/Chron.jl/master?filepath=examples%2Fdemo.ipynb)| Brenhin Keller



