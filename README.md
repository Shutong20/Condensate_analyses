# Condensate_analyses
These MATLAB codes are used in the paper of Shu et al. "Mesoscale molecular assembly is favored by the active, crowded cytoplasm".

Microscope images and time course movies of droplets are first processed through FIJI plugin 'Trackmate', which generate '.csv' files indicating spots properties using particle detection function and '.xlm' files indicating particle tracking trajectories. The codes then use these '.csv' and '.xlm' files as inputs and analyze droplets properties including the mean pixel intensity, total droplet intensity, number of droplets etc, as well as droplets trajectories properties including diffusivity, angle correlation function and normalized velocity autocorrelation functions.
