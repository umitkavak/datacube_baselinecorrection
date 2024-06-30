
# FITS Data Cube Processing for SOFIA Observatory's GREAT Instrument

This repository contains a Python script for processing a FITS data cube from the SOFIA Observatory's GREAT instrument. The script performs baseline correction, removes horizontal striping artifacts, and generates a moment 0 map.

## Requirements

Make sure you have the following Python packages installed:

- `numpy`
- `astropy`
- `matplotlib`
- `scipy`

You can install these packages using pip:

```bash
pip install numpy astropy matplotlib scipy
```

## Usage

1. **Place your FITS file in the same directory as the script**.
2. **Update the `filename` variable in the script to the path of your FITS file**.
3. **Run the script**:

```bash
python fits_data_cube_processing.py
```

## Script Overview

### Load the FITS Data Cube

Reads the FITS file and extracts the data cube and header information.

```python
filename = 'path_to_your_fits_file.fits'
hdu = fits.open(filename)
data_cube = hdu[0].data
header = hdu[0].header
```

### Extract Velocity Axis Information

Calculates the velocity axis using the header information.

```python
crval3 = header.get('CRVAL3', 0)
cdelt3 = header.get('CDELT3', 1)
crpix3 = header.get('CRPIX3', 1)
n_chan = data_cube.shape[0]
vel_axis = crval3 + (np.arange(n_chan) - crpix3) * cdelt3
```

### Baseline Subtraction

Fits and subtracts a polynomial baseline from each spectrum, excluding the specified velocity range of the [OI] line.

```python
def subtract_baseline_protected(cube, vel_axis, protect_min_vel, protect_max_vel, order=5):
    n_chan, n_y, n_x = cube.shape
    baseline_subtracted_cube = np.zeros_like(cube)
    baseline_cube = np.zeros_like(cube)
    
    protect_mask = (vel_axis >= protect_min_vel) & (vel_axis <= protect_max_vel)
    
    for i in range(n_y):
        for j in range(n_x):
            spectrum = cube[:, i, j]
            masked_spectrum = np.ma.masked_array(spectrum, mask=protect_mask)
            coeffs = np.polyfit(vel_axis[~masked_spectrum.mask], masked_spectrum[~masked_spectrum.mask], order)
            baseline = np.polyval(coeffs, vel_axis)
            baseline_subtracted_cube[:, i, j] = spectrum - baseline
            baseline_cube[:, i, j] = baseline
            
    return baseline_subtracted_cube, baseline_cube

baseline_subtracted_data_cube, baseline_cube = subtract_baseline_protected(data_cube, vel_axis, protect_min_vel, protect_max_vel, order=5)
```

### Median Filtering

Applies a median filter to reduce horizontal stripes.

```python
filtered_data_cube = median_filter(baseline_subtracted_data_cube, size=(1, 5, 5))
```

### Moment 0 Map Generation

Integrates the filtered data cube over the velocity axis to create the moment 0 map.

```python
moment0_map = np.sum(filtered_data_cube, axis=0)
```

### Save Outputs

Saves the baseline-subtracted data cube and the moment 0 map as new FITS files.

```python
baseline_subtracted_filename = 'NGC7538_upGREAT_OI63_res15_grid3_baseline_filtered.fits'
fits.writeto(baseline_subtracted_filename, filtered_data_cube, header, overwrite=True)

moment0_filename = 'NGC7538_upGREAT_OI63_res15_grid3_moment0.fits'
hdu_moment0 = fits.PrimaryHDU(moment0_map, header=header)
hdu_moment0.writeto(moment0_filename, overwrite=True)
```

### Plotting the Moment 0 Map

Plots the moment 0 map.

```python
plt.figure(figsize=(8, 8))
plt.imshow(moment0_map, origin='lower', cmap='inferno')
plt.colorbar(label='Integrated Intensity (K km/s)')
plt.xlabel('RA')
plt.ylabel('DEC')
plt.title('Moment 0 Map')
plt.show()
```

### Plotting Spectra Before and After Baseline Subtraction

Plots spectra at random positions to visually inspect the correction.

```python
def plot_spectra(data_cube, baseline_cube, vel_axis, positions):
    for (y, x) in positions:
        spectrum_original = data_cube[:, y, x]
        baseline = baseline_cube[:, y, x]
        spectrum_subtracted = spectrum_original - baseline
        
        plt.figure(figsize=(12, 6))
        plt.plot(vel_axis, spectrum_original, label='Original Spectrum')
        plt.plot(vel_axis, baseline, label='Fitted Baseline')
        plt.plot(vel_axis, spectrum_subtracted, label='Baseline Subtracted Spectrum')
        plt.xlabel('Velocity (m/s)')
        plt.ylabel('Intensity')
        plt.title(f'Spectrum at position (y={y}, x={x})')
        plt.legend()
        plt.show()

n_y = header['NAXIS2']
n_x = header['NAXIS1']
num_positions = 5
np.random.seed(42)
random_positions = [(np.random.randint(0, n_y), np.random.randint(0, n_x)) for _ in range(num_positions)]
plot_spectra(data_cube, baseline_cube, vel_axis, random_positions)
```

## Files

- `fits_data_cube_processing.py`: The main script for processing the FITS data cube.
- `README.md`: This file, describing the repository and usage instructions.

## License

This project is licensed under the MIT License.

## Contact

For any questions or issues, please contact umitkavak34@gmail.com
