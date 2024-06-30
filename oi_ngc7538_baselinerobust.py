import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter

# Load the FITS data cube
filename = '/Volumes/Elements/Projects/NGC7538_SOFIACII/Observations_AllFiles/upGREAT/NGC7538_upGREAT_OI63_res15_grid3.fits'
hdu = fits.open(filename)
data_cube = hdu[0].data
header = hdu[0].header

# Extract velocity axis information from the header with default values if keys are not found
crval3 = header.get('CRVAL3', 0)
cdelt3 = header.get('CDELT3', 1)
crpix3 = header.get('CRPIX3', 1)

# Calculate the velocity axis
n_chan = data_cube.shape[0]
vel_axis = crval3 + (np.arange(n_chan) - crpix3) * cdelt3

# Define the velocity range to protect
protect_min_vel = -75e3  # km/s to m/s
protect_max_vel = -45e3  # km/s to m/s

# Function to fit and subtract the baseline while protecting the specified velocity range
def subtract_baseline_protected(cube, vel_axis, protect_min_vel, protect_max_vel, order=3):
    """
    Subtract a polynomial baseline from each spectrum in the data cube,
    protecting a specified velocity range.
    
    Parameters:
    cube (numpy.ndarray): The data cube with shape (n_chan, n_y, n_x)
    vel_axis (numpy.ndarray): The velocity axis
    protect_min_vel (float): The lower bound of the velocity range to protect
    protect_max_vel (float): The upper bound of the velocity range to protect
    order (int): The order of the polynomial to fit to the baseline
    
    Returns:
    numpy.ndarray: The baseline-subtracted data cube
    numpy.ndarray: The baseline cube
    """
    n_chan, n_y, n_x = cube.shape
    baseline_subtracted_cube = np.zeros_like(cube)
    baseline_cube = np.zeros_like(cube)
    
    # Create a mask for the velocity range to protect
    protect_mask = (vel_axis >= protect_min_vel) & (vel_axis <= protect_max_vel)
    
    for i in range(n_y):
        for j in range(n_x):
            spectrum = cube[:, i, j]
            # Mask the protected velocity range
            masked_spectrum = np.ma.masked_array(spectrum, mask=protect_mask)
            # Fit a polynomial to the unmasked baseline
            coeffs = np.polyfit(vel_axis[~masked_spectrum.mask], masked_spectrum[~masked_spectrum.mask], order)
            baseline = np.polyval(coeffs, vel_axis)
            baseline_subtracted_cube[:, i, j] = spectrum - baseline
            baseline_cube[:, i, j] = baseline
            
    return baseline_subtracted_cube, baseline_cube

# Subtract the baseline
baseline_subtracted_data_cube, baseline_cube = subtract_baseline_protected(data_cube, vel_axis, protect_min_vel, protect_max_vel, order=4)

# Apply a median filter to reduce horizontal stripes
filtered_data_cube = median_filter(baseline_subtracted_data_cube, size=(1, 5, 5))

# Save the filtered data cube to the specified FITS file
baseline_subtracted_filename = '/Volumes/Elements/Projects/NGC7538_SOFIACII/Observations_AllFiles/upGREAT/NGC7538_upGREAT_OI63_res15_grid3_baseline_filtered.fits'
fits.writeto(baseline_subtracted_filename, filtered_data_cube, header, overwrite=True)

# Create the moment 0 map by integrating over the velocity axis
moment0_map = np.sum(filtered_data_cube, axis=0)

# Save the moment 0 map as a new FITS file
moment0_filename = '/Volumes/Elements/Projects/NGC7538_SOFIACII/Observations_AllFiles/upGREAT/NGC7538_upGREAT_OI63_res15_grid3_moment0.fits'
hdu_moment0 = fits.PrimaryHDU(moment0_map, header=header)
hdu_moment0.writeto(moment0_filename, overwrite=True)

# Plot the moment 0 map
plt.figure(figsize=(8, 8))
plt.imshow(moment0_map, origin='lower', cmap='inferno', vmin=0, vmax=150)
plt.colorbar(label='Integrated Intensity (K km/s)')
plt.xlabel('RA')
plt.ylabel('DEC')
plt.title('Moment 0 Map')
plt.savefig("newmoment0_oi63.png")
#plt.show()

# Plot spectra before and after baseline subtraction
def plot_spectra(data_cube, baseline_cube, vel_axis, positions):
    """
    Plot spectra at specified positions before and after baseline subtraction.
    
    Parameters:
    data_cube (numpy.ndarray): The original data cube
    baseline_cube (numpy.ndarray): The baseline cube
    vel_axis (numpy.ndarray): The velocity axis
    positions (list of tuples): List of (y, x) positions to plot
    """
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

# Get pixel dimensions from the header
n_y = header['NAXIS2']
n_x = header['NAXIS1']

# Generate random positions
num_positions = 10
np.random.seed(42)  # For reproducibility
random_positions = [(np.random.randint(0, n_y), np.random.randint(0, n_x)) for _ in range(num_positions)]

# Plot spectra at the defined positions
plot_spectra(data_cube, baseline_cube, vel_axis, random_positions)
