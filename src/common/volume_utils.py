import cupy as  cp
from cupyx.scipy.ndimage import rotate, shift
import copy

def rotate_volume( volume,rot=0,tilt=0,psi=0,order=1):
	newvolume = copy.deepcopy(volume)
	if rot != 0:
		newvolume = rotate(newvolume, angle=rot, axes=(1,2), order=order, mode='constant', reshape=False)
	if tilt != 0:
		newvolume = rotate(newvolume, angle=-tilt, axes=(0,2), order=order, mode='constant', reshape=False)
	if psi != 0:
		newvolume = rotate(newvolume, angle=psi, axes=(1,2), order=order, mode='constant', reshape=False)
	return newvolume


def rotate_image( image, psi=0, x=0, y=0, order=1 ):
	if x != 0 or y != 0:
		image = shift ( image, (y,x), mode='wrap' )
	if psi != 0:
		image = rotate ( image, -psi, axes=(0,1), mode='constant', reshape=False, order=order )
	return image



### for rotating volume projections
def rotate_image_3dto2d ( image, psi=0, x=0, y=0, order=1, rotate_mode='constant' ):
	if psi != 0:
		image = rotate ( image, psi, axes=(0,1), mode=rotate_mode, reshape=False, order=order )
	if x!=0 or y!= 0:
		image = shift ( image, (-y,-x), mode='constant' )
	return image


def project_volume( volume, rot=0, tilt=0, psi=0, x=0, y=0, order=1):
	volume_rotated = rotate_volume ( volume, rot, tilt, 0, order=order)
	projection = cp.sum( volume_rotated, axis=0 )
	projection = rotate_image_3dto2d ( projection, psi, x, y, order)
	return projection



def image_bin( img, bin_factor):
	while bin_factor != 1:
		img = img.reshape(img.shape[0]//2, 2, img.shape[1]//2, 2)
		img = cp.mean(img, axis=3)
		img = cp.mean(img, axis=1)
		bin_factor /= 2
	return img


def create_circular_mask(box_size, mask_radius):
	center = box_size // 2
	y, x = np.ogrid[:box_size, :box_size]

	# Calculate distances from center
	dist_from_center = np.sqrt((x - center)**2 + (y - center)**2)

	# Create the mask
	mask = (dist_from_center <= mask_radius).astype(int)

	return mask

def low_pass_filter( image_array, resolution=20, pixel_size=1):
	box_size = image_array.shape[0]
	mask_radius = box_size * pixel_size / resolution

	mask = create_circular_mask(box_size, mask_radius)
	image_fft = cp.fft.fft2(image_array)
	image_fft = cp.fft.fftshift(image_fft)

	result = image_fft * mask
	result = cp.fft.ifftshift(result)
	result = cp.fft.ifft2(result).real

	return result.astype(cp.float32)