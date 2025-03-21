#! /usr/bin/env python

from src.scripts.diameTR import main
import src.scripts.diameTR


# Find the highest two peaks. min_gap defines the minimum distance of this two peaks
# Best for GalCer tubules
def find_peak_strategy( array_1d, min_gap=80):

	array_1d_np = cp.asnumpy(array_1d)
	array_1d_np = array_1d_np - array_1d_np.min()
	peaks = find_peaks( array_1d_np )[0]

	sorted_peaks = peaks[np.argsort(array_1d_np[peaks])][::-1][:6]
	peak_values = array_1d_np[ sorted_peaks]

	center = len(array_1d_np) // 2
	center = len(array_1d_np) // 2

	small_index_peak = None
	big_index_peak = None

	for peak in sorted_peaks:
		if peak < center - min_gap // 2:
			if small_index_peak is None or abs(peak - (center - min_gap // 2)) < abs(small_index_peak - (center - min_gap // 2)):
				small_index_peak = peak
		elif peak > center + min_gap // 2:
			if big_index_peak is None or abs(peak - (center + min_gap // 2)) < abs(big_index_peak - (center + min_gap // 2)):
				big_index_peak = peak

	if small_index_peak is not None and big_index_peak is not None:
		return small_index_peak, big_index_peak
	else:
		return 0, 0

src.scripts.diameTR.find_peak_strategy = find_peak_strategy

if __name__ == '__main__':
	main()