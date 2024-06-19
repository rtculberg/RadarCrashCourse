import numpy as np 
import scipy

#-------------------------------------------------------------------------------------------------
# surf_ind = Surface2Index(surface, ax_label)
#
# Takes in:
# surface - 1D array of the location of the surface in the radargram in terms of two-way travel 
# time, elevation, etc
# ax_label - 1D array of the axis label for each vertical pixel in a trace (for example, two-way 
# travel time or elevation)
#
# Returns:
# surf_ind - 1D array of the vertical pixel index of the surface in the radargram
#-----------------------------------------------------------------------------------------------------
def Surface2Index(surface,ax_label):
    surf_ind = np.empty(surface.shape)
    for i in range(surf_ind.size):
        surf_ind[i] = np.argmin(np.abs(ax_label - surface[i]))
    return surf_ind.astype(int)

# Helper function to create a gaussian kernel
def gaussian(x,mu,sigma):
    return np.exp(-np.power((x-mu)/sigma, 2.)/2.)

# Helper function to find the first peak after the MacFerrin tracked surface, which tends to sit on the max gradient
def first_maximum(surface_indices, data):
    improved_surface = np.empty(surface_indices.shape, dtype=surface_indices.dtype)
    for i in range(improved_surface.shape[0]):
        start = surface_indices[i].astype(int)
        stop = start + 20
        if stop >= data.shape[1]:
            stop = data.shape[1]-1
        seg = data[start:stop,i]
        peaks,_ = scipy.signal.find_peaks(seg)
        if peaks.shape[0] == 0:
            improved_surface[i] = surface_indices[i]
        else:
            improved_surface[i] = peaks[0] + start
    return improved_surface

# Helper function to remove large jumps in surface tracking
def get_rid_of_false_surface_jumps(surface_indices):
    improved_surface = surface_indices.copy()

    jumps = improved_surface[1:] - improved_surface[:-1]
    # Substitute any large jumps with brightest pixel in a window of original surface.  Do this until large jumps either go away or have all been corrected to original surface.
    for i in range(len(jumps)):

        # Slope windowsize = number of pixels we use to average the previous slope.
        slope_windowsize = 10
        if i < slope_windowsize:
            continue
        mean_slope = np.mean(np.array(jumps[i-slope_windowsize:i], dtype=np.float64))

        # Find the difference of this slope from the last five stops
        difference_from_mean_slope = jumps[i] - mean_slope
        # Ignore if it's jumped less than 3 from the mean recent slope, or less than 50% greater than the mean slope at this time.
        if (difference_from_mean_slope < 5) or (difference_from_mean_slope < (1.5*mean_slope)):
            continue

        # tune settings
        jump_lookahead = 20 # Number of pixels to look ahead to see if we can find a matching down-jump
        if i+jump_lookahead > len(jumps):
            jump_lookahead = len(jumps) - i

        # This is how close the surface on the "other side" of the jump must be to the original slope to be considered for it.
        jump_magnitude_threshold = 1.10

        # See if we can find a point in the near future that would approximate the current slope.
        slopes_ahead = np.cumsum(jumps[i:i+jump_lookahead]) / np.arange(1,jump_lookahead+1)
        opposite_match = np.argmax(slopes_ahead <= (mean_slope * jump_magnitude_threshold))

        if opposite_match > 0:
            # We found a match, onward!
            opposite_match_index = i + opposite_match
            for j in range(i+1,opposite_match_index+1):
                improved_surface[j] = np.round(improved_surface[i] + float(improved_surface[opposite_match_index+1] - improved_surface[i])*(j-i)/(opposite_match_index+1-i))

            # now recompute jumps
            jumps = improved_surface[1:] - improved_surface[:-1]
            continue

        # IF THE ABOVE DIDN'T WORK, TRY THE 'JUMP' TECHNIQUE, SEEING WHETHER AN ANOMALOUS 'JUMP' IS COUNTERBALANCED BY AN
        # OPPOSITE AND (APPROXIMATELY) EQUAL JUMP IN THE OPPOSITE DIRECTION.
        # Don't worry about any trends less than 12 pixels.  Hills do that.
        jump = jumps[i]
        if abs(jump) < 5:
            continue

        # tune settings
        jump_lookahead = 50 # Number of pixels to look ahead to see if we can find a matching down-jump
        jump_magnitude_threshold = 0.50 # What fraction of the original jump the new jump has to be (in the opposite direction) to qualify.

        # see if we can find a jump in the near-future that crosses this threshold in the other direction.  If so, we've found our counter-part
        if jump < 0:
            opposite_jump_index = np.argmax((jumps[i:i+jump_lookahead]) > (-jump*jump_magnitude_threshold))
        elif jump > 0:
            opposite_jump_index = np.argmax((jumps[i:i+jump_lookahead]) < (-jump*jump_magnitude_threshold))

        if opposite_jump_index > 0:
            opposite_jump_index += i
        else: # If we didn't find a partner opposite offset, skip and move along.
            continue

        # Linearly interpolate, get to the closest pixel
        try:
            for j in range(i+1,opposite_jump_index+1):
                improved_surface[j] = np.round(improved_surface[i] + float(improved_surface[opposite_jump_index+1] - improved_surface[i])*(j-i)/(opposite_jump_index+1-i))
        except IndexError:
            print("i", i, "j", j, "opposite_jump_index", opposite_jump_index, improved_surface.shape, jumps.shape)
            # Break the program here.
            100/0

        # now recompute jumps
        jumps = improved_surface[1:] - improved_surface[:-1]
        continue
    return improved_surface

#---------------------------------------------------------------------------------------------
# original_indices, improved_indices = RetrackSurface_MacFerrin(data)
#
# Takes in:
# data - CReSIS radargram dictionary
#
# Returns:
# original_indices - 1D array of the vertical pixel location of the surface at each trace from 
# the original CReSIS surface tracking
# improved_indices - 1D array of the vertical pixel location of the surface at each trace from 
# the new surface retracking
#---------------------------------------------------------------------------------------------
def RetrackSurface_MacFerrin(data, starter):
    traces = 10*np.log10(data['Data'])
    improved_indices = np.zeros(data['Surface'].shape)
    original_indices = np.zeros(data['Surface'].shape)
    ind_count = 0
    for i in data['Surface']:
        original_indices[ind_count] = np.argmin(np.abs(data['Time'] - i))
        ind_count = ind_count + 1

    # 3) Perform surface pick crawling threshold behavior mask (assume a step-change analysis [goes from weak->strong at surface], and continuity of surface in 
    # original file.)
    # Create a step-change mask to optimze where the returns transition from "dark" to "bright"
    MASK_RADIUS = 50
    vertical_span_mask = np.empty([MASK_RADIUS*2,], dtype=np.float64)
    vertical_span_mask[:MASK_RADIUS] = -1.0
    vertical_span_mask[MASK_RADIUS:] = +3.0

    vertical_span_mask = vertical_span_mask * gaussian(np.arange(vertical_span_mask.shape[0]),mu=(MASK_RADIUS-5),sigma=(float(MASK_RADIUS)/3.0))

    # Expand the shape to handle array broadcasting below
    vertical_span_mask.shape = vertical_span_mask.shape[0], 1

    # This is the vertical window size of the extent of the search.  Should be bigger than any jump from one surface pixel to the next.
    MASK_SEARCH_RADIUS = 150

    improved_indices = np.empty(original_indices.shape, dtype=original_indices.dtype)

    # Start at the left suggested vertical pixel starting point
    if starter == 0:
        last_best_index = original_indices[0]
    else:
        last_best_index = starter
    
    # A template graph to use, just have to add in the center vertical index at each point and go from there.
    search_indices_template = np.sum(np.indices((vertical_span_mask.shape[0], 2*MASK_SEARCH_RADIUS)),axis=0) - MASK_SEARCH_RADIUS - MASK_RADIUS
    for i in range(traces.shape[1]):
        # Create an array of indices spanning the top-to-bottom of the MASK_SEARCH_RADIUS, and fanning out MASK_RADIUS above and below that point.
        search_indices = search_indices_template + last_best_index
        # Handle overflow indices if below zero or above max (shouldn't generally happen)... just assign to the top or bottom pixel
        search_indices[search_indices < 0] = 0
        search_indices[search_indices >= traces.shape[0]] = traces.shape[0]-1

        bestfit_sum = np.sum(traces[:,i][search_indices.astype(int)] * vertical_span_mask, axis=0)

        assert bestfit_sum.shape[0] == 2*MASK_SEARCH_RADIUS

        # Get the best fit (with the highest value from the transformation fit)
        last_best_index = search_indices[MASK_RADIUS,np.argmax(bestfit_sum)]
        improved_indices[i] = last_best_index

    # Erase most the little "jump" artifacts in the surface picker.
    improved_indices = get_rid_of_false_surface_jumps(improved_indices)
    improved_indices = first_maximum(improved_indices,traces)
    return original_indices.astype(int), improved_indices.astype(int)
