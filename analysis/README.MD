This folder contains notebooks for analysis and relevant python libraries.

Three libraries exist uner this directory:
csv_decoder.py reads .csv file produced by Josh's showering program, performs clustering and generates images
save_and_load.py saves and loads a fixed list of essential data needed for analysis into .npy files to avoid generating everything fresh every time. These essential data are listed as follow:

	new_background_event_list, new_signal_event_list, (arrays, saved and loaded with save() and load() from .npy)
	new_background_mass_list, new_signal_mass_list, (arrays, saved and from .npy)
    new_background_weight, new_signal_weight, (floats, saved and loaded from .npy)
    new_background_image_list, new_signal_image_list, (arrays, saved and loaded from .npy)
    new_background_recluster_images, new_signal_recluster_images (arrays, saved and loaded from .npy)

	new_background_event_list_clustered, new_signal_event_list_clustered, (list with c++ objects, need to be generated every time. less commonly used.)
    new_background_reclustered, new_signal_reclustered (list with c++ objects, need to be generated every time. less commonly used.)
        
substructure.py contains some methods for substructure variables.

Use of methods are documented in each .py files.