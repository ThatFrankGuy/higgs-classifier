import numpy as np
import csv_decoder
import os

def save(folder_name, background_event_list, signal_event_list, background_mass_list, signal_mass_list,\
        background_weight, signal_weight, background_image_list, signal_image_list,\
        background_recluster_images, signal_recluster_images):
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
        
    # Saving all arrays with np.save (faster)
    np.save(folder_name+'/background_event_list.npy', background_event_list)
    np.save(folder_name+'/signal_event_list.npy', signal_event_list)
    np.save(folder_name+'/background_mass_list.npy', background_mass_list)
    np.save(folder_name+'/signal_mass_list.npy', signal_mass_list)
    np.save(folder_name+'/background_weight.npy', background_weight)
    np.save(folder_name+'/signal_weight.npy', signal_weight)
    np.save(folder_name+'/background_image_list.npy', background_image_list)
    np.save(folder_name+'/signal_image_list.npy', signal_image_list)
    np.save(folder_name+'/background_recluster_images.npy', background_recluster_images)
    np.save(folder_name+'/signal_recluster_images.npy', signal_recluster_images)  
    
def load(folder_name):
        
    # Saving all arrays with np.save (faster)
    new_background_event_list = np.load(folder_name+'/background_event_list.npy', allow_pickle=True)
    new_signal_event_list = np.load(folder_name+'/signal_event_list.npy', allow_pickle=True)
    new_background_mass_list = np.load(folder_name+'/background_mass_list.npy', allow_pickle=True)
    new_signal_mass_list = np.load(folder_name+'/signal_mass_list.npy', allow_pickle=True)
    new_background_weight = float(np.load(folder_name+'/background_weight.npy', allow_pickle=True))
    new_signal_weight = float(np.load(folder_name+'/signal_weight.npy', allow_pickle=True))
    new_background_image_list = np.load(folder_name+'/background_image_list.npy', allow_pickle=True)
    new_signal_image_list = np.load(folder_name+'/signal_image_list.npy', allow_pickle=True)
    new_background_recluster_images = np.load(folder_name+'/background_recluster_images.npy', allow_pickle=True)
    new_signal_recluster_images = np.load(folder_name+'/signal_recluster_images.npy', allow_pickle=True)
    
    return new_background_event_list, new_signal_event_list, new_background_mass_list, new_signal_mass_list,\
        new_background_weight, new_signal_weight, new_background_image_list, new_signal_image_list,\
        new_background_recluster_images, new_signal_recluster_images
        
# load events and recluster
def load_cluster(folder_name):
        
    # Saving all arrays with np.save (faster)
    new_background_event_list = np.load(folder_name+'/background_event_list.npy')
    new_signal_event_list = np.load(folder_name+'/signal_event_list.npy')
    
    new_background_event_list_clustered = cluster_event(background_event_list)
    new_signal_event_list_clustered = cluster_event(signal_event_list)
    
    new_background_reclustered = recluster_event(background_event_list_clustered)
    new_signal_reclustered = recluster_event(signal_event_list_clustered)

    
    return new_background_event_list_clustered, new_signal_event_list_clustered,\
        new_background_reclustered, new_signal_reclustered