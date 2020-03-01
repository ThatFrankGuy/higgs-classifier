import numpy as np
import csv_decoder
import os

def save_binary(name_1='background', name_2 = 'signal', folder_name, process_1_event_list, process_2_event_list, process_1_mass_list, process_2_mass_list,\
        process_1_weight, process_2_weight, process_1_image_list, process_2_image_list,\
        process_1_recluster_images, process_2_recluster_images):
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
        
    # Saving all arrays with np.save (faster)
    np.save(folder_name+'/'+name_1+'_event_list.npy', process_1_event_list)
    np.save(folder_name+'/'+name_2+'_event_list.npy', process_2_event_list)
    np.save(folder_name+'/'+name_1+'_mass_list.npy', process_1_mass_list)
    np.save(folder_name+'/'+name_2+'_mass_list.npy', process_2_mass_list)
    np.save(folder_name+'/'+name_1+'_weight.npy', process_1_weight)
    np.save(folder_name+'/'+name_2+'_weight.npy', process_2_weight)
    np.save(folder_name+'/'+name_1+'_image_list.npy', process_1_image_list)
    np.save(folder_name+'/'+name_2+'_image_list.npy', process_2_image_list)
    np.save(folder_name+'/'+name_1+'_recluster_images.npy', process_1_recluster_images)
    np.save(folder_name+'/'+name_2+'_recluster_images.npy', process_2_recluster_images)  
    
def load_binary(name_1='background', name_2 = 'signal', folder_name):
        
    # Saving all arrays with np.save (faster)
    new_process_1_event_list = np.load(folder_name+'/'+name_1+'_event_list.npy', allow_pickle=True)
    new_process_2_event_list = np.load(folder_name+'/'+name_2+'_event_list.npy', allow_pickle=True)
    new_process_1_mass_list = np.load(folder_name+'/'+name_1+'_mass_list.npy', allow_pickle=True)
    new_process_2_mass_list = np.load(folder_name+'/'+name_2+'_mass_list.npy', allow_pickle=True)
    new_process_1_weight = float(np.load(folder_name+'/'+name_1+'_weight.npy', allow_pickle=True))
    new_process_2_weight = float(np.load(folder_name+'/'+name_2+'_weight.npy', allow_pickle=True))
    new_process_1_image_list = np.load(folder_name+'/'+name_1+'_image_list.npy', allow_pickle=True)
    new_process_2_image_list = np.load(folder_name+'/'+name_2+'_image_list.npy', allow_pickle=True)
    new_process_1_recluster_images = np.load(folder_name+'/'+name_1+'_recluster_images.npy', allow_pickle=True)
    new_process_2_recluster_images = np.load(folder_name+'/'+name_2+'_recluster_images.npy', allow_pickle=True)
    
    return new_process_1_event_list, new_process_2_event_list, new_process_1_mass_list, new_process_2_mass_list,\
        new_process_1_weight, new_process_2_weight, new_process_1_image_list, new_process_2_image_list,\
        new_process_1_recluster_images, new_process_2_recluster_images
        
# load events and recluster
def load_cluster_binary(name_1='background', name_2 = 'signal', folder_name):
        
    # Saving all arrays with np.save (faster)
    new_process_1_event_list = np.load(folder_name+'/'+name_1+'_event_list.npy')
    new_process_2_event_list = np.load(folder_name+'/'+name_2+'_event_list.npy')
    
    new_process_1_event_list_clustered = cluster_event(process_1_event_list)
    new_process_2_event_list_clustered = cluster_event(process_2_event_list)
    
    new_process_1_reclustered = recluster_event(process_1_event_list_clustered)
    new_process_2_reclustered = recluster_event(process_2_event_list_clustered)

    
    return new_process_1_event_list_clustered, new_process_2_event_list_clustered,\
        new_process_1_reclustered, new_process_2_reclustered