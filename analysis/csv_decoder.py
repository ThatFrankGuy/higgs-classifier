# This library reads in the .csv files produced by Josh's pythia code.
# Note that all datasets are sorted BACKGROUND FIRST.

# Importing
import numpy as np
import math
import random
import os
import re
import time
import sys

from itertools import groupby
from io import StringIO
from mpl_toolkits.axes_grid1 import make_axes_locatable

import pyjet
from pyjet import cluster, DTYPE_PTEPM
from pyjet.testdata import get_event

# Reading CSV

# This function reads off the original CSV format: 
# pt, eta, phi, m, id, isCharged \n
# pt, eta, phi, m, id, isCharged \n
# ... (of all particles in one event)
# trimmed mass of identified jet \n
# \n
# and produces a event_list
# [[pt, eta, phi, m, id, isCharged], ...]
def return_event_list(fileName,max_read = float("inf"),weighted=0,pt_cut=1):
    
    printed = 0
    
    event_list = [];mass_list = [];weight_list = [];
    tmp_events = open(fileName).read().split("\n\n")[:-1]
    #print(len(tmp_events))
    for x in tmp_events:
        try:
            if len(event_list) == max_read: #FRANK# limit event size. If violated, interrupt the entire method
                return(event_list,mass_list)
            if weighted == 0:
                                
                #FRANK# rfind finds the last occurance
                #FRANK# [a:b] extract everything b/w a and bth elem, including ath and excluding bth
                #FRANK# why isn't the cut applied to weight and mass list? #ISSUE#
                
                mass_list.append(float(x[x.rfind("\n")+1:-1]))
                to_cut = np.array(np.genfromtxt(x[:x.rfind("\n")].splitlines(), delimiter=","))
                event_list.append(to_cut[[x[0] > pt_cut for x in to_cut]])
            else:
                weight_list.append(float(x[x.rfind("\n")+1:]))
                mass_list.append(float(x[x.rfind("\n",0,x.rfind("\n")-1)+1:x.rfind("\n")+1]))
                to_cut = np.genfromtxt(x[:x.rfind("\n")].splitlines(),delimiter=",")
                event_list.append(to_cut[[x[0] > pt_cut for x in to_cut]])
        except:
            print('Failed to decode CSV into np array. Make sure you have the correct format. ')
            print( sys.exc_info()[0])
            return
    if weighted == 0:
        #print(mass_list[0])
        return(event_list,mass_list)
    else:
        return(event_list,mass_list,weight_list)

#FRANK# Converting event_list to event image
def return_image_list(event_list, width=40, height=40):
    image_list = []
    image_0 = np.zeros((width,height)) #Charged pt #FRANK# I think it's labeled wrong here. Would it matter?
    image_1 = np.zeros((width,height)) #Neutral pt
    image_2 = np.zeros((width,height)) #Charged multiplicity
    
    #FRANK# pt, eta, phi, m, id, isCharged \n
    #FRANK# 0   1    2    3  4   5
    for z in range(len(event_list)):
        image_0 = np.zeros((width,height));image_1 = np.zeros((width,height));image_2 = np.zeros((width,height))
        for x in range(len(event_list[z])):
            phi_index = math.floor(width*event_list[z][x,2]//(2*math.pi)+width//2)
            eta_index = math.floor(height*event_list[z][x,1]//10+height/2) #FRANK# // is integer divide
            eta_index = min(eta_index,height-1)
            eta_index = max(0,eta_index)
            phi_index = int(phi_index);eta_index = int(eta_index)
            if (event_list[z][x,5] == 0):  #FRANK# neutral
                image_0[phi_index,eta_index] = image_0[phi_index,eta_index] + event_list[z][x,0]
            elif (event_list[z][x,5] == 1):  #FRANK# charged
                image_1[phi_index,eta_index] = image_1[phi_index,eta_index] + event_list[z][x,0]
                image_2[phi_index,eta_index] = image_2[phi_index,eta_index] + 1
        image_0 = np.divide(image_0,np.sum(image_0))
        image_1 = np.divide(image_1,np.sum(image_1))
        image_2 = np.divide(image_2,np.sum(image_2))
        image_list.append(np.array([image_0,image_1,image_2]))
    return(image_list)

# This method loads event files to produce:
# event_list, mass_list, image_list, files_Read and weight_list
def load_events(path, contains, debug = 0, max_read = float("inf"), max_files = float("inf"), weighted=0, \
                pt_cut = 0, width=40, height=40):
    print('Loading .csv event files from ' + path + ' containing \"' + contains + '\"')
    reading_event_list,reading_mass_list = [],[]
    reading_image_list = []
    reading_weight_list = []
    files_Read = 0
    print('List of files is: '+ str(os.listdir(path)))
    for i in os.listdir(path):
        print('Currently reading: '+ str(i))
        
        if (files_Read == max_files):
            break
        if len(reading_event_list) >= max_read:
            return(reading_event_list,reading_mass_list,reading_image_list,files_Read)
        if 'swp' in i:
            continue
        if os.path.isfile(os.path.join(path,i)) and (contains) in i:
            if debug==1:
                print(i)
                print(os.path.join(path,i))
            if weighted==0:
                print(i)
                temp_event_list,temp_mass_list = return_event_list(os.path.join(path,i),pt_cut=pt_cut)
            else:
                temp_event_list,temp_mass_list,temp_weight_list = return_event_list(os.path.join(path,i),
                                                                                    weighted=1,pt_cut=pt_cut)
                reading_weight_list = reading_weight_list = temp_weight_list
            temp_image_list = return_image_list(temp_event_list, width, height)
            if (len(temp_image_list) != len(temp_mass_list)):
                print("Image production failure: #image != #weight")
                print("File: " + str(os.path.join(path,i)))
            reading_event_list = reading_event_list + temp_event_list
            reading_mass_list = reading_mass_list + temp_mass_list
            reading_image_list = reading_image_list + temp_image_list
            files_Read = files_Read + 1
            print(str(files_Read) + "files processed.")
    if weighted==0:
        return(reading_event_list,reading_mass_list,reading_image_list,files_Read)
    else:
        return(reading_event_list,reading_mass_list,reading_image_list,files_Read,reading_weight_list)

# This is not used. 
def return_fine_image_list(event_list, event_list_clustered, granularity, which_jet = 0, width=40, height=40):
    image_list = []
    image_0 = np.zeros((width,height)) #Charged pt
    image_1 = np.zeros((width,height)) #Neutral pt
    image_2 = np.zeros((width,height)) #Charged multiplicity

    for z in range(len(event_list)):
        image_0 = np.zeros((width,height))
        image_1 = np.zeros((width,height))
        image_2 = np.zeros((width,height))
        for x in range(len(event_list[z])):
            
            try:
                phi_index = (event_list[z][x,2]-event_list_clustered[z][which_jet].phi)
            except:
                print(z)
            #At this point, phi_index is just delta_phi, which could be anywhere from -2pi to 2pi
            if (phi_index % (2*math.pi) >= (width//2)*granularity) and (phi_index % (2*math.pi) <= 2*math.pi-(width//2)*granularity):
                continue
                #This gets rid of the delta phi's that are far away from the jet
            phi_index = phi_index % (2*math.pi)
            if phi_index > math.pi:
                 phi_index = phi_index - 2*math.pi   
            phi_index = int(math.floor(phi_index/granularity)) #should be good now
            if (phi_index > (width//2)) or (phi_index < -(width//2)):
                print(phi_index)
            phi_index = phi_index + (width//2)
            

            eta_index = int(math.floor((event_list[z][x,1]-event_list_clustered[z][which_jet].eta)/granularity) + height//2)
            if eta_index >= height:
                continue
            if eta_index < 0:
                continue
            
            #finally, lets fill
            if (event_list[z][x,5] == 0):
                image_0[phi_index,eta_index] = image_0[phi_index,eta_index] + event_list[z][x,0]
            elif (event_list[z][x,5] == 1):
                image_1[phi_index,eta_index] = image_1[phi_index,eta_index] + event_list[z][x,0]
                image_2[phi_index,eta_index] = image_2[phi_index,eta_index] + 1

        #Now, lets go through and normalise to 255
        image_0 = np.divide(image_0,np.sum(image_0))
        image_1 = np.divide(image_1,np.sum(image_1))
        image_2 = np.divide(image_2,np.sum(image_2))
        image_list.append(np.array([image_0,image_1,image_2]))
    return(image_list)

def cluster_event(event_list):
    event_list_clustered = []
    for x in range(len(event_list)):
        to_Cluster = np.array([event_list[x][:,0],event_list[x][:,1],event_list[x][:,2],event_list[x][:,3]])
        to_Cluster = np.swapaxes(to_Cluster,0,1)
        to_Cluster = np.core.records.fromarrays(to_Cluster.transpose(), 
                                             names='pT, eta, phi, mass',
                                             formats = 'f8, f8, f8,f8')
        sequence_Cluster = cluster(to_Cluster, R = 0.8,p = -1)
        jets_Cluster = sequence_Cluster.inclusive_jets()
        event_list_clustered.append(jets_Cluster)
    return(event_list_clustered)

def recluster_event(cluster_list):
    reclustered_list= []
    for i in range(len(cluster_list)):
        sequence_Cluster = cluster((cluster_list[i][0]), R=0.2,p=-1)
        jets_Cluster = sequence_Cluster.inclusive_jets()
        reclustered_list.append(jets_Cluster)
    return(reclustered_list)

def return_fine_image_list_reclustered(event_list, event_list_clustered, radius, which_jet = 0,verbose = False, width=40, height=40):
    image_list = []
    image_0 = np.zeros((width,height)) #Neutral pt
    image_1 = np.zeros((width,height)) #Charged pt
    image_2 = np.zeros((width,height)) #Charged multiplicity
    
    no_two = 0

    for z in range(len(event_list)):
        image_0 = np.zeros((width,height))
        image_1 = np.zeros((width,height))
        image_2 = np.zeros((width,height))
        
        if (len(event_list_clustered[z]) > 1):
            #First, let's find the direction of the second-hardest jet relative to the first-hardest subjet
            phi_dir = -(dphi(event_list_clustered[z][1].phi,event_list_clustered[z][0].phi))
            eta_dir = -(event_list_clustered[z][1].eta - event_list_clustered[z][0].eta)
            #Norm difference:
            norm_dir = np.linalg.norm([phi_dir,eta_dir])
            #This is now the y-hat direction. so we can actually find the unit vector:
            y_hat = np.divide([phi_dir,eta_dir],np.linalg.norm([phi_dir,eta_dir]))
            #and we can find the x_hat direction as well
            x_hat = np.array([y_hat[1],-y_hat[0]]) 
        else:
            no_two = no_two + 1
            #continue
            
        if verbose==True:
            print(x_hat,y_hat,norm_dir)
            
        
        for x in range(len(event_list[z])):
            if (len(event_list_clustered[z]) == 1):
                #In the case that the reclustering only found one hard jet (that seems kind of bad, but hey)
                #no_two = no_two+1
                new_coord = [dphi(event_list[z][x,2],event_list_clustered[z][0].phi),event_list[z][x,1]-event_list_clustered[z][0].eta]
                indxs = [math.floor(width*new_coord[0]/(radius*1.5))+width//2,math.floor(height*(new_coord[1])/(radius*1.5))+height//2]
            else:
                #Now, we want to express an incoming particle in this new basis:
                part_coord = [dphi(event_list[z][x,2],event_list_clustered[z][0].phi),event_list[z][x,1]-event_list_clustered[z][0].eta]
                new_coord = np.dot(np.array([x_hat,y_hat]),part_coord)
                #Now, we want to cast these new coordinates into our array
                indxs = [math.floor(width*new_coord[0]/(radius*1.5))+width//2,math.floor(height*(new_coord[1]+norm_dir/1.5)/(radius*1.5))+height//2]
                
            if indxs[0] >= width or indxs[1] >= height or indxs[0] <= 0 or indxs[1] <= 0:
                continue
            phi_index = int(indxs[0]); eta_index = int(indxs[1])
            #finally, lets fill
            if (event_list[z][x,5] == 0):
                image_0[phi_index,eta_index] = image_0[phi_index,eta_index] + event_list[z][x,0]
            elif (event_list[z][x,5] == 1):
                image_1[phi_index,eta_index] = image_1[phi_index,eta_index] + event_list[z][x,0]
                image_2[phi_index,eta_index] = image_2[phi_index,eta_index] + 1

        #Now, lets go through and normalise to 255
        if (np.sum(image_0) == 0 or np.sum(image_1) == 0 or np.sum(image_2) == 0):
            image_list.append(np.array([image_0,image_1,image_2]))
            continue
        image_0 = np.divide(image_0,np.sum(image_0))
        image_1 = np.divide(image_1,np.sum(image_1))
        image_2 = np.divide(image_2,np.sum(image_2))
        image_list.append(np.array([image_0,image_1,image_2]))
    print("Number of events with only one constituent in leading jet: " + str(no_two))
    return(image_list)

# Correct phi range
def fix_phi(phi):
    while phi > math.pi:
        phi = phi - 2*math.pi
    while phi < -math.pi:
        phi = phi + 2*math.pi
    return phi

# Returns the difference in phi between phi, and phi_center
# as a float between (-PI, PI)
def dphi(phi,phi_c):
    
    dphi_temp = phi - phi_c
    while dphi_temp > np.pi:
        dphi_temp = dphi_temp - 2*np.pi
    while dphi_temp < -np.pi:
        dphi_temp = dphi_temp + 2*np.pi
    return (dphi_temp)

# Rapidity
def y(p):
    return ((1/2)*math.log((p.e+p.pz)/(p.e-p.pz)))

def R(con1,con2):
    return (((con1.eta-con2.eta)**2+dphi(con1.phi,con2.phi)**2)**(1/2))

def R_y(con1,con2):
    return (((y(con1)-y(con2))**2+dphi(con1.phi,con2.phi)**2)**(1/2))

def N_2(jcon):
    #Takes jcon, jet constituents
    p_x_total = np.sum([con.px for con in jcon])
    p_y_total = np.sum([con.py for con in jcon])
    p_total = (p_x_total**2+p_y_total**2)**(1/2)
    
    v_1e2 = 0
    for i in range(len(jcon)):
        for j in range(i+1,len(jcon)):
            v_1e2 = v_1e2+ jcon[i].pt*jcon[j].pt*R(jcon[i],jcon[j])/(p_total**2)
    v_2e3 = 0
    for i in range(len(jcon)):
        for j in range(i+1,len(jcon)):
            for k in range(j+1,len(jcon)):
                v_2e3 = v_2e3 + jcon[i].pt*jcon[j].pt*jcon[j].pt*min(R(jcon[i],jcon[j])*R(jcon[i],jcon[k]),
                                                                     R(jcon[j],jcon[k])*R(jcon[i],jcon[j]),
                                                        R(jcon[i],jcon[k])*R(jcon[j],jcon[k]))/(p_total**3)
    return v_2e3/(v_1e2**2)

# Softdrop

class myJet(object):
    def __init__(self,px,py,pz):
        self.px = px; self.py = py; self.pz = pz; self.pt = (px**2+py**2)**(1/2)
        self.phi = math.atan2(py,px); self.eta = -math.log(math.tan(math.atan2(self.pt,self.pz)/2));

class Node(object):
    def __init__(self,data):
        self.data = data
        self.children = []
    def add_child(self,obj):
        self.children.append(obj)

def softdrop(jcon,z=0.1,debug = 0):
    #Takes the constituents of a jet, and softdrops it.
    #First, we need to step through the jet and build the tree of clustering
    #Since we are reclustering the whole thing; just take R = 1; i.e. dont need to think about it
    def distance(con1,con2):
        return R(con1,con2)**2
    pseudojets = []
    nodes = []
    for con in jcon:
        x = Node([con,1])
        #pseudojets.append(x)
        nodes.append(x) #1 means its still a pseudojet; i.e. not been clustered
    def how_many_pseudo(nodes):
        how_many = 0
        for node in nodes:
            if node.data[1] == 1:
                how_many = how_many + 1
        return how_many
    if debug == 1:
        print("len(nodes) : " + str(len(nodes)))
        print(nodes)
    rep = 0
    while how_many_pseudo(nodes) > 1:
        #print(how_many_pseudo(nodes))
        min_distance = float("inf")
        min_index = [0,0]
        for i in range(0,len(nodes)):
            if nodes[i].data[1] == 0: #Its already part of something else
                continue
            for j in range(i+1,len(nodes)):
                if nodes[j].data[1] == 0:
                    continue
                if distance(nodes[i].data[0],nodes[j].data[0]) < min_distance:
                    min_index[0] = i; min_index[1] = j; 
                    min_distance = distance(nodes[i].data[0],nodes[j].data[0])
        i = min_index[0];j=min_index[1];
        new_node = Node([myJet(nodes[i].data[0].px+nodes[j].data[0].px,
                           nodes[i].data[0].py+nodes[j].data[0].py,
                           nodes[i].data[0].pz+nodes[j].data[0].pz),1])
        new_node.add_child(nodes[i])
        new_node.add_child(nodes[j])
        nodes.append(new_node)
        nodes[i].data[1] = 0
        nodes[j].data[1] = 0
    
    #print(nodes)

    to_check = [] #nodes to check
    for i in range(len(nodes)):
        if nodes[i].data[1] == 1:
            to_check.append(nodes[i])
    softcon = []
    while len(to_check) > 0:
        #print(to_check)
        our_childs = to_check[0].children
        #print(our_childs)
        if len(our_childs) == 0:
            softcon.append(to_check[0].data[0])
            to_check.pop(0)
            continue
        if min(our_childs[0].data[0].pt,our_childs[1].data[0].pt)/  \
                        (our_childs[0].data[0].pt+our_childs[1].data[0].pt) > z:
            to_check.append(our_childs[0])
            to_check.append(our_childs[1])
        elif our_childs[0].data[0].pt > our_childs[1].data[0].pt:
            to_check.append(our_childs[0])
        else:
            to_check.append(our_childs[1])
        to_check.pop(0)
    return softcon


# Take in a tuple of image lists, normalze and zero zenter all of them.
def zero_center_and_normalize(image_lists):
    tmp_av = np.average(np.concatenate(image_lists), axis=0)
    tmp_sd = np.std(np.concatenate(image_lists), axis=0)
    for j in range(len(image_lists)):
	    for i in range(len(image_lists[j])):
	        image_lists[j][i] = np.divide((image_lists[j][i] - tmp_av), (tmp_sd+1e-5)) #perhaps add some r to temp_sd to suppress noise
    return image_lists

    # Normalize and zero-center a pair (background and sample) of images
def zero_center_and_normalize_pair(background_images, signal_images):
    tmp_av = np.average(np.concatenate((background_images, signal_images)), axis=0)
    tmp_sd = np.std(np.concatenate((background_images, signal_images)), axis=0)
    for i in range(len(background_images)):
        background_images[i] = np.divide((background_images[i] - tmp_av), (tmp_sd+1e-5)) #perhaps add some r to temp_sd to suppress noise
    for i in range(len(signal_images)):
        signal_images[i] = np.divide((signal_images[i] - tmp_av), (tmp_sd+1e-5))#/tmp_sd
    return background_images, signal_images

