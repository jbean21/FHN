"""
Python script to load data from csv's and then calculate the order parameter as defined in the literature
"""
import numpy as np
import pandas as pd
from scipy.signal import hilbert
import glob
import csv

def expip(data):
	#Subtract the mean from data to ensure that Argand diagram centred on 0
	mean = np.mean(data)
	data = np.subtract(data,mean)
    #calculates the analytic signal and then exp(i*phi)
	analytic_signal = hilbert(data)
	i_p = np.unwrap(np.angle(analytic_signal))
	exp_ip = list(map(lambda x: np.exp(1j*x),i_p))
	return [exp_ip, analytic_signal]

def aggregate(dataframe):
    #This function will load in a dataframe, calculate expip on each time
    #series, and then sum them.
    df = pd.read_csv(dataframe, header = 0)
    N = len(df.columns) 
    data = []
    for i in range(1,N): #start from 1 to negate time column
        data.append(expip(df.iloc[:,i])[0])
    orderParam = [np.sum(x)/N for x in zip(*data)]
    OPavg = np.mean(orderParam)
    orderParam.append(OPavg)
    return orderParam

def globbing():
	#loads in all csv files in a folder and calculates the order parameter
	#for each, and saves the resulting time series as a csv.
    filenames = glob.glob("*.csv")
    if len(filenames) == 0:
        print("No CSV files present!")
    for i in filenames:
        orderParam = aggregate(i)
        np.savetxt("OP"+i,np.asarray(orderParam),delimiter = ",\n")
        #now want to extract the system data from the csv and record the parameters along with average order parameter
        x = i.split("_")
        N = x[1]
        neigh = x[2]
        locy = x[3]
        k = x[4]
        thresh = x[5]
        #calculate average order param
        OPAVG = np.mean(orderParam)
        fields = [N,neigh,locy,k,thresh,OPAVG]
        with open("OrderParameter.csv", 'a') as f:
            writer = csv.writer(f)
            writer.writerow(fields)

def localOP(filename):
	#calculates the order parameter locally
    #read in the csv of data
	df = pd.read_csv(filename, header = 0)
    #read in number of neurons
	N = len(df.columns)
	data = []
	for i in range(1,N): #start from 1 to negate time column
		temp = []
		temp.append(expip(df.iloc[:,i])[0])
		temp.append(expip(df.iloc[:,i-1])[0])
		temp.append(expip(df.iloc[:,i+1])[0])
		data.append(np.mean([sum(x)/3 for x in zip(*temp)])) #average local order parameter for each neuron, calculated only from the nearest neighbours
	x = dataframe.split("_")
	N = x[1]
	neigh = x[2]
	locy = x[3]
	k = x[4]
	thresh = x[5]
	fields = [N,neigh,locy,k,thresh]
	with open("LocalOrderParameter.csv", 'a') as f:
		writer = csv.writer(f)
		writer.writerow(fields)
		writer.writerow(data)






