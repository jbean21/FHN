import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sn
from sklearn import svm
from mlxtend.plotting import plot_decision_regions

#Load in Global Order Parameter csv.
eq = pd.read_csv("OrderParameter.csv")

#Drop and Nan's if they exist
eq.dropna(how = 'all', axis = 1, inplace = True)

f= plt.subplots(figsize=(21,21))
sn.heatmap(eq.corr(),annot=True,fmt='.1f',color='green')

#Create the parameter vector
X = eq[['N','neigh','ndriv','k','thresh']]
#Create vector for order parameter
Y = eq[['p']]
#Create Support Vector Regression object
clf = svm.SVR()
#Fit the data, and then predict Y from given X, rbf is the kernel type
Y_rbf = clf.fit(X,Y).predict(X)
#Y_rbf is now a prediction the model outputs of the global order parameter
#We will want to analyse how good this data is.
#Take the difference between Y_rbf and actual data, and take the statistical moments
differences = np.abs(np.subtract(Y,Y_rbf))
mean_diff = np.mean(differences)
