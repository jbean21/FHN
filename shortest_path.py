import glob
import numpy as np
import FileLoader2

print("Starting Analysis")

progress_counter = 0

output = np.zeros([49,2])

for i in range(49): output[i,0] = i

for ii in sorted(glob.glob("/media/jay/0FD90FF80FD90FF83/PROJECTDATA/2DPERIODIC_2/Data_all/*.csv")):
  df, params = FileLoader2.FileLoader2(ii)
  D = int(params[2])
  N = params[0]
  
  #Track progress of analysis
  #print(progress_counter)
  #Check if first letter of filename is F for FHN, index 38 to remove path from consideration.
  #print ii
  if ii[63] == "d" and D>1 and progress_counter%100==0 and N<101: #N==4  and
    print ii
    v = np.asarray(df[0])
    sqrtN = int(np.sqrt(N))
    
    for i in range(v.shape[0]):
      for k in range(sqrtN):
        if v[i]>=k*sqrtN and v[i]<k*2*sqrtN:
          row_i = k
          break
  
      column_i = v[i]%sqrtN
  
      dists = []
      
      for j in range(v.shape[0]):
        for k in range(sqrtN):
          if v[j]>=k*sqrtN and v[j]<k*2*sqrtN:
            row_j = k
            break
        
        column_j = v[j]%sqrtN
        
        row_j_2 = 9999
        column_j_2 = 9999
        
        #print row_j, column_j
        if row_j>=(sqrtN)/2.0:    row_j_2    -= sqrtN
        if column_j>=(sqrtN)/2.0: column_j_2 -= sqrtN
        
        if (row_j_2-row_i)<(row_j-row_i):             row_j=row_j_2
        if (column_j_2-column_i)<(column_j-column_i): column_j=column_j_2
        
        #print row_j, column_j
        dists.append(np.sqrt((row_j-row_i)**2+(column_j-column_i)**2))
    
    dists = filter(lambda a: a != 0, dists)
    print dists
    shortest_dist = min(dists)
    print "Distance: ", shortest_dist
    output[D-1,1] += shortest_dist
    print "\n----"
  
  progress_counter += 1

np.savetxt("/media/jay/0FD90FF80FD90FF83/PROJECTDATA/2DPERIODIC_2/SHORTEST_PATH_"+str(N)+"_DIAG.csv", output, delimiter=',', newline='\n', fmt="%0.3f")







