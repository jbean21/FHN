import glob
import FileLoader
import OrderParameter
#import CorrelationPlots
#import ContourPlots
#import OrderVsDisorder

print("Starting Analysis")

FileLoader.OPinit()
progress_counter = 0
for i in sorted(glob.glob("/media/jay/0FD90FF80FD90FF83/PROJECTDATA/2DPERIODIC_2/Data/*.csv")):
  #Track progress of analysis
  print(progress_counter)
  #Check if first letter of filename is F for FHN, index 38 to remove path from consideration.
  if i[59] == "v":
    df, params = FileLoader.FileLoader(i)
    print(params)
    OrderParameter.OP(i,df,params[0],params[1],params[2],params[3],params[4])
 
  progress_counter += 1
op_csv = FileLoader.OPread()
#CorrelationPlots.Correlation(op_csv)
#ContourPlots.ContourPlots(op_csv)
#OrderVsDisorder.OrderVsDisorder(op_csv)







