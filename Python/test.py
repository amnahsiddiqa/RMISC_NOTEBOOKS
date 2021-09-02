
## Random snippets for my short term memory loss
## exporting nested list 

import csv

data = [[74.9818, 0.09],
 [75.0021, 0.18],
 [75.0022, 0.17],
 [75.0026, 0.07],
 [75.0027, 0.07],]


pd.DataFrame(data).to_csv("output.csv",index=False)


## printing nested list
print(statistics[:3])  
