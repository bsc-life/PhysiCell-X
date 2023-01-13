# parses the output slurm file and plots the time between simulation timepoints
import re
import matplotlib.pyplot as plt
import pandas as pd
output_file = "/home/thalia/BSC/benchmark_x_results/output-27005269/output-27005269"
file = open(output_file, 'r')
file = file.readlines()
list_dicts = []
obj = {}
for line in file:
    interval = {}
    total = {}
    if len(obj)==4:
        list_dicts.append(obj)
        obj = {}
    if "current simulated time" in line:
        timepoint= re.findall('\d+\.?\d*',line)[0]
        obj['timepoint'] = int(timepoint)
    if "total agents" in line:
        agents= re.findall('\d+\.?\d*',line)[0]
        obj['agents'] = int(agents)
    if "interval wall time" in line:
        res = re.findall('\d+\.?\d*',line)
        interval = 60*60*24*float(res[0]) + 60*60*float(res[1]) + 60*float(res[2]) + float(res[3])
        obj['interval'] = float(interval)
        print(interval)
        
    if "total wall time" in line:
        res = re.findall('\d+\.?\d*',line)
        total = 60*60*24*float(res[0]) + 60*60*float(res[1]) + 60*float(res[2]) + float(res[3])
        obj['total'] = float(total)
        
# turn interval  interval to seconds
df = pd.DataFrame(list_dicts)
df['sum_intervals'] = df['interval'].cumsum()
ax =df.plot(x="timepoint", y=["interval" ,"total"])
ax.set_ylabel("seconds")
print(df)
plt.savefig("/home/thalia/BSC/benchmark_x_results/output-27005269/27005269-times.png")
plt.show()
