# parses the output slurm file and plots the time between simulation timepoints
import re
import matplotlib.pyplot as plt
import pandas as pd
output_files = ["/gpfs/projects/bsc08/shared_projects/PerMedCoE/debug_physicell_X/results_1000/output_vanila_1_1_48/output-27086249",
"/gpfs/projects/bsc08/shared_projects/PerMedCoE/debug_physicell_X/results_1000/output_n1_ntasks1_t48_1000/output-27073163",
"/gpfs/projects/bsc08/shared_projects/PerMedCoE/debug_physicell_X/results_1000/output_n1_ntasks2_t24_1000/output-27073164"]

dfs= []
for output_file in output_files:
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
    
    dfs.append(df)
# dfs[0]['total'] = dfs[0]['total']/60
# dfs[1]['total'] = dfs[1]['total']/60
# dfs[0]['total'] = dfs[0]['total']/60
# for i in range(len(output_files)):
#     dfs[i]['total'] = dfs[i]['total']/60

# ax =dfs[0].plot(x="total", y=["agents"])

# # ax =dfs[1].plot(x="timepoint", y=["agents"])
# dfs[1].plot(x="total", y=["agents"],ax=ax)
# dfs[2].plot(x="total", y=["agents"],ax=ax,title="Physicell Vanilla vs PhysiCell-X using 1 node")

# ax.set_ylabel("Number of Agents")
# ax.set_xlabel("Minutes (Computational)")
# plt.legend(['PhysiCell with 48 OpenMP Threads','PhysiCell-X 1 MPI 48 OpenMP Threads','PhysiCell-X 2 MPI 24 OpenMP Threads'])
# plt.savefig("agents-real-time.png")
# plt.show()

ax =dfs[0].plot(x="timepoint", y=["agents"])

dfs[1].plot(x="timepoint", y=["agents"],ax=ax)
dfs[2].plot(x="timepoint", y=["agents"],ax=ax,title="Physicell Vanilla vs PhysiCell-X using 1 node")
ax.set_ylabel("Number of Agents")
ax.set_xlabel("Minutes (Simulation)")
plt.legend(['PhysiCell with 48 OpenMP Threads','PhysiCell-X 1 MPI 48 OpenMP Threads','PhysiCell-X 2 MPI 24 OpenMP Threads'])
plt.savefig("agents-sim-time.png")
plt.show()