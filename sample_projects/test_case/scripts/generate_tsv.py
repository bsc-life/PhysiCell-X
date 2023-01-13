import numpy as np
import pandas as pd 

# generate 1000,50.000,500.000
# i will build up on the file o 1000 to keep consitency
# use uniform random to produce the positions 

no_points = [1000,50000,500000]
cell_radius = 0.0001
x_min = -1000
y_min = -1000
z_min = -1000
x_max = 1000
y_max = 1000
z_max = 1000

# points =[[]]
points = np.random.uniform(x_min+cell_radius,x_max-cell_radius,(3,1))
overlap = False
cells = 0
# for n in no_points:
# while cells<500000-1:
#     # print(cells)
#     # print(i)
#     # generate a point
#     new_point = np.random.uniform(x_min+cell_radius,x_max-cell_radius,(3,1))
#     # y dimension is number of alredy added points
#     # print(points.shape)
#     for i in range(0,points.shape[1]):
#         distance = np.sqrt(((points[0,i]-new_point[0,0])*(points[0,i]-new_point[0,0]))+
#                             ((points[1,i]-new_point[1,0])*(points[1,i]-new_point[1,0]))+
#                             ((points[2,i]-new_point[2,0])*(points[2,i]-new_point[2,0]))
#                             )

#         # break
#         if distance<=0.95 * 2.0 *cell_radius:            
#             break
#     if not overlap:
#         # print(new_point)
#         cells=cells+1
#     else:
#         print("overlap")    
#     points = np.append(points,new_point,axis=1)
# # pd.DataFrame(points).to_csv("/home/thalia/BSC/physicell_x/sample_projects/test_case/config/cells1000.csv")
# np.savetxt("~/physicell_x/sample_projects/test_case/config/cells500000.csv", np.transpose(points), delimiter="," ,fmt='%f')

# print(points.shape)


# implementation 2
points = np.random.uniform(x_min+cell_radius,x_max-cell_radius,(3,500000))
overlap = True
for i in range(0,points.shape[1]-1):
    distance = np.sqrt(((points[0,i]-points[0,i+1])*(points[0,i]-points[0,i+1]))+
                        ((points[1,i]-points[1,i+1])*(points[1,i]-points[1,i+1]))+
                        ((points[2,i]-points[2,i+1])*(points[2,i]-points[2,i+1]))
                        )
    if distance<=0.95 * 2.0 *cell_radius:
        print("overlap")
np.savetxt("cells500000_2.csv", np.transpose(points), delimiter="," ,fmt='%f')
        