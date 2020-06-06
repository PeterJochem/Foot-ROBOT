## README

### File structure

This repository contains hundreds of simulation results for plate intruding and extracting in granular materials.
I divide those data into three parts:

1. **Extract**:

This part includes data when the plate extract from deeper area of the granular materials onto the shallow area.

2. **Intrude**

This part includes data when the plate intrude from the top of granular materials into the deeper area.

3. **Small angles**

This part includes data when the plaate intrude or extract with a small gamma angles (-18.36,18.36,-33.7,33.7). Since simulation with small gamma angles require wider x dimension, I put small angles cases out of the extract and the intrude.

### How to understand the output csv files

The output csv files have two types **output_plate_forces.csv** and **output_plate_positions.csv** which respectively write down the position and resistive forces information of the plate in simulations. 

For the ouput forces csv files, before the simualtion starts, I firstly write down the gamma and beta information and then write down the time and three directions of resistive forces. So the first column records the time and the second, third and fourth columns record the x direction, y direction and z direction forces respectively. For example, we can read from pictures below: 

> Let's look at the first two rows. The gamma angle is -90 degrees and the beta is 0. The forces of three directions are zero N at 0s.

![IMG](https://github.com/HappyLamb123/Foot-ROBOT/blob/master/img/forces.PNG?raw=true)

Similarly with the output forces csv files, the output position files firtly write down the gamma and beta information and then and then write down the time and three directions of the positions of the center of mass of the plate. The only difference is that I also write down the z positions of the top sphere of the granular materials. So the first column records the time and the second, third and fourth columns record the x direction, y direction and z direction positions of the center of mass of the plate respectively. And also, the fifth colunm records the z positions of the top granular sphere. For example, we can read from pictures below: 

> Let's look at the first two rows. The gamma is -90 degrees and the beta is 0. The x direction position, y direction position and z direction position of the center of mass of the plate are 0, 0 and 12cm respectively. The z position of the top granular sphere is 8.73585cm.

![IMG](https://github.com/HappyLamb123/Foot-ROBOT/blob/master/img/positions.PNG?raw=true)












