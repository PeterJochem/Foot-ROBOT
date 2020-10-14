# How to Generate an Intrusion Dataset 
This folder has all the files neccesary to generate a dataset of intrusions. Every 20 timesteps of the simulation, we record Gamma, Beta, Depth, Position_x, Position_z, Velocity x, Velocity y, Velocity Z, GRF X, GRF Y, GRF Z into data_sets/dset.csv where Position is the foot's position, velocity is the foot's velocity, depth is the top of the plate (in z direction), and GRF is the ground reaction force exerted by the granular material on the plate. 

# Building
Juntao wrote a great [description](https://github.com/HappyLamb123/Foot-ROBOT) of how to configure the build for our purposes. It contains everything to get our simulations compiled and running. 

# intrude.cpp
This is the Chrono code that simulates the foot (or plate) interacting with the soft ground. It has a few command line arguments. The first argument is gamma and the second is the magnitude of the velocity with which the plate moves through the bed of granular material. Gamma is measured in radians. The speed argument should always be positive. It is in cm/s. To run it, do ```./myexe <gamma> <speed>``` For example, ```./myexe 3.14 3.0``` means that gamma is 3.14 rads and the plate is moving at a speed of 3cm/s. Each run of intrude.cpp has a constant gamma and intrusion speed and then varies the beta value over a range of values. As of now, each run uses each of beta = {PI/2.0, PI/3.0, PI/6.0, 0, -PI/6.0, -PI/3.0}. 

# setup.json
This is read by intrude.cpp and specifies some of the parameters of the simulation. It determines physical properties of the soft ground particles. 

# createDataSet.sh
This is a bash script for generating a dataset over a range of diffrent (gamma, beta, speed) triples. For both of gamma and speed, it takes 3 parameters. It takes a min value, max value, and numItems parameter. The min and max values specify the interval over which to vary the parameter. The numItems parameter specefies how many groups to split the above interval into and sample from. For example, if the speed min value is 0.1, it's max value is 1.0, and its numItems value is 5, then the bash script generates data with 5 seperate speeds evenly spaced in the interval of [0.1, 1.0] (cm/s). <br /> 
Unfortunately, we have combinatorial explosion. Each run of the intrude.cpp file takes a fixed gamma and speed value. Within a single run of the program, we vary the beta value. This allows us to generate data over the space of (gamma, beta, speed) triples. But, if the numItems parameter is 5 for speed and 10 for gamma, then we need to run 50 seperate instances of the program with new command line arguments. Adding another item to the input space multiplies the time to generate the data. For example, if we varied the dimensions of the plate, we would want to generate data evenly from the quadruple of (gamma, beta, speed, plate_shape). The number of trials would be the product of the discretization factor of each parameter. <br />
Remember to make this executable with ```chmod +x createDataSet.sh```. 


# Transferring Dataset from Kronos Your Local Machine
```scp -r kronos@10.105.107.29:<path to data folder>/data_sets </some/local/directory/on/your/cpu>```