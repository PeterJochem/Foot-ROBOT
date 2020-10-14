# Demo code for foot intrude and extact

### 1. In both code, the foot (rigid plate) moves at constant speed with fixed body orientations.
### 2. To get the datasets for different intrusion and extract speed, you only need to change the  intrusion and extract velocities where I comment in those codes.
### 3. You should also choose the time_end carefully to make sure the plate is always move inside of the granular media. In addition, you can also change the initial positions of the plate to keep it inside the granular box

# Creating a Dataset with Bash Scripts 
We want to avoid having our program stop because of our ssh connection breaking. The simplest way to do this is to use tmux. If we create a tmux session on Kronos and then run our bash script within the tmux session, then even if our ssh connection breaks, the bash script/it's children processes will still run. <br />

To view if there are any tmux sessions active, run ```tmux ls```. If one is active, we can use it. Attatch to the tmux session by running ```tmux attach-session -t <session_number/id>```. If there are no active tmux session, create one by running ```tmux new -s <session_name>```. Then attatch to it. A great, short intro guide/cheatsheet to tmux can be found [here](https://linuxize.com/post/getting-started-with-tmux/) <br />


There is a seperate bash script for each of intrude and extract. To generate the full dataset, run these sequentially in tmux sessions. 
