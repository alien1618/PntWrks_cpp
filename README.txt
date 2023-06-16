========================================================================
To run in Debian:
========================================================================
1. Ensure the  required packages are installed:
sudo apt-get install bash gcc g++ make cmake libomp-dev openmpi-bin libgomp1 build-essential python3 python3-numpy python3-matplotlib python3-scipy

2. Edit permissions to shell script:
chmod u+x run.sh

3. Inside the FPM2D folder directory run compiler script:
./run.sh

4. Exporting jpg images of solution is done by:
python3 plotmesh.py OR python3 plotpoints.py

NOTE: Ensure the name of the variable you want to plot (defined in line 11 of plotmesh.py or  plotpoints.py)
matches the name of the results in the outputs folder. You also may need to edit the scale value to scale the plot as needed

5. To export to video:
ffmpeg -i phi%d.jpg -vcodec mpeg4 phi.avi

where phi is the name of the image sequence you want to export


========================================================================
To run in FreeBSD:
========================================================================
1. Ensure the required packages are installed:
sudo pkg install gcc openmpi openmpi3 gmake cmake python python3 py38-numpy py38-scipy py38-matplotlib

2. Edit permissions to shell script:
chmod u+x run.sh

3. Run the compiler:
bash run.sh

4. Exporting jpg images of solution is done by:
python3 plotmesh.py OR python3 plotpoints.py

NOTE: Ensure the name of the variable you want to plot (defined in line 11 of plotmesh.py or  plotpoints.py)
matches the name of the results in the outputs folder. You also may need to edit the scale value to scale the plot as needed

5. To export to video:
ffmpeg -i phi%d.jpg -vcodec mpeg4 phi.avi

where phi is the name of the image sequence you want to export
