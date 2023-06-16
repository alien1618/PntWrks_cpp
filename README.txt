========================================================================
To run in Debian:
========================================================================
1. Ensure the  required packages are installed:
sudo apt-get install bash gcc g++ make cmake libomp-dev openmpi-bin libgomp1 build-essential python3 python3-numpy python3-matplotlib python3-scipy

2. Edit permissions to shell script:
chmod u+x run.sh
chmod u+x cln.sh

3. Inside the "examples" folder, copy a case study and paste inside the PntWrks folder and rename the file to main.cpp

3. Inside the PntWrks folder directory run compiler script:
./run.sh

4. When the simulation is running, a folder called "out" is generated which has all the output data from the program. You can visualize the results in VTK format in Paraview, VisIt, or Salome-platform

4. If you have data output in TXT format, you can use the included python script "pltpnts.py" or "pltmsh.py" to export JPG images of the solution as such:
python3 pltmsh.py OR python3 pltpnts.py

NOTE: Ensure the name of the variable you want to plot (defined in line 11 of pltmsh.py or  pltpints.py)
matches the name of the results in the outputs folder. You also may need to edit the scale value to scale the plot as needed

5. To export the created JPG images inside the "pics" folder to video go inside the "pics" folder and in the terminal type:
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
