Siemens MAMMOMAT Revelation inputs for the GATE 10 room-planning study

Keep gate_mammography and L2_container_v2 together inside Revelation_Scanner_Inputs.
They contain the native scanner modules, 28 and 49 kVp spectra, material tables,
L2 surface meshes, sphere positions and insert positions used in the room study.

The student launcher reads the file exported from the room viewer. To run L2,
choose the L2 phantom in the viewer, edit the room and exposure, then use
File → Save geometry settings (.json). In Terminal, enter the extracted
project folder and pass the saved file to the launcher. For example, if you
named it my_L2_room.json in Downloads:

cd Mammography/radiation-protection   # the folder of this README
bash run_cases.sh "$HOME/Downloads/my_L2_room.json"

The default is one batch of 10,000 histories for the saved exposure. Terminal
prints the selected settings and the location of the results. See the project
README.md for batches, multiple saved files and full DBT sweeps.
