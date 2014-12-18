cd Irreg_3D
if [ ! -d "Results" ]; then
 mkdir Results
fi
(../../bin/./HOS-NWT 1> /dev/null ; tec360 comp_3D.lay) &
cd Regular_3D
if [ ! -d "Results" ]; then
 mkdir Results
fi
(../../bin/./HOS-NWT 1> /dev/null ; tec360 regular_3D.lay) &
