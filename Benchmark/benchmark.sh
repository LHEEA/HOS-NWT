cd Irreg_2D
if [ ! -d "Results" ]; then
 mkdir Results
fi
(../../bin/./HOS-NWT 1> /dev/null ; tec360 comp_Max.lay) &
cd ../Regular
if [ ! -d "Results" ]; then
 mkdir Results
fi
(../../bin/./HOS-NWT 1> /dev/null ; tec360 regular.lay) &
cd ../Sloshing
if [ ! -d "Results" ]; then
 mkdir Results
fi
(../../bin/./HOS-NWT 1> /dev/null ; tec360 sloshing.lay) &
