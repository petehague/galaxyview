cd v1-1-0-0.3-0
for halo in `seq 1 12`
do
  dir=v1-$halo-0-0.3-0
  cd ../$dir
  ~/ApostleView/pollux -j -c
  ~/ApostleView/pvplot.r
  ~/Durham/RainfallMCMC/rainfall ~/Durham/ABGHalo/AlphaBetaGammaHalo.ini data=curvedata.txt path=obscurve
  ~/Durham/RainfallMCMC/rainfall ~/Durham/ABGHalo/AlphaBetaGammaHalo.ini data=realcurvedata.txt path=realcurve
done

for halo in `seq 1 12`
do
  dir=v1-$halo-0-0.3-0.5
  cd ../$dir
  ~/ApostleView/pollux -j -c
  ~/ApostleView/pvplot.r
  ~/Durham/RainfallMCMC/rainfall ~/Durham/ABGHalo/AlphaBetaGammaHalo.ini data=curvedata.txt path=obscurve
done

cd ..
~/ApostleView/maketable.r /Volumes/Workspace/ApostlesHI/surveymrhi
