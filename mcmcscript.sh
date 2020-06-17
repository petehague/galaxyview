cd v1-1-0-0.3-0
for halo in `seq 1 12`
do
  dir=v1-$halo-0-0.3-0
  cd ../$dir
  $1/RainfallMCMC/rainfall $1/AlphaBetaGammaHalo.ini data=curvedata.txt path=obscurve
  $1/RainfallMCMC/rainfall $1/AlphaBetaGammaHalo.ini data=realcurvedata.txt path=realcurve
  dir=v1-$halo-0-0.3-0.5
  cd ../$dir
  $1/RainfallMCMC/rainfall $1/AlphaBetaGammaHalo.ini data=curvedata.txt path=obscurve
  $1/RainfallMCMC/rainfall $1/AlphaBetaGammaHalo.ini data=realcurvedata.txt path=realcurve
done

cd ..
