# 31P-siml
## Short script(s) to quickly simulate a 31P NMR spectrum with isotropic, lamellar, hexagonal phases

In ```notebook``` you will find the quick and dirty jupyter notebook version. It if faster to work with if you know what you are doing.


For ```standard```:
Get the spectrum you want to simulate from TopSpin and save it as ```31P_INP``` in the standard folder.


Next change ```settings.txt``` to fit your spectrometer. Values for 300/600/750MHz are in there.

Now change the values in ```parameters.txt```. You can simulate up to two different components.
The parameters are:


`LB(axi/hex/iso)(1/2)`  :  this is the LB(line broadening) applied to the simulated FID of the phase of each component

`csa1/2`  :  CSA(chemical shift anisotropy of the component)

`c(1/2) and a(1/2)`  :  ellipsis axes ratio of your vesicle, sets the ratio of vesicle deformation/orientation within the magnetic field

`axi/hex/iso(1/2)`  :  percentage of the specific phase you have in the spectrum, axi=lamellar, hex=hecagonal, iso=isotropic

`drift(1/2)` :  drift values you can apply to the entire spectrum of each component

`plotdrift` :  drift value you can apply to your real data if it is shifted off-center

`isoshift(1/2)`  :  shift of the isotropic peaks if that has occured

`SHOW`  :  whether you want the matplotlib widget to pop out at every run of the program
