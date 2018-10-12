#define PH_PI 3.1415927
#define PH_TWOPI 6.283185308
double phrng(double angle) // return in range (0,2pi)
{
  double val=angle;
  while (val<0.0)      val += PH_TWOPI;
  while (val>PH_TWOPI) val -= PH_TWOPI;
  return(val);
}
/* --------------------------------- */
double delphrng(double angle) // return in range (-2pi,2pi )
{
  double val=angle;
  while (val<-PH_PI) val += PH_TWOPI;
  while (val> PH_PI) val -= PH_TWOPI;
  return(val);
}
