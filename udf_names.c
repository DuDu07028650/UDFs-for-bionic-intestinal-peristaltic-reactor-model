/* This file generated automatically. */
/*          Do not modify.            */
#include "udf.h"
#include "prop.h"
#include "dpm.h"
extern DEFINE_GRID_MOTION(segmentation_upperwall,domain,dt,time,dtime);
extern DEFINE_GRID_MOTION(segmentation_lowerwall,domain,dt,time,dtime);
extern DEFINE_SOURCE(enzyme_source,c,c_thread,dS,eqn);
extern DEFINE_EXCHANGE_PROPERTY(custom_drag,cell,mix_thread,s_col,f_col);
extern DEFINE_PROFILE(backflow_mass_fraction,t,i);
UDF_Data udf_data[] = {
{"segmentation_upperwall", (void (*)(void))segmentation_upperwall, UDF_TYPE_GRID_MOTION},
{"segmentation_lowerwall", (void (*)(void))segmentation_lowerwall, UDF_TYPE_GRID_MOTION},
{"enzyme_source", (void (*)(void))enzyme_source, UDF_TYPE_SOURCE},
{"custom_drag", (void (*)(void))custom_drag, UDF_TYPE_EXCHANGE_PROPERTY},
{"backflow_mass_fraction", (void (*)(void))backflow_mass_fraction, UDF_TYPE_PROFILE},
};
int n_udf_data = sizeof(udf_data)/sizeof(UDF_Data);
#include "version.h"
void UDF_Inquire_Release(int *major, int *minor, int *revision)
{
  *major = RampantReleaseMajor;
  *minor = RampantReleaseMinor;
  *revision = RampantReleaseRevision;
}
