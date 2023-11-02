#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "../../allvars.h"
#include "mergers_io.h"

#ifdef OUTPUT_HOST_PROPERTIES_FOR_BH_MERGERS
void add_merger_event(int ThisTask, float Merger_Time, long long Merging_ID1 , MyFloat Merging_BH_Mass1, long long Merging_ID2, MyFloat Merging_BH_Mass2 , MyFloat Merging_BH_Mdot1 , MyFloat Merging_BH_Mdot2 , MyFloat Merging_BH_HostHaloMass1 , MyFloat Merging_BH_HostHaloMass2 , MyFloat Merging_BH_HostStellarMass1 , MyFloat Merging_BH_HostStellarMass2 , MyFloat Merging_BH_HostGasMass1 , MyFloat Merging_BH_HostGasMass2 , MyFloat Merging_BH_HostSFR1 , MyFloat Merging_BH_HostSFR2)
#else
void add_merger_event(int ThisTask, float Merger_Time, long long Merging_ID1 , MyFloat Merging_BH_Mass1, long long Merging_ID2, MyFloat Merging_BH_Mass2)
#endif
{
   memset(&MergerEvents[Nmergers], 0, sizeof(struct merger_properties));
   MergerEvents[Nmergers].ID1 = Merging_ID1; 
   MergerEvents[Nmergers].ID2 = Merging_ID2;
   MergerEvents[Nmergers].TaskID = ThisTask;
   MergerEvents[Nmergers].Time = Merger_Time;
   MergerEvents[Nmergers].BH_Mass1 = Merging_BH_Mass1;
   MergerEvents[Nmergers].BH_Mass2 = Merging_BH_Mass2;
#ifdef OUTPUT_HOST_PROPERTIES_FOR_BH_MERGERS
   MergerEvents[Nmergers].BH_Mdot1 = Merging_BH_Mdot1;
   MergerEvents[Nmergers].BH_Mdot2 = Merging_BH_Mdot2;
   MergerEvents[Nmergers].BH_HostHaloMass1 = Merging_BH_HostHaloMass1;
   MergerEvents[Nmergers].BH_HostHaloMass2 = Merging_BH_HostHaloMass2;
   MergerEvents[Nmergers].BH_HostStellarMass1 = Merging_BH_HostStellarMass1;
   MergerEvents[Nmergers].BH_HostStellarMass2 = Merging_BH_HostStellarMass2;
   MergerEvents[Nmergers].BH_HostGasMass1 = Merging_BH_HostGasMass1;
   MergerEvents[Nmergers].BH_HostGasMass2 = Merging_BH_HostGasMass2;
   MergerEvents[Nmergers].BH_HostSFR1 = Merging_BH_HostSFR1;
   MergerEvents[Nmergers].BH_HostSFR2 = Merging_BH_HostSFR2;
#endif
   Nmergers++;
}
