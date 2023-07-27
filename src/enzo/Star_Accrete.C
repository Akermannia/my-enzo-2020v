/***********************************************************************
/
/  ADD ANY MASS MARKED FOR ACCRETION ONTO THE STAR PARTICLE
/
/  written by: John Wise
/  date:       March, 2009
/  modified1: Ji-hoon Kim
/             September, 2009
/
************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "phys_constants.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int Star::Accrete(void)
{
   if (this->CurrentGrid == NULL || 
      (this->naccretions == 0 && fabs(this->DeltaMass) < tiny_number))
    return SUCCESS;
  

  int dim, i, n, count;
  FLOAT time = CurrentGrid->Time;
  float dt = CurrentGrid->dtFixed;
  float this_dt;
  
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
    VelocityUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, time);

  /* Sum up prescribed mass accretion, then add to mass */

  n = 0;
  while (n < naccretions && accretion_time[n] <= time+dt) {
    if (n+1 == naccretions)
      this_dt = time - accretion_time[n];
    else
      this_dt = accretion_time[n+1] - accretion_time[n];
    DeltaMass += accretion_rate[n++] * this_dt * TimeUnits;
  }

//  printf("star::Accrete: old_Mass = %lf, DeltaMass = %f\n", Mass, DeltaMass); 
  double old_mass = Mass;
  Mass += (double)(DeltaMass);
//  FinalMass += (double)(DeltaMass);
//  printf("star::Accrete: new_Mass = %lf, DeltaMass = %f\n", Mass, DeltaMass); 


  /* Conserve momentum: change star particle velocity due to accreted material */

  /* [1] For BlackHole, 
     it is now accurately done in Star_SubtractAccretedMassFromCell */

  /* [2] For Star Formation, 
     We can still do this in approximate way (JHW, Jan10) */

  double ratio1, ratio2, new_vel;

  if (type != MBH || type != BlackHole) {
    ratio2 = DeltaMass / Mass;
    ratio1 = 1.0 - ratio2;
    Metallicity = ratio1 * Metallicity + ratio2 * deltaZ;
    deltaZ = 0.0;
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      vel[dim] = ratio1 * vel[dim] + ratio2 * delta_vel[dim];
      delta_vel[dim] = 0.0;
    }
  } else {
    deltaZ = 0.0;
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      delta_vel[dim] = 0.0;
    }
  }

  /* [3] For MBH, 
     because gas mass is added to MBH from many cells with zero net momentum,
     just decrease the particle's velocity accordingly. */

  if (type == MBH) {
    const int particle_count_threshold = 20;
    int particle_count = 0;
    float particle_velocity[MAX_DIMENSION] = {0.};
    const float particle_damping_timescale = 0.001f;
    const float acceleration_damping_timescale = 0.0001f;
    const float alpha = -std::expm1(-this_dt / particle_damping_timescale);
    const float beta = -std::expm1(-this_dt / acceleration_damping_timescale);
    const float distance2 = pow(7 * CurrentGrid->GetCellWidth(0, 0), 2.0f);

    for (int i = 0; i < CurrentGrid->NumberOfParticles; ++i)
    {
      if (CurrentGrid->ParticleType[i] != PARTICLE_TYPE_DARK_MATTER &&
          CurrentGrid->ParticleType[i] != PARTICLE_TYPE_STAR)
          continue;
      float dx = CurrentGrid->ParticlePosition[0][i] - pos[0];
      float dy = CurrentGrid->ParticlePosition[1][i] - pos[1];
      float dz = CurrentGrid->ParticlePosition[2][i] - pos[2];
      if (dx * dx + dy * dy + dz * dz > distance2)
        continue;
      particle_velocity[0] += CurrentGrid->ParticleVelocity[0][i];
      particle_velocity[1] += CurrentGrid->ParticleVelocity[1][i];
      particle_velocity[2] += CurrentGrid->ParticleVelocity[2][i];
      ++particle_count;
    }
    for (size_t d = 0; d < MAX_DIMENSION; ++d) {
      particle_velocity[d] /= max(1.0f, float(particle_count));
    }

    float acceleration[MAX_DIMENSION] = {0.};
    float acceleration_magnitude = 0.0f;
    for (size_t d = 0; d < MAX_DIMENSION; ++d) {
      acceleration[d] = vel[d] - last_vel[d];
      acceleration_magnitude += acceleration[d] * acceleration[d];
    }

    if (acceleration_magnitude > 0.0f) {
      acceleration_magnitude = std::sqrt(acceleration_magnitude);
      for (size_t d = 0; d < MAX_DIMENSION; ++d) {
        acceleration[d] /= acceleration_magnitude;
      }
      float inner_acceration_velocity = 0.0f;
      for (size_t d = 0; d < MAX_DIMENSION; ++d) {
        inner_acceration_velocity += acceleration[d] * vel[d];
      }

      // damp the velocity perpendicular to the acceleration
      for (size_t d = 0; d < MAX_DIMENSION; ++d) {
        vel[d] = (1.f - beta) * vel[d] +
          beta * (inner_acceration_velocity * acceleration[d]);
      }
    }

    // damp the velocity relative to surrounding particles
    if (particle_count >= particle_count_threshold) {
      for (size_t d = 0; d < MAX_DIMENSION; ++d) {
        vel[d] = (1.f - alpha) * vel[d] + alpha * particle_velocity[d];
      }
    }

    printf("StarAccrete[%"ISYM"]: n_drag = %"ISYM", vel = (%"GSYM", %"GSYM", %"GSYM"), pos = (%"GSYM", %"GSYM", %"GSYM"), mass = %"GSYM"\n",
      Identifier, particle_count, vel[0], vel[1], vel[2], pos[0], pos[1], pos[2], FLOAT(Mass));

    // Keep the last velocity for future use
    for (size_t d = 0; d < MAX_DIMENSION; ++d) {
      last_vel[d] = vel[d];
    }
  }


  /* Keep the last accretion_rate for future use */

  if (n > 0)  last_accretion_rate = accretion_rate[n-1]; 

  //fprintf(stdout, "star::Accrete:  last_accretion_rate = %"GOUTSYM
	//  " SolarMass/yr, time = %"GOUTSYM", "
	//  "accretion_time[0] = %"GOUTSYM", this_dt = %"GOUTSYM
	//  ", DeltaMass = %"GOUTSYM", Mass = %lf\n",
	//  last_accretion_rate*yr_s, time, accretion_time[0], this_dt, DeltaMass, Mass);

  /* Remove these entries in the accretion table */

  count = 0;
  naccretions -= n;
  if (naccretions > 0) {
    FLOAT *temp_time = new FLOAT[naccretions];
    float *temp_rate = new float[naccretions];
    for (i = 0; i < naccretions+n; i++)
      if (accretion_time[i] <= time+dt) {
	temp_rate[count] = accretion_rate[i];
	temp_time[count] = accretion_time[i];
	count++;
      }
  
    delete [] accretion_rate;
    delete [] accretion_time;
    accretion_rate = temp_rate;
    accretion_time = temp_time;
  } else {  
    // No more accretion data
    accretion_rate = NULL;
    accretion_time = NULL;
  }
  
  return SUCCESS;
}

