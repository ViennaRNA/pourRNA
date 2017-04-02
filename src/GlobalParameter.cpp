/*
 * GlobalParameter.cpp
 *
 *  Created on: 12.08.2014
 *      Author: Gregor Entzian
 */

#include "GlobalParameter.h"

GlobalParameter::GlobalParameter (double temperature)
{
  model_detailsT md;
  set_model_details (&md);
  EnergyParameter = get_scaled_parameters (temperature, md);
}

GlobalParameter * GlobalParameter::Instance = NULL;

GlobalParameter *
GlobalParameter::getInstance ()
{
  if (Instance == NULL)
    {
      Instance = new GlobalParameter (37);
    }
  return (Instance);
}

void
GlobalParameter::Initialize (double temperature)
{
  if (Instance != NULL)
    {
      delete Instance;
      Instance = NULL;
    }
  Instance = new GlobalParameter (temperature);
}

GlobalParameter::~GlobalParameter ()
{
  if (EnergyParameter != NULL)
    {
      free (EnergyParameter);
      EnergyParameter = NULL;
    }
}
