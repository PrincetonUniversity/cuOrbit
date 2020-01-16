/*  Copyright 2019, 2020 Garrett Wright, Princeton Plasma Physic Lab

    This file is part of CuOrbit.

    CuOrbit is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CuOrbit is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CuOrbit.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "orbit_config_api.h"

// test
int main(int argc, char **argv){
  char default_config_file[] = "INPUT/config.ini";
  char* config_file = default_config_file;

  if(argc > 2) {
    printf("\nError, too many arguments"
           "\nusage: %s <optional_config_file>\n"
           "\t Defaults to INPUT/config.ini\n\n", argv[0]);
  }
  if(argc == 2){
    config_file = argv[1];
    printf("User supplied config file path %s\n", config_file);
  }

  /* Contruct a config object,
   which conttructs Perturb, Equlib, and Particles */
  orbit_Config_t* Cfg=orbit_Config_ctor();

  printf("Initialize Config\n");
  orbit_initialize_Config(Cfg, config_file);

  orbit_main_loop(Cfg);

  return 0;
}
