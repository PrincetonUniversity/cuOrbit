#ifndef SET_ORBIT_CONFIG_H_
#define SET_ORBIT_CONFIG_H_

/* these are set in the .c file for now */
extern const int IDP;
extern const int IDT;
extern const int NTOR;

typedef struct Config Config_t;

void initialize_Config(Config_t*);


#endif
