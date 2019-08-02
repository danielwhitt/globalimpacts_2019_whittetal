T=-2:.25:32;
S=20:.25:38;
[Tg,Sg]=meshgrid(T,S);
[RHOOUT,DRHODT,DRHODS] = mjwfstate(0,Tg,Sg);
ALPHA=-DRHODT./1026;
BETA=DRHODS./1026;

xidx=[318 197 281 382 33];
