#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define MAXLEN 256

int main(int argc, char ** argv) {
  FILE * fin;
  int ngx, ngy, ngz;
  int ndata;
  double *v;
  int nskip;

  char line[MAXLEN];
  char * pch;
  int ii, ix, iy;
  double avg, t1, t2, t3, lat_c;

  fin=fopen(argv[1], "r");
  printf("# Reading data from file %s\n", argv[1]);

  for(ii=0; ii<4; ii++) {
    fgets(line, MAXLEN, fin);
  }
  fgets(line, MAXLEN, fin);
  sscanf(line, " %lf %lf %lf", &t1, &t2, &t3);
  lat_c=sqrt(t1*t1+t2*t2+t3*t3);
  fgets(line, MAXLEN, fin);
  fgets(line, MAXLEN, fin);
  pch=strtok(line, " ");
  while (pch) {
    sscanf(pch, " %d", &ii);
    nskip+=ii;
    pch = strtok(NULL, " ");
  }

  printf("#  Skipping %5d lines...\n", nskip);

  for(ii=0; ii<nskip+2; ii++) {
    fgets(line, MAXLEN, fin);
  }

  fgets(line, MAXLEN, fin);
  sscanf(line, " %d %d %d", &ngx, &ngy, &ngz);
  v=(double *)malloc(sizeof(double)*ngx*ngy*ngz);

  printf("# Data dimension: %4dx%4dx%4d\n", ngx, ngy, ngz);
  for(ii=0; ii<ngx*ngy*ngz; ii++) {
    if (ii%5==0) {
      fgets(line, MAXLEN, fin);
      pch=line;
    }
    sscanf(pch, " %lf", v+ii);
    pch+=18;
  }

  for(ii=0; ii<ngz; ii++) {
    avg=0.0;
    for(ix=0; ix<ngx; ix++) {
      for(iy=0; iy<ngy; iy++) {
        avg+=v[ii*ngx*ngy+iy*ngx+ix];
      }
    }
    printf("%14.9f%22.16f\n", (ii*lat_c)/ngz, avg/(ngx*ngy));
  }
  avg=0.0;
  for(ix=0; ix<ngx; ix++) {
    for(iy=0; iy<ngy; iy++) {
      avg+=v[iy*ngx+ix];
    }
  }
  printf("%14.9f%22.16f\n", lat_c, avg/(ngx*ngy));


  free(v);
  return 0;
}
