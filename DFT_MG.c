#include <stdio.h>
#include <stdlib.h>
#include "Hartree.h"
#include "rtwtypes.h"
#include "string.h"

/*This program is free software: you can redistribute it and/or modify it 
under the terms of the GNU General Public License as published by the Free 
Software Foundation, either version 3 of the License, or(at your option) 
any later version.

This program is distributed in the hope that it will be useful,but WITHOUT 
ANY WARRANTY;  without even the implied warranty of  MERCHANTABILITY or 
FITNESS FOR A PARTICULAR PURPOSE.   See the GNU General Public License 
for more details. See <http://www.gnu.org/licenses/>.*/



void dispHelp() {
  printf("---------------------------------------------------------------\n");
  printf("---ATOMS CALCULATION by Luca  Massarelli & Valentina Gregori---\n");
  printf("---------------------------------------------------------------\n"); 
  printf("This program run the calculation for atom from 0 to 4 electrons\n"); 
  printf("you can set seguent argument from command line:\n");
  printf("-Z  ATOMIC NUMBER (default:4)(max:20)\n"); 
  printf("-S1 OCCUPATION NUMBER OF S1 ORBITAL (default: 2) (max:2)n");
  printf("-S2 OCCUPATION NUMBER OF S2 ORBITAL (default: 2) (max:2)\n");
  printf("-C  TYPE OF CORREATION [\"wigner\"or \"cep-ald\"(default)]\n"); 
  printf("-Q  \"filename\" SAVE THE CHARGE DENSITY ON THE FILE NAME\n");
  printf("-H  SHOW THIS SCREEN\n");
}


                  

int main(int argc, char *argv[]) {

  int i;
  double temp;
  real_T z = 4;
  real_T f[5] = {2,2,0,0,0};
  char_T type[8] = "cep-ald";
  int32_T type_size[2] = {1,7}; 
  real_T r[500] = {0};
  real_T q[2500] = {0};
  real_T  autov;
  real_T eps[5];
  real_T eei;
  real_T eel;
  real_T exc;
  real_T dexc;
  real_T ecc;
  real_T dec;
  real_T ekin;
  real_T b_etotal;
  int fileFlag = 0;
  char filename[100] = {0};
  FILE* fp;

  for(i = 1; i < argc; i++) {                     //VALIDAZIONE INPUT
    if(strcmp(argv[i],"-Z") == 0) {
      temp = atof(argv[i+1]);
      if(temp >= 0 && temp <= 4){
	z = temp;
	i++;
      } else {
	printf("Z input not valid: 0<=z<=20\n");
	exit(0);
      }
    } else if(strcmp(argv[i],"-S1") == 0) {
        temp = atof(argv[i+1]);
      if(temp >= 0 && temp <= 2){
	  f[0] = temp;
	  i++;
      } else {
	  printf("S1 input not valid: 0<=S1<=2\n");
	  exit(0);
      }
    }else if(strcmp(argv[i],"-S2") == 0) {
	temp = atof(argv[i+1]);
	if(temp >= 0 && temp <= 2){
	  f[1] = temp;
	  i++;
	} else {
	  printf("S2 input not valid: 0<=S2<=4\n");
	  exit(0);
	}
    }else if(strcmp(argv[i],"-C") == 0) {
        if(strcmp(argv[i+1],"wigner") == 0){
	  strcpy(type,"wigner");
	  type_size[1] = 6;
	  i++;
        }else if(strcmp(argv[i+1],"cep-ald")) {
	  strcpy(type,"cep-ald");
	  type_size[1] = 7;
	  i++;
	}else {
	    printf("Correlation not valid: \"wigner\" or \"cep-ald\"");
	    exit(0);
	}
    }else if(strcmp(argv[i],"-Q") == 0) {
        if(strlen(argv[i+1]) > 100){
	  printf("filename too long");
	  exit(0);
	}else {
	  fileFlag = 1;
	  strcpy(filename,argv[i+1]);
	  i++;
	}
    }else {
       dispHelp();
       exit(0);
    }
    
  }
    

  printf("\n\n-----ATOMIC CALULATION BY Luca Massarelli & Valentina Gregori-----\n");
  printf("\nPass command line argument -H for help\n");
  printf("\n------------------------START CALCULATION-------------------------\n\n");
  printf("z = %2.2lf\n",z);
  printf("occupazione s1: %2.2lf\n",f[0]);
  printf("occupazione s2: %2.2lf\n",f[1]);
  printf("type = %s\n", type);

  Hartree_initialize();
  Hartree(z,f,type,type_size,&b_etotal,r,q, &autov,eps, &eei, &eel,&exc, &dexc, &ecc, &dec, &ekin);

  printf("\n--------------------------END CALCULATION-------------------------\n");
  printf("\n");

  printf("|total energy:           %6.4lf\n",b_etotal);
  printf("|eigenvalue 1s:          %6.4lf\n",eps[0]);
  printf("|eigenvalue 2s:          %6.4lf\n",eps[1]);
  printf("|sum eigenvalue:         %6.4lf\n",autov);
  printf("|electron-nucleo energy: %6.4lf\n",eei);
  printf("|electrostatic energy:   %6.4lf\n",eel);
  printf("|exchange energy:        %6.4lf\n",exc - ecc);
  printf("|correlation energy:     %6.4lf\n",ecc);
  printf("|correction exchange:    %6.4lf\n",dexc - dec);
  printf("|correction correlation: %6.4lf\n",dec); 

  printf("\n------------------------------------------------------------------\n");

  if(fileFlag == 1) {
    fp = fopen(filename,"w");
    if(fp < 0) {
      printf("ERROR ON OPENING FILE");
      exit(0);
    }
    fprintf(fp,"r,charge_density_1s, charge_density_2s, total_charge_density\n");
    for(i = 0; i < 500; i++) {
      fprintf(fp,"%lf,%lf,%lf,%lf \n",r[i],q[i],q[500+i],2*q[i]+ 2*q[500+i]);
    }
    fclose(fp);
    printf("created csv file: %s", filename);
  }

  printf("\n------------------------------------------------------------------\n");

}

