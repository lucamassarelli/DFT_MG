# DFT_MG
PROGRAM FOR ATOMIC CALCULATION WITH NUMBER OF ELECTRON FROM 0 TO 4

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or(at your option) any later version.

This program is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty ofMERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
See <http://www.gnu.org/licenses/>.

OPERATING SYSTEM:
	-WINDOWS (TESTED)
	-LINUX (TESTED)
	-MAC OS (NOT TESTED)

INSTALL:
	UNPACK THE ARCHIVE AND RUN FROM TERMINAL:
	make install clean
RUN:
	./DFT_MG

COMMAND LINE ARGUMENT:
-Z  ATOMIC NUMBER (default:4)(max:20) 
-S1 OCCUPATION NUMBER OF S1 ORBITAL (default: 2) (max:2)
-S2 OCCUPATION NUMBER OF S2 ORBITAL (default: 2) (max:2)
-C  TYPE OF CORREATION "wigner"or "cep-ald"(default) 
-Q  "filename" SAVE THE CHARGE DENSITY ON THE FILE NAME;
-H  SHOW THIS SCREEN

EXAMPLE
./DFT_MG -Z 4 -S1 2 -S2 2 -C wigner -Q data.txt

MADE BY
	-LUCA MASSARELLI 	lmass@hotmail.it
	-VALENTINA GREGORI	gregori.vale@gmail.com
	
