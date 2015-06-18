/*
 * Hartree.c
 *
 * Code generation for function 'Hartree'
 *
 * C source code generated on: Sun Jun 07 18:43:36 2015
 *
 */

/* Include files */
#include <stdio.h>
#include "rt_nonfinite.h"
#include "Hartree.h"

/* Function Declarations */
static void Griglia(real_T r1, real_T r[500], real_T nx);
static void Schroed(real_T n, real_T *eps, real_T z, const real_T r[500], const
                    real_T v[500], real_T mmax, real_T u[500]);
static real_T ThomasFermi(real_T r, real_T z);
static boolean_T b_eml_strcmp(const char_T a_data[20], const int32_T a_size[2]);
static real_T dxc1(real_T x, const char_T type_data[20], const int32_T
                   type_size[2]);
static void eml_li_find(const boolean_T x[500], int32_T y_data[500], int32_T
  y_size[1]);
static boolean_T eml_strcmp(const char_T a_data[20], const int32_T a_size[2]);
static void etotal(const real_T r[500], const real_T rhot[500], const real_T
                   vhe[500], const real_T vio[500], const real_T vxc[500], const
                   real_T vc[500], const real_T dxc[500], const real_T dc[500],
                   const real_T e[5], const real_T f[5], real_T nx, real_T *etot,
                   real_T *exc, real_T *ec, real_T *ekin, real_T *eel, real_T
                   *eei, real_T *dec, real_T *dexc);
static real_T fxc1(real_T x, const char_T type_data[20], const int32_T
                   type_size[2]);
static void poisson(const real_T r[500], real_T z, real_T r0, const real_T coef
                    [5], const real_T chi[2500], real_T nx, const real_T eps[5],
                    const real_T f[5], const char_T type_data[20], const int32_T
                    type_size[2], real_T vtot[500], real_T *etot, real_T *exc,
                    real_T *ec, real_T *ekin, real_T *eel, real_T *eei, real_T
                    *dec, real_T *dexc);
static real_T rdivide(real_T x, real_T y);
static real_T rt_powd_snf(real_T u0, real_T u1);

/* Function Definitions */
static void Griglia(real_T r1, real_T r[500], real_T nx)
{
  int32_T i;
  r[0] = r1;

  /* primo punto della griglia */
  for (i = 0; i < (int32_T)(nx + -1.0); i++) {
    /* ciclo: costruzione griglia */
    r[i + 1] = 1.05 * r[(int32_T)((2.0 + (real_T)i) - 1.0) - 1];
  }
}

static void Schroed(real_T n, real_T *eps, real_T z, const real_T r[500], const
                    real_T v[500], real_T mmax, real_T u[500])
{
  real_T up[500];
  real_T upp[500];
  real_T cf[500];
  int32_T i;
  real_T mch;
  real_T amesh;
  real_T al;
  real_T als;
  real_T emax;
  real_T emin;
  real_T node;
  int32_T iteraz;
  int32_T exitg1;
  boolean_T bv2[500];
  int32_T tmp_size[1];
  int32_T tmp_data[500];
  int32_T loop_ub;
  int32_T i0;
  real_T c3;
  real_T sn;
  real_T uout;
  real_T upout;
  real_T nin;
  for (i = 0; i < 500; i++) {
    u[i] = 0.0;
    up[i] = 0.0;
    upp[i] = 0.0;
    cf[i] = 0.0;
  }

  mch = 0.0;

  /* parametri griglia, errore relativo energia, l(l+1), energia max e min */
  amesh = r[1] / r[0];
  al = log(amesh);
  als = al * al;
  emax = v[(int32_T)mmax - 1] + 0.0 / (r[(int32_T)mmax - 1] * r[(int32_T)mmax -
    1]);
  emin = -0.0;
  for (i = 0; i < (int32_T)mmax; i++) {
    node = v[i] + 0.0 / (r[i] * r[i]);
    if ((emin <= node) || rtIsNaN(node)) {
    } else {
      emin = node;
    }
  }

  /* riporta l'energia (negativa) dell'input nella finestra emin-emax */
  if (*eps > emax) {
    *eps = 1.25 * emax;
  }

  if (*eps < emin) {
    *eps = 0.75 * emin;
  }

  if (*eps > emax) {
    *eps = 0.5 * (emax + emin);
  }

  /*  azzera il contatore delle iterazioni */
  iteraz = 0;

  /* *********************************************************************** */
  /* *** FINE FASE A ******************************************************* */
  /*  */
  /* *** FASE B: iterazioni successive che migliorano eps rispetto al valore */
  /* *** di input fino alla precisione voluta******************************* */
  /* *********************************************************************** */
  do {
    exitg1 = 0;
    iteraz++;

    /*  iteraz conta iterazioni */
    /*  stop (con messaggio) se si superano 100 iterazioni */
    if (iteraz > 100) {
      exitg1 = 1;
    } else {
      /*  definisci array coefficienti nell'eq. differ. per u (cioe' per chi) */
      for (i = 0; i < (int32_T)mmax; i++) {
        cf[i] = als * 0.0 + 2.0 * als * (v[i] - *eps) * (r[i] * r[i]);
      }

      for (i = 0; i < 500; i++) {
        bv2[i] = rtIsNaN(cf[i]);
      }

      eml_li_find(bv2, tmp_data, tmp_size);
      loop_ub = tmp_size[0];
      for (i0 = 0; i0 < loop_ub; i0++) {
        cf[tmp_data[i0] - 1] = 0.0;
      }

      /* trova il punto d'inversione classica sapendo l'array cf e l'energia eps, */
      /*  e' il punto di raccordo della fz. d'onda radiale (cioe' l'indice mch) */
      i = (int32_T)mmax;
      while (((cf[i - 2] > 0.0) || (cf[i - 1] <= 0.0)) && (i > 2)) {
        mch = i;
        i--;
      }

      mch--;

      /* stampa messaggio d'errore se non si trova il punto */
      if (i == 3) {
        exitg1 = 1;
      } else {
        /* c>>> FASE B1: integrazione da r=0 fino al punto di raccordo r(mch)>>>>>> */
        /*  funzione d'onda a piccoli r (i=1,4): serie di taylor */
        c3 = z * z / 3.0 - *eps / 3.0;
        for (i = 0; i < 4; i++) {
          node = (r[i] + -z * (r[i] * r[i])) + c3 * rt_powd_snf(r[i], 3.0);
          sn = al * ((r[i] + -z * 2.0 * (r[i] * r[i])) + c3 * 3.0 * rt_powd_snf
                     (r[i], 3.0));
          upp[i] = al * sn + cf[i] * node;
          u[i] = node;
          up[i] = sn;
        }

        /*  integrazione verso l'esterno, cioe' da r=r(5) fino al raggio r(mch) */
        /*  usando predictor una volta, corrector due volte */
        node = 0.0;
        for (i = 0; i < (int32_T)((mch - 1.0) + -3.0); i++) {
          /* adams extrapolation outward */
          u[(int32_T)((4.0 + (real_T)i) + 1.0) - 1] = u[i + 3] + 0.041666667 *
            (((55.0 * up[(int32_T)(4.0 + (real_T)i) - 1] - 59.0 * up[(int32_T)
               ((4.0 + (real_T)i) - 1.0) - 1]) + 37.0 * up[(int32_T)((4.0 +
                (real_T)i) - 2.0) - 1]) - 9.0 * up[(int32_T)((4.0 + (real_T)i) -
              3.0) - 1]);

          /* adams extrapolation outward */
          up[(int32_T)((4.0 + (real_T)i) + 1.0) - 1] = up[(int32_T)(4.0 +
            (real_T)i) - 1] + 0.041666667 * (((55.0 * upp[(int32_T)(4.0 +
            (real_T)i) - 1] - 59.0 * upp[(int32_T)((4.0 + (real_T)i) - 1.0) - 1])
            + 37.0 * upp[(int32_T)((4.0 + (real_T)i) - 2.0) - 1]) - 9.0 * upp
            [(int32_T)((4.0 + (real_T)i) - 3.0) - 1]);
          for (loop_ub = 0; loop_ub < 2; loop_ub++) {
            upp[(int32_T)((4.0 + (real_T)i) + 1.0) - 1] = al * up[(int32_T)((4.0
              + (real_T)i) + 1.0) - 1] + cf[(int32_T)((4.0 + (real_T)i) + 1.0) -
              1] * u[(int32_T)((4.0 + (real_T)i) + 1.0) - 1];

            /* adams interpolation outward */
            up[(int32_T)((4.0 + (real_T)i) + 1.0) - 1] = up[(int32_T)(4.0 +
              (real_T)i) - 1] + 0.041666667 * (((9.0 * upp[(int32_T)((4.0 +
              (real_T)i) + 1.0) - 1] + 19.0 * upp[(int32_T)(4.0 + (real_T)i) - 1])
              - 5.0 * upp[(int32_T)((4.0 + (real_T)i) - 1.0) - 1]) + upp
              [(int32_T)((4.0 + (real_T)i) - 2.0) - 1]);

            /* adams interpolation outward */
            u[(int32_T)((4.0 + (real_T)i) + 1.0) - 1] = u[(int32_T)(4.0 +
              (real_T)i) - 1] + 0.041666667 * (((9.0 * up[(int32_T)((4.0 +
              (real_T)i) + 1.0) - 1] + 19.0 * up[(int32_T)(4.0 + (real_T)i) - 1])
              - 5.0 * up[(int32_T)((4.0 + (real_T)i) - 1.0) - 1]) + up[(int32_T)
              ((4.0 + (real_T)i) - 2.0) - 1]);
          }

          if (u[(int32_T)((4.0 + (real_T)i) + 1.0) - 1] * u[(int32_T)(4.0 +
               (real_T)i) - 1] <= 0.0) {
            node++;
          }
        }

        if ((node - n) + 1.0 < 0.0) {
          /* pochi nodi: inutile andare avanti, provare subito nuova energia */
          emin = *eps;
          *eps *= 0.75;
          if (*eps > emax) {
            *eps = 0.5 * (emin + emax);
          }
        } else if ((node - n) + 1.0 > 0.0) {
          /*  troppi nodi: inutile andare avanti, provare subito nuova energia */
          emax = *eps;
          *eps *= 1.25;
          if (*eps < emin) {
            *eps = 0.5 * (emin + emax);
          }
        } else {
          /*  numero di nodi giusto: ok, procedere all'integrazione verso l'interno */
          uout = u[(int32_T)mch - 1];
          upout = up[(int32_T)mch - 1];

          /*  >>> FINE FASE B1 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
          /*  >>> FASE B2: integrazione da "infinito" al punto di raccordo r(mch)>>>> */
          /*  ultimo punto r(nin): 10 volte raggio del punto d'inversione classico */
          nin = floor(mch + 2.3 / al);
          if (nin + 4.0 > mmax) {
            nin = mmax - 4.0;
          }

          /*  funzione d'onda a grandi r (ultimi 4 punti): esponenziale semplice */
          node = v[(int32_T)nin - 1] - *eps;
          node = sqrt(0.0 / (r[(int32_T)nin - 1] * r[(int32_T)nin - 1]) + 2.0 *
                      node);
          i0 = (int32_T)(((real32_T)nin + 4.0F) + (1.0F - (real32_T)nin));
          for (loop_ub = -1; loop_ub + 1 < i0; loop_ub++) {
            i = (int32_T)nin + loop_ub;
            u[i] = exp(-node * (r[i] - r[(int32_T)nin - 1]));
            up[i] = -r[i] * al * node * u[i];
            upp[i] = al * up[i] + cf[i] * u[i];
          }

          /*  integrazione verso l'interno, cioe' da un raggio grande r(nin) verso */
          /*  raggi piu' piccoli, fino al raggio di raccordo r(mch) */
          for (i = (int32_T)nin - 2; i + 2 > (int32_T)mch; i--) {
            /* adams extrapolation inward */
            u[i] = u[i + 1] + -0.041666667 * (((55.0 * up[i + 1] - 59.0 * up[i +
              2]) + 37.0 * up[i + 3]) - 9.0 * up[i + 4]);

            /* adams extrapolation inward */
            up[i] = up[i + 1] + -0.041666667 * (((55.0 * upp[i + 1] - 59.0 *
              upp[i + 2]) + 37.0 * upp[i + 3]) - 9.0 * upp[i + 4]);
            for (loop_ub = 0; loop_ub < 2; loop_ub++) {
              upp[i] = al * up[i] + cf[i] * u[i];

              /* adams interpolation inward */
              up[i] = up[i + 1] + -0.041666667 * (((9.0 * upp[i] + 19.0 * upp[i
                + 1]) - 5.0 * upp[i + 2]) + upp[i + 3]);

              /* adams interpolation inward */
              u[i] = u[i + 1] + -0.041666667 * (((9.0 * up[i] + 19.0 * up[i + 1])
                - 5.0 * up[i + 2]) + up[i + 3]);
            }
          }

          /* >>> FINE FASE B2 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
          /* >>> FASE B3: raccordo la funzione in r(mch) e ne calcolo la norma>>>>>>       */
          /*  raccordo con continuita': riscalo u(i>mch) */
          node = uout / u[(int32_T)mch - 1];
          i0 = (int32_T)((real32_T)nin + (1.0F - (real32_T)mch));
          for (loop_ub = -1; loop_ub + 1 < i0; loop_ub++) {
            i = (int32_T)mch + loop_ub;
            up[i] *= node;
            u[i] *= node;
          }

          /*  calcolo norma */
          node = r[0] / sqrt(amesh);
          sn = (rt_powd_snf(node, 3.0) / 3.0 + 2.0 * -z * rt_powd_snf(node, 4.0)
                / 4.0) + (2.0 * c3 + -z * -z) * rt_powd_snf(node, 5.0) / 5.0;
          for (i = 0; i < (int32_T)nin; i++) {
            sn += al * r[i] * (u[i] * u[i]);
          }

          /*  >>> FINE FASE B3 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
          /*  >>> FINE FASE B3 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
          /*  >>> FASE B4: valutazione errore energia e decisione se andare avanti>>> */
          /*  teoria perturbativa per energia della prossima iterazione********* */
          /*  basata su relazione fra mismatch derivata e errore energia******** */
          node = 0.5 * uout * (upout - up[(int32_T)mch - 1]) / (sn * al * r
            [(int32_T)mch - 1]);

          /*  se l'errore relativo fra questa e la precedente energia e' minore di */
          /*  relerr (vedi input), normalizza u ed esci dalla subroutine (FASE C) */
          /*  se invece l'errore relativo rispetto alla precedente energia e' */
          /*  maggiore o uguale di relerr (vedi input), correggi l'energia e fai */
          /*  un'altra iterazione (cioe' torna all'inizio della FASE B) */
          if (fabs(node / *eps) < 1.0E-5) {
            /* >>> FINE FASE B3 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
            /* *********************************************************************** */
            /* *** FINE FASE B ******************************************************* */
            /* *** FASE C: normalizza funzione, pulisci array, metti nella variabile */
            /* *** relerr l'effettivo errore relativo ed esci dalla subroutine******** */
            /* *********************************************************************** */
            /*  normalizza la funzione d'onda */
            node = 1.0 / sqrt(sn);
            for (i = 0; i < (int32_T)nin; i++) {
              u[i] *= node;
            }

            /*  pulizia: azzera la funzione d'onda al di la' di r(nin) */
            i0 = (int32_T)((real32_T)mmax + (1.0F - ((real32_T)nin + 1.0F)));
            for (i = 0; i < i0; i++) {
              u[(int32_T)nin + i] = 0.0;
            }

            /*  metti in relerr (output) l'effettivo errore relativo */
            /*  esci dalla subroutine */
            /* ritorna nel programma principale */
            /* *********************************************************************** */
            /* *** FINE FASE c ******************************************************* */
            exitg1 = 1;
          } else {
            if (fabs(node / *eps) > 0.25) {
              node = 0.25 * node * fabs(*eps / node);
            }

            if (node > 0.0) {
              emin = *eps;
            }

            if (node < 0.0) {
              emax = *eps;
            }

            *eps += node;
            if ((*eps > emax) || (*eps < emin)) {
              *eps = 0.5 * (emax + emin);
            }
          }
        }
      }
    }
  } while (exitg1 == 0);
}

static real_T ThomasFermi(real_T r, real_T z)
{
  real_T x;

  /*  tfapot */
  /*  generalized thomas fermi atomic potential */
  x = r / rt_powd_snf(0.69395656 / z, 0.33333333333333331);
  x = z / ((1.0 + sqrt(x) * (0.02747 - x * (0.1486 - 0.007298 * x))) + x *
           (1.243 + x * (0.2302 + 0.006944 * x)));
  if (x < 1.0) {
    x = 1.0;
  }

  return -x / r;
}

static boolean_T b_eml_strcmp(const char_T a_data[20], const int32_T a_size[2])
{
  boolean_T b_bool;
  int32_T k;
  int32_T exitg2;
  int32_T exitg1;
  static const char_T cv1[6] = { 'w', 'i', 'g', 'n', 'e', 'r' };

  b_bool = FALSE;
  k = 0;
  do {
    exitg2 = 0;
    if (k < 2) {
      if (a_size[k] != 1 + 5 * k) {
        exitg2 = 1;
      } else {
        k++;
      }
    } else {
      k = 0;
      exitg2 = 2;
    }
  } while (exitg2 == 0);

  if (exitg2 == 1) {
  } else {
    do {
      exitg1 = 0;
      if (k <= a_size[1] - 1) {
        if (a_data[k] != cv1[k]) {
          exitg1 = 1;
        } else {
          k++;
        }
      } else {
        b_bool = TRUE;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }

  return b_bool;
}

static real_T dxc1(real_T x, const char_T type_data[20], const int32_T
                   type_size[2])
{
  real_T b_dxc1;
  real_T xs;
  real_T rs;
  real_T ec;
  b_dxc1 = 0.0;
  if ((!eml_strcmp(type_data, type_size)) && (!b_eml_strcmp(type_data, type_size)))
  {
  } else {
    /*  correzione energia scambio e correlazione (wigner o ceperley-alder) */
    /*  x=densita**1/3 in a. u. */
    xs = 12.566370614359172 * x;
    b_dxc1 = rt_powd_snf(x, 4.0) * (0.24618625546067416 + 0.2359 / ((1.0 + xs) *
                                     (1.0 + xs))) * 12.566370614359172;
    if (b_eml_strcmp(type_data, type_size)) {
    } else {
      /*  ceperley correlation */
      if (x < 0.0001) {
        b_dxc1 = 0.0;
      } else {
        rs = 0.6203504908994 / x;
        xs = sqrt(rs);
        if (rs < 1.0) {
          xs = log(rs);
          ec = ((0.0311 * xs + -0.048) + 0.002 * rs * xs) + -0.0116 * rs;
          xs = ((0.0311 * xs + -0.058366666666666664) + 0.004 * rs * xs / 3.0) +
            -0.0252 * rs / 3.0;
        } else {
          ec = -0.1423 / ((1.0 + 1.0529 * xs) + 0.3334 * rs);
          xs = ec * ((1.0 + 1.2283833333333334 * xs) + 0.44453333333333328 * rs)
            / ((1.0 + 1.0529 * xs) + 0.3334 * rs);
        }

        b_dxc1 = (0.24618625546067416 * rt_powd_snf(x, 4.0) + (ec - xs) *
                  rt_powd_snf(x, 3.0)) * 12.566370614359172;
      }
    }
  }

  return b_dxc1;
}

static void eml_li_find(const boolean_T x[500], int32_T y_data[500], int32_T
  y_size[1])
{
  int32_T k;
  int32_T i;
  k = 0;
  for (i = 0; i < 500; i++) {
    if (x[i]) {
      k++;
    }
  }

  y_size[0] = k;
  k = 0;
  for (i = 0; i < 500; i++) {
    if (x[i]) {
      y_data[k] = i + 1;
      k++;
    }
  }
}

static boolean_T eml_strcmp(const char_T a_data[20], const int32_T a_size[2])
{
  boolean_T b_bool;
  int32_T k;
  int32_T exitg2;
  int32_T exitg1;
  static const char_T cv0[7] = { 'c', 'e', 'p', '-', 'a', 'l', 'd' };

  b_bool = FALSE;
  k = 0;
  do {
    exitg2 = 0;
    if (k < 2) {
      if (a_size[k] != 1 + 6 * k) {
        exitg2 = 1;
      } else {
        k++;
      }
    } else {
      k = 0;
      exitg2 = 2;
    }
  } while (exitg2 == 0);

  if (exitg2 == 1) {
  } else {
    do {
      exitg1 = 0;
      if (k <= a_size[1] - 1) {
        if (a_data[k] != cv0[k]) {
          exitg1 = 1;
        } else {
          k++;
        }
      } else {
        b_bool = TRUE;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }

  return b_bool;
}

static void etotal(const real_T r[500], const real_T rhot[500], const real_T
                   vhe[500], const real_T vio[500], const real_T vxc[500], const
                   real_T vc[500], const real_T dxc[500], const real_T dc[500],
                   const real_T e[5], const real_T f[5], real_T nx, real_T *etot,
                   real_T *exc, real_T *ec, real_T *ekin, real_T *eel, real_T
                   *eei, real_T *dec, real_T *dexc)
{
  real_T r0;
  real_T corr;
  real_T rhot0;
  int32_T i;

  /*  rhot=4*pi*densita, vio=-Z/r, vhe=pot.hartree, e=autoval, f=occupaz */
  /* parametri della griglia radiale */
  r0 = r[0] / 1.02469507659596;
  corr = (r0 - r[0]) / (r[1] - r[0]);
  rhot0 = (rhot[1] - rhot[0]) * corr;
  r0 = rt_powd_snf(r0, 3.0) / 3.0;
  *eel = 0.5 * rhot0 * ((vhe[1] - vhe[0]) * corr) * r0;
  *eei = rhot0 * ((vio[1] - vio[0]) * corr) * r0;
  *exc = rhot0 * ((vxc[1] - vxc[0]) * corr) * r0;
  *dexc = (dxc[1] - dxc[0]) * corr * r0;
  *ec = rhot0 * ((vc[1] - vc[0]) * corr) * r0;
  *dec = (dc[1] - dc[0]) * corr * r0;
  r0 = 0.0;
  for (i = 0; i < 2; i++) {
    r0 += f[i] * e[i];
  }

  for (i = 0; i < (int32_T)nx; i++) {
    rhot0 = rt_powd_snf(r[(int32_T)(1.0 + (real_T)i) - 1], 3.0);
    *eel += 0.5 * rhot[(int32_T)(1.0 + (real_T)i) - 1] * vhe[(int32_T)(1.0 +
      (real_T)i) - 1] * 0.048790164169432049 * rhot0;
    *eei += rhot[(int32_T)(1.0 + (real_T)i) - 1] * vio[(int32_T)(1.0 + (real_T)i)
      - 1] * 0.048790164169432049 * rhot0;
    *ec += rhot[(int32_T)(1.0 + (real_T)i) - 1] * vc[(int32_T)(1.0 + (real_T)i)
      - 1] * 0.048790164169432049 * rhot0;

    /* integrale del potenziale * n */
    *dec += dc[(int32_T)(1.0 + (real_T)i) - 1] * 0.048790164169432049 * rhot0;
    *dexc += dxc[(int32_T)(1.0 + (real_T)i) - 1] * 0.048790164169432049 * rhot0;
    *exc += rhot[(int32_T)(1.0 + (real_T)i) - 1] * vxc[(int32_T)(1.0 + (real_T)i)
      - 1] * 0.048790164169432049 * rhot0;
  }

  *etot = (r0 + *dexc) - *eel;
  *exc += *dexc;
  *ec += *dec;
  *ekin = ((*etot - *eei) - *eel) - *exc;
}

static real_T fxc1(real_T x, const char_T type_data[20], const int32_T
                   type_size[2])
{
  real_T b_fxc1;
  real_T xs;
  real_T rs;
  b_fxc1 = 0.0;
  if ((!eml_strcmp(type_data, type_size)) && (!b_eml_strcmp(type_data, type_size)))
  {
	printf("nessuna correlazione del tipo trovata");
  } else {
    /*  x=densita**1/3 in a. u. */
    /*  -x*exfac e’ il pot di scambio vx, il resto e’ vc */
    xs = 12.566370614359172 * x;
    b_fxc1 = -x * (0.98474502184269663 + (0.943656 + 8.8963 * x) / ((1.0 + xs) *
                    (1.0 + xs)));
    if (b_eml_strcmp(type_data, type_size)) {
    } else if (x < 0.0001) {
      b_fxc1 = 0.0;

      /* ??? */
    } else {
      rs = 0.6203504908994 / x;
      xs = sqrt(rs);
      if (rs < 1.0) {
        xs = log(rs);
        xs = ((0.0311 * xs + -0.058366666666666664) + 0.004 * rs * xs / 3.0) +
          -0.0252 * rs / 3.0;
      } else {
        xs = -0.1423 / ((1.0 + 1.0529 * xs) + 0.3334 * rs) * ((1.0 +
          1.2283833333333334 * xs) + 0.44453333333333328 * rs) / ((1.0 + 1.0529 *
          xs) + 0.3334 * rs);
      }

      /*  -x*exfac e’ il pot di scambio vx, emuc=vc e’ il pot di correlazione cep-ald */
      b_fxc1 = -x * 0.98474502184269663 + xs;
    }
  }

  return b_fxc1;
}

static void poisson(const real_T r[500], real_T z, real_T r0, const real_T coef
                    [5], const real_T chi[2500], real_T nx, const real_T eps[5],
                    const real_T f[5], const char_T type_data[20], const int32_T
                    type_size[2], real_T vtot[500], real_T *etot, real_T *exc,
                    real_T *ec, real_T *ekin, real_T *eel, real_T *eei, real_T
                    *dec, real_T *dexc)
{
  real_T q[2500];
  real_T aux[2500];
  int32_T i;
  real_T vh[500];
  real_T vio[500];
  real_T vxc[500];
  real_T dxc[500];
  real_T dc[500];
  real_T vc[500];
  real_T n[500];
  real_T auxx[5];
  int32_T j;
  real_T r13;
  boolean_T bv1[500];
  int32_T tmp_size[1];
  int32_T tmp_data[500];
  for (i = 0; i < 2500; i++) {
    q[i] = 0.0;
    aux[i] = 0.0;
  }

  for (i = 0; i < 500; i++) {
    vh[i] = 0.0;
    vio[i] = 0.0;
    vxc[i] = 0.0;
    dxc[i] = 0.0;
    dc[i] = 0.0;
    vc[i] = 0.0;
    n[i] = 0.0;
  }

  for (i = 0; i < 5; i++) {
    auxx[i] = 0.0;
  }

  for (j = 0; f[j] > 0.0; j++) {
    r13 = coef[j] * coef[j] * rt_powd_snf(r0, 3.0) / 3.0;

    /* serve da r=0 al primo punto della griglia */
    for (i = 0; i < (int32_T)nx; i++) {
      r13 += chi[((int32_T)(1.0 + (real_T)i) + 500 * j) - 1] * chi[((int32_T)
        (1.0 + (real_T)i) + 500 * j) - 1] * 0.048790164169432049 * r[(int32_T)
        (1.0 + (real_T)i) - 1];
      q[((int32_T)(1.0 + (real_T)i) + 500 * j) - 1] = r13;

      /* calcolo carica radiale */
    }
  }

  for (j = 0; f[j] > 0.0; j++) {
    r13 = coef[j] * coef[j] * rt_powd_snf(r0, 3.0) / 2.0;
    for (i = 0; i < (int32_T)nx; i++) {
      r13 += chi[((int32_T)(1.0 + (real_T)i) + 500 * j) - 1] * chi[((int32_T)
        (1.0 + (real_T)i) + 500 * j) - 1] * 0.048790164169432049;
      aux[((int32_T)(1.0 + (real_T)i) + 500 * j) - 1] = r13;

      /* integrale in dr di (chi*chi/r) */
    }
  }

  for (j = 0; f[j] > 0.0; j++) {
    auxx[j] = aux[((int32_T)nx + 500 * j) - 1];
    for (i = 0; i < (int32_T)nx; i++) {
      aux[i + 500 * j] -= auxx[j];

      /* sottraggo la costante di integrazione auxx per ottenere vh=1/r a grandi r */
    }
  }

  for (i = 0; i < (int32_T)nx; i++) {
    j = 0;
    vh[i] = 0.0;
    while (f[j] > 0.0) {
      vh[(int32_T)(1.0 + (real_T)i) - 1] += f[j] * (q[((int32_T)(1.0 + (real_T)i)
        + 500 * j) - 1] / r[(int32_T)(1.0 + (real_T)i) - 1] - aux[((int32_T)(1.0
        + (real_T)i) + 500 * j) - 1]);
      j++;
    }

    vio[(int32_T)(1.0 + (real_T)i) - 1] = -(z / r[(int32_T)(1.0 + (real_T)i) - 1]);

    /* calcolo potenziale di hatree e totale    */
  }

  for (i = 0; i < (int32_T)nx; i++) {
    for (j = 0; f[j] > 0.0; j++) {
      n[i] += f[j] * (chi[((int32_T)(1.0 + (real_T)i) + 500 * j) - 1] * chi
                      [((int32_T)(1.0 + (real_T)i) + 500 * j) - 1] /
                      12.566370614359172) / (r[(int32_T)(1.0 + (real_T)i) - 1] *
        r[(int32_T)(1.0 + (real_T)i) - 1]);
    }
  }

  for (i = 0; i < 500; i++) {
    bv1[i] = rtIsNaN(n[i]);
  }

  eml_li_find(bv1, tmp_data, tmp_size);
  j = tmp_size[0];
  for (i = 0; i < j; i++) {
    n[tmp_data[i] - 1] = 0.0;
  }

  for (i = 0; i < 500; i++) {
    n[i] = n[i] * 4.0 * 3.1415926535897931;
  }

  for (i = 0; i < (int32_T)nx; i++) {
    r13 = 0.0;

    /*  0.43012701 la radice cubica dell’inverso di 4*pi */
    if (n[i] > 1.0E-10) {
      r13 = 0.43012701 * rt_powd_snf(n[(int32_T)(1.0 + (real_T)i) - 1],
        0.33333333);
    }

    vxc[(int32_T)(1.0 + (real_T)i) - 1] = fxc1(r13, type_data, type_size);
    dxc[(int32_T)(1.0 + (real_T)i) - 1] = dxc1(r13, type_data, type_size);

    /*  qui si ricava correlaz sottraendo da scambio+correlaz lo scambio */
    /*  attenzione dc contiene un fattore 4*pi (per questo 3.091 anziche’ 0.246) */
    vc[(int32_T)(1.0 + (real_T)i) - 1] = vxc[(int32_T)(1.0 + (real_T)i) - 1] +
      0.984 * r13;
    dc[(int32_T)(1.0 + (real_T)i) - 1] = dxc[(int32_T)(1.0 + (real_T)i) - 1] -
      3.094 * rt_powd_snf(r13, 4.0);
  }

  for (i = 0; i < 500; i++) {
    vtot[i] = (vio[i] + vh[i]) + vxc[i];
  }

  etotal(r, n, vh, vio, vxc, vc, dxc, dc, eps, f, nx, etot, exc, ec, ekin, eel,
         eei, dec, dexc);
}

static real_T rdivide(real_T x, real_T y)
{
  return x / y;
}

static real_T rt_powd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T d0;
  real_T d1;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    d0 = fabs(u0);
    d1 = fabs(u1);
    if (rtIsInf(u1)) {
      if (d0 == 1.0) {
        y = rtNaN;
      } else if (d0 > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d1 == 0.0) {
      y = 1.0;
    } else if (d1 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

void Hartree(real_T z, const real_T f[5], const char_T type_data[20], const
             int32_T type_size[2], real_T *b_etotal, real_T r[500], real_T q
             [2500], real_T *autov, real_T eps[5], real_T *eei, real_T *eel,
             real_T *exc, real_T *dexc, real_T *ec, real_T *dec, real_T *ekin)
{
  real_T chi[2500];
  real_T vtot[500];
  real_T coef[5];
  real_T deltae[5];
  int32_T i;
  real_T r1;
  real_T nx;
  real_T vold[500];
  boolean_T exitg1;
  boolean_T bv0[500];
  int32_T b_i;
  int32_T tmp_size[1];
  int32_T tmp_data[500];
  int32_T ix;
  real_T epsold[5];
  real_T mtmp;
  boolean_T exitg2;

  /*  Questa funzione esegue il metodo di Hartree per un atomo di Helio con due elettroni. */
  /*  paramentri da passare all funzione: */
  /*  z = Numero atomico; */
  /*  f = arary occupazione; */
  /*  type =  tipo correlazione "cep-ald","wigner" */
  /*  la funzione ritorna: */
  /*  -EnergyStory : l'andamento dell'energia totale calcolata fino alla convergenza */
  /*  -r : il raggio della griglia logaritmica su cui viene eseguito il calcolo */
  /*  -chiFinal : la funzione d'onda calcolata */
  /*   */
  /* serve nelle dichiarazioni successive */
  *b_etotal = 0.0;
  *eei = 0.0;
  *eel = 0.0;
  *exc = 0.0;
  *dexc = 0.0;
  *ec = 0.0;
  *dec = 0.0;
  *ekin = 0.0;
  *autov = 0.0;
  memset(&chi[0], 0, 2500U * sizeof(real_T));
  memset(&vtot[0], 0, 500U * sizeof(real_T));
  for (i = 0; i < 5; i++) {
    eps[i] = 0.0;
    coef[i] = 0.0;
    deltae[i] = 0.0;
  }

  memset(&q[0], 0, 2500U * sizeof(real_T));
  memset(&r[0], 0, 500U * sizeof(real_T));

  r1 = 0.01 / z;

  /* parametri della griglia radiale */
  nx = floor(log(rdivide(z * 75.0, 0.005)) / 0.048790164169432049);
  if (nx > 500.0) {
    /* controllo dimensioni: ho sfondato? */
  } else {
    /* raggio per primo pezzetto integrali */
    Griglia(r1, r, nx);

    /* costruisco griglia radiale r(i) */
    /*  scegli Z* ottimale = Z-5/16 e calcola risultati variazionali */
    /* comincia ciclo autoconsistente partendo da  fz. variazionale (full feedback)                                                           */
    vtot[0] = ThomasFermi(r[0], z);
    for (i = 1; i - 1 < (int32_T)(nx + -1.0); i++) {
      vtot[i] = ThomasFermi(r[i], z);
    }

    for (i = 0; f[i] > 0.0; i++) {
      eps[i] = -((z - 0.3125) * (z - 0.3125)) / (2.0 * (real_T)((i + 1) * (i + 1)));

      /* valore iniziale degli autovalori da cui partire */
      if (eps[i] > vtot[(int32_T)nx - 1]) {
        eps[i] = 2.0 * vtot[(int32_T)nx - 1];
      }
    }

    memcpy(&vold[0], &vtot[0], 500U * sizeof(real_T));
    i = 0;
    exitg1 = FALSE;
    while ((exitg1 == FALSE) && (i < 100)) {
      /* ciclo iterativo, max numero iterazioni=100 */
      for (b_i = 0; b_i < 500; b_i++) {
        bv0[b_i] = rtIsNaN(vtot[b_i]);
      }

      eml_li_find(bv0, tmp_data, tmp_size);
      b_i = tmp_size[0];
      for (ix = 0; ix < b_i; ix++) {
        vtot[tmp_data[ix] - 1] = 0.0;
      }

      for (ix = 0; ix < 500; ix++) {
        vtot[ix] = 0.5 * vtot[ix] + 0.5 * vold[ix];
      }

      for (b_i = 0; b_i < 5; b_i++) {
        epsold[b_i] = eps[b_i];
      }

      /* precisione relativa (criterio accuratezza) per schroed */
      for (b_i = 0; f[b_i] > 0.0; b_i++) {
        Schroed(b_i + 1, &eps[b_i], z, r, vtot, nx, *(real_T (*)[500])&chi[500 *
                b_i]);

        /* Risolvo l'equazione di Hartree per l'orbitale 1s */
      }

      for (b_i = 0; f[b_i] > 0.0; b_i++) {
        coef[b_i] = (chi[500 * b_i] * r[1] / r[0] - chi[1 + 500 * b_i] * r[0] /
                     r[1]) / (r[1] - r[0]);

        /* coef = derivata prima radiale (approx. numerica) di chi in r=0 */
        deltae[b_i] = fabs(eps[b_i] - epsold[b_i]);

        /* differenza di energia da un'iteraz. all'altra */
      }

      memcpy(&vold[0], &vtot[0], 500U * sizeof(real_T));
      poisson(r, z, r1 / 1.02469507659596, coef, chi, nx, eps, f, type_data,
              type_size, vtot, b_etotal, exc, ec, ekin, eel, eei, dec, dexc);
      b_i = 1;
      mtmp = deltae[0];
      if (rtIsNaN(deltae[0])) {
        ix = 2;
        exitg2 = FALSE;
        while ((exitg2 == FALSE) && (ix < 6)) {
          b_i = ix;
          if (!rtIsNaN(deltae[ix - 1])) {
            mtmp = deltae[ix - 1];
            exitg2 = TRUE;
          } else {
            ix++;
          }
        }
      }

      if (b_i < 5) {
        while (b_i + 1 < 6) {
          if (deltae[b_i] > mtmp) {
            mtmp = deltae[b_i];
          }

          b_i++;
        }
      }

      if (mtmp < 1.0E-8) {
        exitg1 = TRUE;
      } else {
        i++;
      }
    }

    /* fine ciclo iterativo */
    /*  se arriva qui vuol dire che deltar non e' minore di 0.00000001 au: avvertire */
    for (b_i = 0; f[b_i] > 0.0; b_i++) {
      *autov += f[b_i] * eps[b_i];
    }

    for (b_i = 0; f[b_i] > 0.0; b_i++) {
      for (i = 0; i < (int32_T)nx; i++) {
        q[i + 500 * b_i] = chi[((int32_T)(1.0 + (real_T)i) + 500 * b_i) - 1] *
          chi[((int32_T)(1.0 + (real_T)i) + 500 * b_i) - 1];
      }
    }
  }
}

void Hartree_initialize(void)
{
  rt_InitInfAndNaN(8U);
}

void Hartree_terminate(void)
{
  /* (no terminate code required) */
}

/* End of code generation (Hartree.c) */
