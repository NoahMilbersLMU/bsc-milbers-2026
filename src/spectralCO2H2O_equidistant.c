#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "ascii.h"
#include "schwarzschild.h"
#include "runlog.h"

#define G 9.80665
#define R 287.0
#define CP 1005.0

#define K_BOLTZMANN 1.380649E-23
#define H_PLANCK 6.62607015E-34
#define C_LIGHT 299792458

double planck (double T, double wavelength)
{
  wavelength*=1e-9;  /* convert from nm to m */

  return 2.0*H_PLANCK*C_LIGHT*C_LIGHT/wavelength/wavelength/wavelength/wavelength/wavelength/(exp(H_PLANCK*C_LIGHT/wavelength/K_BOLTZMANN/T)-1.0)/1.0e9;  /* W/(m2 nm sterad) */
}


static int sortlayer (const void *x1, const void *x2);

int main (int argc, char **argv)
{
  //  double Q0=1361.0/4.0*(1.0-0.3);
  double Q0=240;
  double dt=10000; //86400/4.0;
  double dTmax=5.0;   /* maximum temperature increment per time step */
  double dtmax=40000;
  double p0=1000;  /* required for potential temperature */
  int nlev=0, nlyr=0, nwvl=0, ny=0;

  int adaptive=1;    /* adaptive time step */
  int convection=1;  /* convection or no convection */

  /* ==========================================================
   * METHOD-SPECIFIC DECLARATIONS  (Equidistant, rotating)
   * No RNG: this method is fully deterministic.
   * ========================================================== */
  int K=0;             /* requested number of equidistant samples */
  int *idx=NULL;       /* selected wavelength indices */
  int stride=0;        /* index spacing between selected wavelengths */
  int offset=0;        /* current rotation offset (0..stride-1) */
  int nwvl_use=0;      /* number of wavelengths used this timestep */
  /* ========================================================== */

  char atmfilename[FILENAME_MAX]="lbl.arts/fpda.lbl.atm";
  char tauCO2filename[FILENAME_MAX]="./lbl.arts/lbl.co2.asc";
  char tauH2Ofilename[FILENAME_MAX]="./lbl.arts/lbl.h2o.asc";
  char tauCH4filename[FILENAME_MAX]="./lbl.arts/lbl.ch4.asc";
  char tauO3filename [FILENAME_MAX]="./lbl.arts/lbl.o3.asc";
  char tauN2Ofilename[FILENAME_MAX]="./lbl.arts/lbl.n2o.asc";

  int status = 0;

  double t=0;
  double *plev=NULL, *Tlev=NULL;
  double *wvl=NULL, **tauCO2H2O=NULL;
  double **tauCO2_raw=NULL, **tauH2O_raw=NULL;
  double **tauCH4_raw=NULL, **tauO3_raw=NULL, **tauN2O_raw=NULL;
  double *plyr=NULL, *Tlyr=NULL, *theta=NULL, *Blyr=NULL;
  double *tmp=NULL;
  double *hr=NULL;
  double Bg=0;
  double hrmax=0;

  int iv=0, ilyr=0;

  double *edn_lambda=NULL, *eup_lambda=NULL;
  double *edn_lambda_save=NULL, *eup_lambda_save=NULL;
  double *edn=NULL, *eup=NULL;


  char *dummy=NULL;
  int h2oabs=1;
  double co2mult=1.0;


  if (argc<3 || argc>4) {
    fprintf (stderr, "Usage: %s <h2o 0/1> <CO2 multiplier> [K]\n", argv[0]);
    fprintf (stderr, "  K: number of equidistant wavelength samples (0 = full LBL)\n");
    return -1;
  }

  /* read tau from input parameter */
  h2oabs  = strtol (argv[1], &dummy, 0);
  co2mult = strtod (argv[2], &dummy);

  if (argc>=4)
    K = strtol (argv[3], &dummy, 0);

  /* allow the atmosphere path to be overridden via LBL_ATMFILE */
  snprintf (atmfilename, FILENAME_MAX, "%s", lbl_atmfile (atmfilename));


  /* first read atmospheric profile */
  status = read_3c_file (atmfilename, &tmp, &plev, &Tlev, &nlev);
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, atmfilename);
    return status;
  }

  fprintf (stderr, " ... read %d levels from %s\n", nlev, atmfilename);

  /* setting temperature to 288 K */
  for (int ilev=0; ilev<nlev; ilev++)
    Tlev[ilev] = 288;

  nlyr = nlev-1;

  /* temperature and pressure of layers */
  Tlyr  = calloc (nlyr, sizeof(double));
  plyr  = calloc (nlyr, sizeof(double));
  theta = calloc (nlyr, sizeof(double));
  Blyr  = calloc (nlyr, sizeof(double));
  hr    = calloc (nlyr, sizeof(double));

  eup   = calloc (nlev, sizeof(double));
  edn   = calloc (nlev, sizeof(double));

  eup_lambda_save = calloc (nlev, sizeof(double));
  edn_lambda_save = calloc (nlev, sizeof(double));


  for (ilyr=0; ilyr<nlyr; ilyr++)  {
    Tlyr[ilyr] = 288 - (nlyr-1-ilyr)*2; //(Tlev[ilyr] + Tlev[ilyr+1])/2.0;
    plyr[ilyr] = (plev[ilyr] + plev[ilyr+1])/2.0;
  }

  for (int ilev=0; ilev<nlev-1; ilev++)
    fprintf (stderr, "%2d %7.1f %5.1f   %7.1f %5.1f\n",
	     ilev, plev[ilev], Tlev[ilev], plyr[ilev], Tlyr[ilev]);

  /* read CO2 optical thickness profile */
  status = ASCII_file2xy2D (tauCO2filename,
			    &nwvl, &ny,
			    &wvl, &tauCO2_raw);
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, tauCO2filename);
    return status;
  }

  fprintf (stderr, " ... read %d wavelengths and %d layers from %s\n", nwvl, ny, tauCO2filename);

  if (ny!=nlyr) {
    fprintf (stderr, "Error, number of layers in %s (%d) and %s (%d) differ!\n",
	     atmfilename, nlyr, tauCO2filename, ny);
    return -1;
  }

  /* read H2O optical thickness profile */
  if (h2oabs!=0) {
    status = ASCII_file2xy2D (tauH2Ofilename,
			      &nwvl, &ny,
			      &wvl, &tauH2O_raw);
    if (status!=0) {
      fprintf (stderr, "Error %d reading %s\n", status, tauCO2filename);
      return status;
    }

    fprintf (stderr, " ... read %d wavelengths and %d layers from %s\n", nwvl, ny, tauH2Ofilename);

    if (ny!=nlyr) {
      fprintf (stderr, "Error, number of layers in %s (%d) and %s (%d) differ!\n",
	       atmfilename, nlyr, tauH2Ofilename, ny);
      return -1;
    }
  }


  /* read CH4 optical thickness profile */
  status = ASCII_file2xy2D (tauCH4filename,
			    &nwvl, &ny,
			    &wvl, &tauCH4_raw);
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, tauCH4filename);
    return status;
  }

  fprintf (stderr, " ... read %d wavelengths and %d layers from %s\n", nwvl, ny, tauCH4filename);

  if (ny!=nlyr) {
    fprintf (stderr, "Error, number of layers in %s (%d) and %s (%d) differ!\n",
	     atmfilename, nlyr, tauCH4filename, ny);
    return -1;
  }

  /* read N2O optical thickness profile */
  status = ASCII_file2xy2D (tauN2Ofilename,
			    &nwvl, &ny,
			    &wvl, &tauN2O_raw);
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, tauN2Ofilename);
    return status;
  }

  fprintf (stderr, " ... read %d wavelengths and %d layers from %s\n", nwvl, ny, tauN2Ofilename);

  if (ny!=nlyr) {
    fprintf (stderr, "Error, number of layers in %s (%d) and %s (%d) differ!\n",
	     atmfilename, nlyr, tauN2Ofilename, ny);
    return -1;
  }

  /* read O3 optical thickness profile */
  status = ASCII_file2xy2D (tauO3filename,
			    &nwvl, &ny,
			    &wvl, &tauO3_raw);
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, tauO3filename);
    return status;
  }

  fprintf (stderr, " ... read %d wavelengths and %d layers from %s\n", nwvl, ny, tauO3filename);

  if (ny!=nlyr) {
    fprintf (stderr, "Error, number of layers in %s (%d) and %s (%d) differ!\n",
	     atmfilename, nlyr, tauO3filename, ny);
    return -1;
  }


  /* add water vapor and CO2 after applying scaling factors */
  tauCO2H2O=calloc(nwvl, sizeof(double));
  for (iv=0; iv<nwvl; iv++) {
    tauCO2H2O[iv] = calloc(nlyr, sizeof(double));

    for (ilyr=0; ilyr<nlyr; ilyr++) {
      if (h2oabs!=0)
	tauCO2H2O[iv][ilyr] = co2mult * tauCO2_raw[iv][ilyr] + tauH2O_raw[iv][ilyr];
      else
	tauCO2H2O[iv][ilyr] = co2mult * tauCO2_raw[iv][ilyr];

      tauCO2H2O[iv][ilyr] +=  tauCH4_raw[iv][ilyr] +  tauN2O_raw[iv][ilyr] +  tauO3_raw[iv][ilyr];
    }
  }


  /* free memory of tau_updn */
  ASCII_free_double (tauCO2_raw, nwvl);
  if (h2oabs!=0)
    ASCII_free_double (tauH2O_raw, nwvl);

  /* ==========================================================
   * METHOD-SPECIFIC PRE-LOOP SETUP  (Equidistant, rotating)
   * Compute stride = floor(nwvl/K). The index array is rebuilt
   * every time step from the current offset, so allocate it
   * with enough headroom for that rebuild.
   * ========================================================== */
  if (K>0 && K<nwvl) {
    stride = nwvl / K;
    /* allocate full spectrum: the per-step rebuild can write up to
       (nwvl/stride + 3) indices (leading 0 + stride loop + trailing nwvl-1),
       so size to nwvl to be immune to stride arithmetic. */
    idx = calloc(nwvl, sizeof(int));
    fprintf (stderr, " ... Rotating equidistant mode: K=%d requested, stride=%d\n", K, stride);
  } else {
    fprintf (stderr, " ... Full LBL mode (no sampling)\n");
    stride = 0;
    nwvl_use = nwvl;
    idx = calloc(nwvl, sizeof(int));
    for (int i=0; i<nwvl; i++)
      idx[i] = i;
  }
  /* ========================================================== */

  double tmax = lbl_tmax_seconds ();
  FILE *logf = lbl_open_log ();
  int step = 0;

  while (t<tmax) {

    /* ==========================================================
     * METHOD-SPECIFIC TIME-LOOP BODY  (Equidistant, rotating)
     * Rebuild the equidistant index list at the current offset,
     * always append the last grid point so the trapezoidal sum
     * covers the full spectral range, then advance the offset by
     * one for the next step (wrapping inside [0, stride)).
     * Accumulation uses the same trapezoidal rule as in the
     * Shuffle methods over the sorted idx[] array.
     * ========================================================== */
    if (stride > 0) {
      nwvl_use = 0;
      /* always include the first wavelength for full spectral coverage */
      if (offset != 0) {
	idx[nwvl_use] = 0;
	nwvl_use++;
      }
      for (int i=offset; i<nwvl; i+=stride) {
	idx[nwvl_use] = i;
	nwvl_use++;
      }
      /* always include the last wavelength for full spectral coverage */
      if (idx[nwvl_use-1] != nwvl-1) {
	idx[nwvl_use] = nwvl-1;
	nwvl_use++;
      }
      offset = (offset + 1) % stride;
    }
    /* ========================================================== */

    /* reset edn and eup */
    for (int ilev=0; ilev<nlev; ilev++) {
      edn[ilev] = 0;
      eup[ilev] = 0;
    }

    /* spectral integration over the (rebuilt) idx[] subset */
    for (int ik=0; ik<nwvl_use; ik++) {
      iv = idx[ik];

      for (ilyr=0; ilyr<nlyr; ilyr++)
	Blyr[ilyr] = planck (Tlyr[ilyr], wvl[iv]);

      Bg = Blyr[nlyr-1];

      status = schwarzschild (tauCO2H2O[iv], Blyr, nlev, Bg, &edn_lambda, &eup_lambda);
      if (status!=0) {
	fprintf (stderr, "Error %d returned by schwarzschild()\n", status);
	return -1;
      }

      /* trapezoidal sum over adjacent selected wavelengths */
      if (ik>0) {
	int iv_prev = idx[ik-1];
	for (int ilev=0; ilev<nlev; ilev++) {
	  edn[ilev] += (edn_lambda_save[ilev] + edn_lambda[ilev]) / 2.0 * (wvl[iv] - wvl[iv_prev]);
	  eup[ilev] += (eup_lambda_save[ilev] + eup_lambda[ilev]) / 2.0 * (wvl[iv] - wvl[iv_prev]);
	}
      }

      for (int ilev=0; ilev<nlev; ilev++) {
	edn_lambda_save[ilev] = edn_lambda[ilev];
	eup_lambda_save[ilev] = eup_lambda[ilev];
      }

      free(eup_lambda); free(edn_lambda);
    }

    /* calculate heating rates */
    for (int ilev=0; ilev<nlev-1; ilev++)
      hr[ilev] = - G / CP * (edn[ilev+1]-eup[ilev+1] - (edn[ilev]-eup[ilev])) / ((plev[ilev+1]-plev[ilev]) * 100.0);

    /* add surface radiation budget to lowest layer */
    hr[nlev-2] += G / CP * (Q0 + edn[nlev-1] - eup[nlev-1]) / ((plev[nlev-1]-plev[nlev-2]) * 100.0);


    /* determine maximum heating rate */
    if (adaptive) {
      hrmax = 0;
      for (ilyr=0;ilyr<nlyr;ilyr++)
	if (fabs(hr[ilyr]) > hrmax)
	  hrmax = fabs(hr[ilyr]);

      /* time step */
      dt = dTmax/hrmax;
      if (dt>dtmax)
	dt=dtmax;
    }

    t+=dt;

    /* add temperature increment */
    for (int ilev=0; ilev<nlev-1; ilev++)
      Tlyr[ilev] += hr[ilev] * dt;

    for (int ilev=0; ilev<nlev-1; ilev++)
      theta[ilev] = Tlyr[ilev] * pow (p0 / plyr[ilev], 2.0/7);

    if (convection) {
      /* Fix B: energy-conserving dry convective adjustment (enthalpy-conserving
         pairwise mixing; theta_mix = sum(dm*T)/sum(dm*pi), pi=(p/p0)^kappa). */
      int conv_changed = 1, conv_iter = 0;
      while (conv_changed && conv_iter++ < 10000) {
        conv_changed = 0;
        for (int i=0; i<nlev-2; i++) {
          double th_i  = Tlyr[i]   * pow (p0/plyr[i],   2.0/7);
          double th_i1 = Tlyr[i+1] * pow (p0/plyr[i+1], 2.0/7);
          if (th_i < th_i1 - 1.0e-9) {
            double pi_i  = pow (plyr[i]/p0,   2.0/7);
            double pi_i1 = pow (plyr[i+1]/p0, 2.0/7);
            double dm_i  = plev[i+1]-plev[i];
            double dm_i1 = plev[i+2]-plev[i+1];
            double thm = (dm_i*Tlyr[i] + dm_i1*Tlyr[i+1])
                       / (dm_i*pi_i + dm_i1*pi_i1);
            Tlyr[i]   = thm * pi_i;
            Tlyr[i+1] = thm * pi_i1;
            conv_changed = 1;
          }
        }
      }
    }

    fprintf (stdout, "%7.2f %6.1f %6.1f  %6.1f\n", t/86400.0, Tlyr[nlyr-1], Tlyr[0], dt/3600);
    fflush(stdout);

    lbl_log_step (logf, step, t, dt, Tlyr[nlyr-1], Tlyr[0]);
    step++;
  }

  if (logf)
    fclose (logf);

  for (int ilev=0; ilev<nlev-1; ilev++)
    fprintf (stderr, "%2d %7.1f %5.1f %5.1f %5.1f\n", ilev, plyr[ilev], Tlyr[ilev], theta[ilev], hr[ilev] * dt);

  free(idx);

  return 0;
}


static int sortlayer (const void *x1, const void *x2)
{
  double *theta1 = (double *) x1;
  double *theta2 = (double *) x2;

  if (*theta1 < *theta2)
    return 1;

  if (*theta1 > *theta2)
    return -1;

  return 0;
}
