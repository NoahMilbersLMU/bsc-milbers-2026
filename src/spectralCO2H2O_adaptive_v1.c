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
   * METHOD-SPECIFIC DECLARATIONS  (Adaptive-Memory Sampling v1)
   *
   * Each wavelength carries a pick weight wgt[iv] in (0,1].
   * When a wavelength is drawn its weight is knocked down to
   * `suppress`; every time step all weights recover linearly
   * back towards 1 by `recover` per step. Recovery is scaled
   * so a picked wavelength is back at full weight after about
   * rectime * (nwvl/K) time steps, i.e. roughly the time the
   * sampler needs to visit the rest of the spectrum once.
   *
   * v1: probability-corrected (Horvitz-Thompson) estimator -
   * every sample is weighted by 1/(pick probability), so each
   * individual time step is an unbiased flux estimate.
   * (v2 uses the plain dwvl*nwvl/K weight instead.)
   * ========================================================== */
  int K=0;               /* number of draws per time step (0 = full LBL) */
  unsigned int seed=0;
  double *dwvl=NULL;     /* per-grid-point bin widths */
  double *wgt=NULL;      /* adaptive pick weight (the memory) */
  double *cdf=NULL;      /* cumulative sum of wgt*dwvl, rebuilt each step */
  int *drawn=NULL;       /* indices drawn in the current time step */
  double suppress=0.05;  /* weight assigned right after being picked */
  double rectime=1.0;    /* recovery time in units of nwvl/K time steps */
  double recover=0.0;    /* additive weight recovery per time step */
  double W_total=0.0;    /* current normalisation sum(wgt) */
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


  if (argc<3 || argc>7) {
    fprintf (stderr, "Usage: %s <h2o 0/1> <CO2 multiplier> [K] [seed] [suppress] [rectime]\n", argv[0]);
    fprintf (stderr, "  K:        number of adaptive draws per time step (0 = full LBL)\n");
    fprintf (stderr, "  suppress: pick weight right after being drawn, 0<s<=1 (default 0.05)\n");
    fprintf (stderr, "  rectime:  weight recovery time in units of nwvl/K steps (default 1.0)\n");
    return -1;
  }

  /* read tau from input parameter */
  h2oabs  = strtol (argv[1], &dummy, 0);
  co2mult = strtod (argv[2], &dummy);

  if (argc>=4)
    K = strtol (argv[3], &dummy, 0);
  if (argc>=5)
    seed = (unsigned int) strtol (argv[4], &dummy, 0);
  else
    seed = (unsigned int) time(NULL);
  if (argc>=6)
    suppress = strtod (argv[5], &dummy);
  if (argc>=7)
    rectime = strtod (argv[6], &dummy);

  if (suppress<=0.0 || suppress>1.0) {
    fprintf (stderr, "Error, suppress must be in (0,1], got %g\n", suppress);
    return -1;
  }
  if (rectime<=0.0) {
    fprintf (stderr, "Error, rectime must be positive, got %g\n", rectime);
    return -1;
  }

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
   * METHOD-SPECIFIC PRE-LOOP SETUP  (Adaptive-Memory Sampling)
   * All weights start at 1 (first draw is plain uniform MC over
   * the bin widths). The per-step recovery increment is chosen
   * so that suppress -> 1 takes rectime * nwvl/K time steps.
   * ========================================================== */
  dwvl = calloc (nwvl, sizeof(double));
  dwvl[0] = (wvl[1] - wvl[0]) / 2.0;
  for (iv=1; iv<nwvl-1; iv++)
    dwvl[iv] = (wvl[iv+1] - wvl[iv-1]) / 2.0;
  dwvl[nwvl-1] = (wvl[nwvl-1] - wvl[nwvl-2]) / 2.0;

  if (K > nwvl) {
    fprintf (stderr, " ... K capped to nwvl=%d\n", nwvl);
    K = nwvl;
  }

  if (K>0) {
    srand(seed);
    wgt   = calloc (nwvl, sizeof(double));
    cdf   = calloc (nwvl, sizeof(double));
    drawn = calloc (K,    sizeof(int));

    for (iv=0; iv<nwvl; iv++)
      wgt[iv] = 1.0;

    recover = (1.0 - suppress) * (double)K / ((double)nwvl * rectime);

    fprintf (stderr, " ... Adaptive-memory sampling v1 (probability-corrected): K=%d, seed=%u, suppress=%.3f, recovery over %.0f steps\n",
	     K, seed, suppress, (1.0 - suppress) / recover);
  } else {
    fprintf (stderr, " ... Full LBL mode (no sampling)\n");
  }
  /* ========================================================== */

  double tmax = lbl_tmax_seconds ();
  FILE *logf = lbl_open_log ();
  int step = 0;

  while (t<tmax) {

    /* reset edn and eup */
    for (int ilev=0; ilev<nlev; ilev++) {
      edn[ilev] = 0;
      eup[ilev] = 0;
    }

    /* ==========================================================
     * METHOD-SPECIFIC TIME-LOOP BODY  (Adaptive-Memory Sampling)
     * Rebuild the CDF over wgt (weights changed last step), draw
     * K grid indices via inverse CDF (binary search), and
     * reweight each sample by dwvl[iv] / (K * p[iv]) with
     * p[iv] = wgt[iv]/W_total, i.e. dwvl[iv]*W_total/(K*wgt[iv]),
     * keeping the estimator unbiased. With all weights equal this
     * reduces exactly to the uniform-MC estimator
     * E_lambda * dwvl[iv] * nwvl/K. Afterwards all weights
     * recover, the drawn ones are suppressed.
     * ========================================================== */
    if (K>0) {
      W_total = 0.0;
      for (iv=0; iv<nwvl; iv++) {
	W_total += wgt[iv];
	cdf[iv] = W_total;
      }

      for (int isamp=0; isamp<K; isamp++) {
	/* draw u ~ U(0, W_total), invert CDF via binary search */
	double u = (double)rand() / ((double)RAND_MAX + 1.0) * W_total;
	int lo=0, hi=nwvl-1;
	while (lo < hi) {
	  int mid = lo + (hi - lo) / 2;
	  if (cdf[mid] <= u)
	    lo = mid + 1;
	  else
	    hi = mid;
	}
	iv = lo;
	drawn[isamp] = iv;

	for (ilyr=0; ilyr<nlyr; ilyr++)
	  Blyr[ilyr] = planck (Tlyr[ilyr], wvl[iv]);

	Bg = Blyr[nlyr-1];

	status = schwarzschild (tauCO2H2O[iv], Blyr, nlev, Bg, &edn_lambda, &eup_lambda);
	if (status!=0) {
	  fprintf (stderr, "Error %d returned by schwarzschild()\n", status);
	  return -1;
	}

	double sample_weight = dwvl[iv] * W_total / ((double)K * wgt[iv]);
	for (int ilev=0; ilev<nlev; ilev++) {
	  edn[ilev] += edn_lambda[ilev] * sample_weight;
	  eup[ilev] += eup_lambda[ilev] * sample_weight;
	}

	free(eup_lambda); free(edn_lambda);
      }

      /* memory update: everyone recovers a bit, drawn ones are suppressed */
      for (iv=0; iv<nwvl; iv++) {
	wgt[iv] += recover;
	if (wgt[iv] > 1.0)
	  wgt[iv] = 1.0;
      }
      for (int isamp=0; isamp<K; isamp++)
	wgt[drawn[isamp]] = suppress;

    } else {
      /* Full LBL: trapezoidal integration over all wavelengths */
      for (iv=1; iv<nwvl; iv++) {
	for (ilyr=0; ilyr<nlyr; ilyr++)
	  Blyr[ilyr] = planck (Tlyr[ilyr], wvl[iv]);

	Bg = Blyr[nlyr-1];

	status = schwarzschild (tauCO2H2O[iv], Blyr, nlev, Bg, &edn_lambda, &eup_lambda);
	if (status!=0) {
	  fprintf (stderr, "Error %d returned by schwarzschild()\n", status);
	  return -1;
	}

	if (iv>0) {
	  for (int ilev=0; ilev<nlev; ilev++) {
	    edn[ilev] += (edn_lambda_save[ilev] + edn_lambda[ilev]) / 2.0 * (wvl[iv] - wvl[iv-1]);
	    eup[ilev] += (eup_lambda_save[ilev] + eup_lambda[ilev]) / 2.0 * (wvl[iv] - wvl[iv-1]);
	  }
	}

	for (int ilev=0; ilev<nlev; ilev++) {
	  edn_lambda_save[ilev] = edn_lambda[ilev];
	  eup_lambda_save[ilev] = eup_lambda[ilev];
	}

	free(eup_lambda); free(edn_lambda);
      }
    }
    /* ========================================================== */

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

  free(dwvl);
  free(wgt);
  free(cdf);
  free(drawn);

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
