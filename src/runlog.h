/* runlog.h - shared run-logging + I/O configuration for the spectral RT model.
 *
 * Header-only (static functions). Include after the standard headers.
 * Behaviour is controlled entirely by environment variables, so the same
 * binary serves a single interactive run and a large batch campaign without
 * touching the positional command-line arguments.
 *
 *   LBL_ATMFILE  override the atmospheric profile path
 *                (default: whatever the caller passes to lbl_atmfile()).
 *   LBL_SIMDAYS  simulated duration in days (default 1095 = 3*365).
 *                Short values let the full pipeline be validated quickly.
 *   LBL_OUTFILE  if set, a per-step CSV trajectory is written here.
 *   LBL_META     free-text metadata written as a comment header in the CSV
 *                (e.g. "profile=05 method=montecarlo K=1000 seed=42").
 *
 * The CSV has one row per time step:
 *     step,day,dt_s,Tsfc,Ttoa
 * Tsfc is the lowest-layer temperature, Ttoa the top-layer temperature.
 */

#ifndef RUNLOG_H
#define RUNLOG_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Return the atmosphere file path: LBL_ATMFILE if set, else the default. */
static const char *lbl_atmfile (const char *deflt)
{
  const char *e = getenv ("LBL_ATMFILE");
  return (e && e[0]) ? e : deflt;
}

/* Return the simulated-time bound in seconds (LBL_SIMDAYS days, default 1095). */
static double lbl_tmax_seconds (void)
{
  const char *e = getenv ("LBL_SIMDAYS");
  double days = 3.0 * 365.0;
  if (e && e[0]) {
    double v = strtod (e, NULL);
    if (v > 0.0)
      days = v;
  }
  return days * 86400.0;
}

/* Open the trajectory log if LBL_OUTFILE is set; write metadata + column
 * header. Returns NULL when no logging is requested (original behaviour). */
static FILE *lbl_open_log (void)
{
  const char *path = getenv ("LBL_OUTFILE");
  if (!path || !path[0])
    return NULL;

  FILE *f = fopen (path, "w");
  if (!f) {
    fprintf (stderr, "Warning: could not open LBL_OUTFILE '%s' for writing\n", path);
    return NULL;
  }

  const char *meta = getenv ("LBL_META");
  if (meta && meta[0])
    fprintf (f, "# %s\n", meta);
  fprintf (f, "step,day,dt_s,Tsfc,Ttoa\n");
  fflush (f);
  return f;
}

/* Write one trajectory row. No-op when f is NULL. */
static void lbl_log_step (FILE *f, int step, double t_seconds, double dt_seconds,
                          double Tsfc, double Ttoa)
{
  if (!f)
    return;
  fprintf (f, "%d,%.5f,%.2f,%.5f,%.5f\n",
           step, t_seconds / 86400.0, dt_seconds, Tsfc, Ttoa);
}

#endif /* RUNLOG_H */
