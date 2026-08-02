/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "math.h"

#include "an-bool.h"
#include "anqfits.h"
#include "astrometry/index_shard.h"
#include "engine.h"
#include "engine_internal.h"
#include "engine_private.h"
#include "errors.h"
#include "fileutils.h"
#include "fitsioutils.h"
#include "index_shard_config.h"
#include "ioutils.h"
#include "log.h"
#include "mathutil.h"
#include "onefield.h"
#include "sip-utils.h"
#include "solver.h"

static job_t* job_new() {
    job_t* job = calloc(1, sizeof(job_t));
    if (!job) {
        SYSERROR("Failed to allocate a new job_t.");
        return NULL;
    }
    job->scales = dl_new(8);
    job->depths = il_new(8);
    job->index_shard_workers_override = INDEX_SHARD_WORKERS_UNSET;
    return job;
}

void job_free(job_t* job) {
    if (!job)
        return;
    dl_free(job->scales);
    il_free(job->depths);
    free(job);
}

double engine_job_imagew(job_t* job) {
    return job->bp.solver.field_maxx;
}

double engine_job_imageh(job_t* job) {
    return job->bp.solver.field_maxy;
}

static void parse_sip_coeffs(const qfits_header* hdr, const char* prefix, sip_t* wcs) {
    char key[64];
    int order, i, j;
    sprintf(key, "%sSAO", prefix);
    order = qfits_header_getint(hdr, key, -1);
    if (order >= 2) {
        if (order > 9)
            order = 9;
        wcs->a_order = order;
        wcs->b_order = order;
        for (i=0; i<=order; i++) {
            for (j=0; (i+j)<=order; j++) {
                if (i+j < 1)
                    continue;
                sprintf(key, "%sA%i%i", prefix, i, j);
                wcs->a[i][j] = qfits_header_getdouble(hdr, key, 0.0);
                sprintf(key, "%sB%i%i", prefix, i, j);
                wcs->b[i][j] = qfits_header_getdouble(hdr, key, 0.0);
            }
        }
    }
    sprintf(key, "%sSAPO", prefix);
    order = qfits_header_getint(hdr, key, -1);
    if (order >= 2) {
        if (order > 9)
            order = 9;
        wcs->ap_order = order;
        wcs->bp_order = order;
        for (i=0; i<=order; i++) {
            for (j=0; (i+j)<=order; j++) {
                if (i+j < 1)
                    continue;
                sprintf(key, "%sAP%i%i", prefix, i, j);
                wcs->ap[i][j] = qfits_header_getdouble(hdr, key, 0.0);
                sprintf(key, "%sBP%i%i", prefix, i, j);
                wcs->bp[i][j] = qfits_header_getdouble(hdr, key, 0.0);
            }
        }
    }
}

static anbool parse_job_from_qfits_header(const qfits_header* hdr, job_t* job) {
    onefield_t* bp = &(job->bp);
    solver_t* sp = &(bp->solver);

    double dnil = -LARGE_VAL;
    char *pstr;
    int n;
    anbool run;

    anbool default_tweak = TRUE;
    int default_tweakorder = 2;
    double default_odds_toprint = 1e6;
    double default_odds_tokeep = 1e9;
    double default_odds_tosolve = 1e9;
    double default_odds_totune = 1e6;
    //double default_image_fraction = 1.0;
    char* fn;
    double val;
    char pretty[FITS_LINESZ+1];

    onefield_init(bp);
    // must be in this order because init_parameters handily zeros out sp
    solver_set_default_values(sp);

    // Here we assume that the field's pixel coordinataes go from zero to IMAGEW,H.
    sp->field_maxx = qfits_header_getdouble(hdr, "IMAGEW", dnil);
    sp->field_maxy = qfits_header_getdouble(hdr, "IMAGEH", dnil);
    if ((sp->field_maxx == dnil) || (sp->field_maxy == dnil) ||
        (sp->field_maxx <= 0.0) || (sp->field_maxy <= 0.0)) {
        logerr("Must specify positive \"IMAGEW\" and \"IMAGEH\".\n");
        goto bailout;
    }

    sp->verify_uniformize = qfits_header_getboolean(hdr, "ANVERUNI", sp->verify_uniformize);
    sp->verify_dedup = qfits_header_getboolean(hdr, "ANVERDUP", sp->verify_dedup);

    val = qfits_header_getdouble(hdr, "ANPOSERR", 0.0);
    if (val > 0.0)
        sp->verify_pix = val;
    val = qfits_header_getdouble(hdr, "ANCTOL", 0.0);
    if (val > 0.0)
        sp->codetol = val;
    val = qfits_header_getdouble(hdr, "ANDISTR", 0.0);
    if (val > 0.0)
        sp->distractor_ratio = val;

    onefield_set_solvedout_file  (bp, fn=fits_get_long_string(hdr, "ANSOLVED"));
    free(fn);
    onefield_set_solvedin_file  (bp, fn=fits_get_long_string(hdr, "ANSOLVIN"));
    free(fn);
    onefield_set_match_file   (bp, fn=fits_get_long_string(hdr, "ANMATCH" ));
    free(fn);
    onefield_set_rdls_file    (bp, fn=fits_get_long_string(hdr, "ANRDLS"  ));
    free(fn);
    onefield_set_scamp_file   (bp, fn=fits_get_long_string(hdr, "ANSCAMP" ));
    free(fn);
    onefield_set_wcs_file     (bp, fn=fits_get_long_string(hdr, "ANWCS"   ));
    free(fn);
    onefield_set_corr_file    (bp, fn=fits_get_long_string(hdr, "ANCORR"  ));
    free(fn);
    onefield_set_cancel_file  (bp, fn=fits_get_long_string(hdr, "ANCANCEL"));
    free(fn);

    onefield_set_xcol(bp, fn=fits_get_dupstring(hdr, "ANXCOL"));
    free(fn);
    onefield_set_ycol(bp, fn=fits_get_dupstring(hdr, "ANYCOL"));
    free(fn);

    bp->timelimit = qfits_header_getdouble(hdr, "ANTLIM", 0.0);
    bp->cpulimit = qfits_header_getdouble(hdr, "ANCLIM", 0.0);
    if (qfits_header_getstr(hdr, "ANSHWRK")) {
        int requested_workers =
            qfits_header_getint(hdr, "ANSHWRK", INT_MIN);

        if (requested_workers == INT_MIN) {
            logerr("Invalid ANSHWRK worker value in augmented job header.\n");
            goto bailout;
        }

        job->index_shard_workers_override = requested_workers;
        job->index_shard_workers_override_set = TRUE;
    }
    bp->logratio_tosolve = log(qfits_header_getdouble(hdr, "ANODDSSL", default_odds_tosolve));
    logverb("Set odds ratio to solve to %g (log = %g)\n", exp(bp->logratio_tosolve), bp->logratio_tosolve);


    sp->logratio_toprint = log(qfits_header_getdouble(hdr, "ANODDSPR", default_odds_toprint));
    sp->logratio_tokeep = log(qfits_header_getdouble(hdr, "ANODDSKP", default_odds_tokeep));
    sp->logratio_totune = log(qfits_header_getdouble(hdr, "ANODDSTU", default_odds_totune));
    sp->logratio_bail_threshold = log(qfits_header_getdouble(hdr, "ANODDSBL", DEFAULT_BAIL_THRESHOLD));
    val = qfits_header_getdouble(hdr, "ANODDSST", 0.0);
    if (val > 0.0)
        sp->logratio_stoplooking = log(val);
    bp->best_hit_only = TRUE;

    // gotta keep it to solve it!
    sp->logratio_tokeep = MIN(sp->logratio_tokeep, bp->logratio_tosolve);
    // gotta print it to keep it (so what if that doesn't make sense)!
    sp->logratio_toprint = MIN(sp->logratio_toprint, sp->logratio_tokeep);

    // job->image_fraction = qfits_header_getdouble(hdr, "ANIMFRAC", job->image_fraction);
    job->include_default_scales = qfits_header_getboolean(hdr, "ANAPPDEF", 0);

    sp->parity = PARITY_BOTH;
    pstr = qfits_pretty_string_r(qfits_header_getstr(hdr, "ANPARITY"), pretty);
    if (pstr && streq(pstr, "NEG"))
        sp->parity = PARITY_FLIP;
    else if (pstr && streq(pstr, "POS"))
        sp->parity = PARITY_NORMAL;

    sp->set_crpix_center = qfits_header_getboolean(hdr, "ANCRPIXC", FALSE);
    sp->crpix[0] = qfits_header_getdouble(hdr, "ANCRPIX1", sp->crpix[0]);
    sp->crpix[1] = qfits_header_getdouble(hdr, "ANCRPIX2", sp->crpix[1]);
    sp->set_crpix = (sp->set_crpix_center ||
                     // were the values set?
                     qfits_header_getstr(hdr, "ANCRPIX1") ||
                     qfits_header_getstr(hdr, "ANCRPIX2"));

    if (qfits_header_getboolean(hdr, "ANTWEAK", default_tweak)) {
        int order = qfits_header_getint(hdr, "ANTWEAKO", default_tweakorder);
        //bp->do_tweak = TRUE;
        sp->do_tweak = TRUE;
        sp->tweak_aborder = order;
        sp->tweak_abporder = order;
    }

    if (!sp->do_tweak) {
        // No tweak: set tweak order to linear, because the tweak alg
        // can still be invoked via tune-up.
        sp->tweak_aborder = sp->tweak_abporder = 1;
    }

    val = qfits_header_getdouble(hdr, "ANQSFMIN", 0.0);
    if (val > 0.0)
        bp->quad_size_fraction_lo = val;
    val = qfits_header_getdouble(hdr, "ANQSFMAX", 0.0);
    if (val > 0.0)
        bp->quad_size_fraction_hi = val;

    job->ra_center = qfits_header_getdouble(hdr, "ANERA", dnil);
    job->dec_center = qfits_header_getdouble(hdr, "ANEDEC", dnil);
    job->search_radius = qfits_header_getdouble(hdr, "ANERAD", dnil);
    job->use_radec_center = ((job->ra_center     != dnil) &&
                             (job->dec_center    != dnil) &&
                             (job->search_radius != dnil));

    // tag-along columns
    bp->rdls_tagalong_all = qfits_header_getboolean(hdr, "ANTAGALL", FALSE);
    if (!bp->rdls_tagalong_all) {
        n = 1;
        while (1) {
            char key[64];
            char* val;
            sprintf(key, "ANTAG%i", n);
            val = fits_get_dupstring(hdr, key);
            if (!val)
                break;
            if (!bp->rdls_tagalong)
                bp->rdls_tagalong = sl_new(16);
            sl_append_nocopy(bp->rdls_tagalong, val);
            n++;
        }
    }

    // sort RDLS column
    bp->sort_rdls = fits_get_dupstring(hdr, "ANRDSORT");

    n = 1;
    while (1) {
        char key[64];
        double lo, hi;
        sprintf(key, "ANAPPL%i", n);
        lo = qfits_header_getdouble(hdr, key, 0.);
        sprintf(key, "ANAPPU%i", n);
        hi = qfits_header_getdouble(hdr, key, 0.);
        if ((hi == 0.) && (lo == 0.))
            break;
        if ((lo != 0.) && (hi != 0.)) {
            if ((lo < 0) || (lo > hi)) {
                logerr("Scale range %g to %g is invalid: min must be >= 0, max must be >= min.\n", lo, hi);
                goto bailout;
            }
        }
        dl_append(job->scales, lo);
        dl_append(job->scales, hi);
        n++;
    }

    n = 1;
    while (1) {
        char key[64];
        int dlo, dhi;
        sprintf(key, "ANDPL%i", n);
        dlo = qfits_header_getint(hdr, key, 0);
        sprintf(key, "ANDPU%i", n);
        dhi = qfits_header_getint(hdr, key, 0);
        if (dlo == 0 && dhi == 0)
            break;
        if ((dlo < 1) || (dlo > dhi)) {
            logerr("Depth range %i to %i is invalid: min must be >= 1, max must be >= min.\n", dlo, dhi);
            goto bailout;
        }
        il_append(job->depths, dlo);
        il_append(job->depths, dhi);
        n++;
    }

    n = 1;
    while (1) {
        char lokey[64];
        char hikey[64];
        int lo, hi;
        sprintf(lokey, "ANFDL%i", n);
        lo = qfits_header_getint(hdr, lokey, -1);
        if (lo == -1)
            break;
        sprintf(hikey, "ANFDU%i", n);
        hi = qfits_header_getint(hdr, hikey, -1);
        if (hi == -1)
            break;
        if ((lo <= 0) || (lo > hi)) {
            char pretty1[FITS_LINESZ+1];
            char pretty2[FITS_LINESZ+1];
            logerr("Field range %i to %i is invalid: min must be >= 1, max must be >= min.\n", lo, hi);
            qfits_pretty_string_r(qfits_header_getstr(hdr, lokey), pretty1);
            qfits_pretty_string_r(qfits_header_getstr(hdr, hikey), pretty2);
            logmsg("  (FITS headers: \"%s = %s\", \"%s = %s\")\n",
                   lokey, pretty1, hikey, pretty2);
            goto bailout;
        }

        onefield_add_field_range(bp, lo, hi);
        n++;
    }

    n = 1;
    while (1) {
        char key[64];
        int fld;
        sprintf(key, "ANFD%i", n);
        fld = qfits_header_getint(hdr, key, -1);
        if (fld == -1)
            break;
        if (fld <= 0) {
            qfits_pretty_string_r(qfits_header_getstr(hdr, key), pretty);
            logerr("Field %i is invalid: must be >= 1.  (FITS header: \"%s = %s\")\n", fld, key, pretty);
            goto bailout;
        }

        onefield_add_field(bp, fld);
        n++;
    }

    n = 1;
    while (1) {
        char key[64];
        sip_t wcs;
        char* keys[] = { "ANW%iPIX1", "ANW%iPIX2", "ANW%iVAL1", "ANW%iVAL2",
                         "ANW%iCD11", "ANW%iCD12", "ANW%iCD21", "ANW%iCD22" };
        double* vals[] = { &(wcs.wcstan. crval[0]), &(wcs.wcstan.crval[1]),
                           &(wcs.wcstan.crpix[0]), &(wcs.wcstan.crpix[1]),
                           &(wcs.wcstan.cd[0][0]), &(wcs.wcstan.cd[0][1]),
                           &(wcs.wcstan.cd[1][0]), &(wcs.wcstan.cd[1][1]) };
        int j;
        int bail = 0;
        memset(&wcs, 0, sizeof(wcs));
        for (j = 0; j < 8; j++) {
            sprintf(key, keys[j], n);
            *(vals[j]) = qfits_header_getdouble(hdr, key, dnil);
            if (*(vals[j]) == dnil) {
                bail = 1;
                break;
            }
        }
        if (bail)
            break;

        // SIP terms
        sprintf(key, "ANW%i", n);
        parse_sip_coeffs(hdr, key, &wcs);

        sip_ensure_inverse_polynomials(&wcs);

        onefield_add_verify_wcs(bp, &wcs);
        n++;
    }

    // Distortion to apply before matching...
    do {
        sip_t dsip;
        double p0, p1;
        memset(&dsip, 0, sizeof(sip_t));
        p0 = qfits_header_getdouble(hdr, "ANDPIX0", dnil);
        if (p0 == dnil)
            break;
        p1 = qfits_header_getdouble(hdr, "ANDPIX1", dnil);
        if (p1 == dnil)
            break;
        dsip.wcstan.crpix[0] = p0;
        dsip.wcstan.crpix[1] = p1;
        parse_sip_coeffs(hdr, "AND", &dsip);
        if ((dsip.a_order > 1 && dsip.b_order > 1) ||
            (dsip.ap_order > 1 && dsip.bp_order > 1)) {
            sp->predistort = malloc(sizeof(sip_t));
            memcpy(sp->predistort, &dsip, sizeof(sip_t));
        }
    } while (0);

    sp->pixel_xscale = qfits_header_getdouble(hdr, "ANPXSCAL", 0.);

    run = qfits_header_getboolean(hdr, "ANRUN", FALSE);

    // Default: solve first field.
    if (run && !il_size(bp->fieldlist)) {
        onefield_add_field(bp, 1);
    }

    return TRUE;

 bailout:
    return FALSE;
}

job_t* engine_read_job_file(engine_t* engine, const char* jobfn) {
    qfits_header* hdr;
    job_t* job;
    onefield_t* bp;

    // Read primary header.
    hdr = anqfits_get_header2(jobfn, 0);
    if (!hdr) {
        ERROR("Failed to parse FITS header from file \"%s\"", jobfn);
        return NULL;
    }
    job = job_new();
    if (!parse_job_from_qfits_header(hdr, job)) {
        job_free(job);
        qfits_header_destroy(hdr);
        return NULL;
    }
    qfits_header_destroy(hdr);

    bp = &(job->bp);

    if (engine_resolve_index_shard_workers(engine, job)) {
        solver_cleanup(&bp->solver);
        onefield_cleanup(bp);
        job_free(job);
        return NULL;
    }

    onefield_set_field_file(bp, jobfn);

    // If the job has no scale estimate, search everything provided
    // by the engine
    if (!dl_size(job->scales) || job->include_default_scales) {
        double arcsecperpix;
        arcsecperpix = deg2arcsec(engine->minwidth) /
            engine_job_imagew(job);
        dl_append(job->scales, arcsecperpix);
        arcsecperpix = deg2arcsec(engine->maxwidth) /
            engine_job_imagew(job);
        dl_append(job->scales, arcsecperpix);
    }

    // SECTION INDEX-SHARD: limit-precedence
    /*
     * Wall time is the primary user-facing deadline.  A job may reduce the
     * backend wall limit but may not increase it.  CPU time remains aggregate
     * process user+system time; a job value overrides the optional config
     * default because the site-level hard ceiling is now the wall limit.
     */
    {
      engine_limit_policy_t limits;
      double job_walllimit = bp->timelimit;
      double cfg_walllimit = engine->walllimit;
      double job_cpulimit = bp->cpulimit;
      double cfg_cpulimit = engine->cpulimit;

      engine_limit_policy_resolve(
          job_walllimit, cfg_walllimit,
          job_cpulimit, cfg_cpulimit, &limits);

      bp->timelimit = limits.wall_seconds;
      bp->total_timelimit = limits.wall_seconds;
      if (index_shard_pthread_enabled(bp)) {
        bp->total_cpulimit = limits.cpu_seconds;
        bp->cpulimit = 0.0;
      } else {
        bp->cpulimit = limits.cpu_seconds;
        bp->total_cpulimit = limits.cpu_seconds;
      }

      if (limits.wall_job_clamped) {
        logmsg("Requested wall limit %g s reduced to backend limit %g s.\n",
               job_walllimit, limits.wall_seconds);
      }
      if (limits.wall_seconds > 0.0 && limits.cpu_seconds > 0.0) {
        logmsg("Limits: wall=%g s elapsed engine time; "
               "CPU=%g aggregate process seconds; workers=%i.\n",
               limits.wall_seconds, limits.cpu_seconds,
               bp->index_shard_workers);
      } else if (limits.wall_seconds > 0.0) {
        logmsg("Limits: wall=%g s elapsed engine time; "
               "CPU=unlimited aggregate process time; workers=%i.\n",
               limits.wall_seconds, bp->index_shard_workers);
      } else if (limits.cpu_seconds > 0.0) {
        logmsg("Limits: wall=unlimited elapsed engine time; "
               "CPU=%g aggregate process seconds; workers=%i.\n",
               limits.cpu_seconds, bp->index_shard_workers);
      } else {
        logmsg("Limits: wall=unlimited elapsed engine time; "
               "CPU=unlimited aggregate process time; workers=%i.\n",
               bp->index_shard_workers);
      }
      logverb("[limits] wall_effective=%g wall_job=%g wall_config=%g "
              "wall_source=%s wall_clamped=%i "
              "cpu_effective=%g cpu_job=%g cpu_config=%g cpu_source=%s\n",
              limits.wall_seconds, job_walllimit, cfg_walllimit,
              limits.wall_from_job ? "job" :
                  (cfg_walllimit > 0.0 ? "config" : "none"),
              limits.wall_job_clamped ? 1 : 0,
              limits.cpu_seconds, job_cpulimit, cfg_cpulimit,
              limits.cpu_from_job ? "job" :
                  (cfg_cpulimit > 0.0 ? "config" : "none"));
    }

    logverb("[index-shard] engine limits after setup: "
            "cpulimit=%f total_cpulimit=%f timelimit=%g total_timelimit=%g\n",
            bp->cpulimit, bp->total_cpulimit, bp->timelimit,
            bp->total_timelimit);

    // If the job didn't specify depths, set defaults.
    if (il_size(job->depths) == 0) {
        if (il_size(engine->default_depths) != 0) {
            il_append_list(job->depths, engine->default_depths);
        } else {
            /*
             * An empty site default means the original unbounded depth
             * interval. Keep this scientific search space independent of
             * worker count and of the legacy "inparallel" token.
             */
            il_append(job->depths, 0);
            il_append(job->depths, 0);
        }
    }

    if (engine->cancelfn)
        onefield_set_cancel_file(bp, engine->cancelfn);
    if (engine->solvedfn)
        onefield_set_solved_file(bp, engine->solvedfn);

    return job;
}

void job_set_cancel_file(job_t* job, const char* fn) {
    onefield_set_cancel_file(&(job->bp), fn);
}

void job_set_solved_file(job_t* job, const char* fn) {
    onefield_set_solved_file(&(job->bp), fn);
}

// Modify all filenames to be relative to "dir".
int job_set_base_dir(job_t* job, const char* dir) {
    return job_set_output_base_dir(job, dir) ||
        job_set_input_base_dir(job, dir);
}

int job_set_input_base_dir(job_t* job, const char* dir) {
    char* path;
    onefield_t* bp = &(job->bp);
    logverb("Changing input file base dir to %s\n", dir);
    if (bp->fieldfname) {
        path = resolve_path(bp->fieldfname, dir);
        logverb("Changing %s to %s\n", bp->fieldfname, path);
        onefield_set_field_file(bp, path);
    }
    return 0;
}

int job_set_output_base_dir(job_t* job, const char* dir) {
    char* path;
    onefield_t* bp = &(job->bp);
    logverb("Changing output file base dir to %s\n", dir);
    if (bp->cancelfname) {
        path = resolve_path(bp->cancelfname, dir);
        logverb("Cancel file was %s, changing to %s.\n", bp->cancelfname, path);
        onefield_set_cancel_file(bp, path);
    }
    if (bp->solved_in) {
        path = resolve_path(bp->solved_in, dir);
        logverb("Changing %s to %s\n", bp->solved_in, path);
        onefield_set_solvedin_file(bp, path);
    }
    if (bp->solved_out) {
        path = resolve_path(bp->solved_out, dir);
        logverb("Changing %s to %s\n", bp->solved_out, path);
        onefield_set_solvedout_file(bp, path);
    }
    if (bp->matchfname) {
        path = resolve_path(bp->matchfname, dir);
        logverb("Changing %s to %s\n", bp->matchfname, path);
        onefield_set_match_file(bp, path);
    }
    if (bp->indexrdlsfname) {
        path = resolve_path(bp->indexrdlsfname, dir);
        logverb("Changing %s to %s\n", bp->indexrdlsfname, path);
        onefield_set_rdls_file(bp, path);
    }
    if (bp->scamp_fname) {
        path = resolve_path(bp->scamp_fname, dir);
        logverb("Changing %s to %s\n", bp->scamp_fname, path);
        onefield_set_scamp_file(bp, path);
    }
    if (bp->corr_fname) {
        path = resolve_path(bp->corr_fname, dir);
        logverb("Changing %s to %s\n", bp->corr_fname, path);
        onefield_set_corr_file(bp, path);
    }
    if (bp->wcs_template) {
        path = resolve_path(bp->wcs_template, dir);
        logverb("Changing %s to %s\n", bp->wcs_template, path);
        onefield_set_wcs_file(bp, path);
    }
    return 0;
}
