/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#ifndef QUADFILE_H
#define QUADFILE_H

#include <sys/types.h>
#include <stdint.h>

#include "astrometry/qfits_header.h"
#include "astrometry/fitsbin.h"
#include "astrometry/anqfits.h"

typedef struct {
    unsigned int numquads;
    unsigned int numstars;
    int dimquads;
    // upper bound of AB distance of quads in this index
    double index_scale_upper;
    // lower bound
    double index_scale_lower;
    // unique ID of this index
    int indexid;
    // healpix covered by this index
    int healpix;
    // Nside of the healpixelization
    int hpnside;

    fitsbin_t* fb;
    // when reading:
    uint32_t* quadarray;
} quadfile_t;

quadfile_t* quadfile_open(const char* fname);
quadfile_t* quadfile_open_fits(anqfits_t* fits);
quadfile_t* quadfile_open_fits_metadata(anqfits_t* fits);

char* quadfile_get_filename(const quadfile_t* qf);

quadfile_t* quadfile_open_for_writing(const char* quadfname);

quadfile_t* quadfile_open_in_memory(void);

int quadfile_switch_to_reading(quadfile_t* qf);

int quadfile_close(quadfile_t* qf);

// Look at each quad, and ensure that the star ids it contains are all
// less than the number of stars ("numstars").  Returns 0=ok, -1=problem
int quadfile_check(const quadfile_t* qf);

// Copies the star ids of the stars that comprise quad "quadid".
// There will be qf->dimquads such stars.
// (this will be less than starutil.h : DQMAX, for ease of static
// allocation of arrays that will hold quads of stars)
int quadfile_get_stars(const quadfile_t* qf, unsigned int quadid,
                       unsigned int* stars);

// Compatibility advisory wrapper. Returns zero unless validation or delivery
// fails. The subsequent mapped lookup remains authoritative.
int quadfile_prefetch_stars(const quadfile_t* qf,
                            const unsigned int* quadids,
                            int nquads);

// Strict internal preparation contract: positive means every requested row is
// ready, zero means inapplicable, and -1 means complete preparation failed.
int quadfile_prepare_stars(const quadfile_t* qf,
                           const unsigned int* quadids,
                           int nquads);

/*
 * Submit a complete bounded set of Quad rows to the payload loader. The
 * returned ticket must be waited or cancelled before the quadfile is closed.
 * FITSBIN_PAYLOAD_IO_SUBMIT_READY returns without a ticket when the exact
 * live-mapping completion record already covers every requested page. Refusal
 * returns zero and leaves the original mapped lookup authoritative.
 */
int quadfile_prefetch_stars_submit(
    const quadfile_t* qf,
    const unsigned int* quadids,
    int nquads,
    fitsbin_payload_io_ticket_t** ticket);

/*
 * Advise the mapped pages containing the selected Quad rows. This does not
 * read the rows and never uses the compatibility payload descriptor. It
 * returns the number of advised spans, zero when inapplicable, or -1 on
 * invalid input or advice failure.
 */
int quadfile_advise_rows(const quadfile_t* qf,
                         const unsigned int* quadids,
                         int nquads);

int quadfile_write_quad(quadfile_t* qf, unsigned int* stars);

int quadfile_dimquads(const quadfile_t* qf);

int quadfile_nquads(const quadfile_t* qf);

int quadfile_fix_header(quadfile_t* qf);

int quadfile_write_header(quadfile_t* qf);

double quadfile_get_index_scale_upper_arcsec(const quadfile_t* qf);

double quadfile_get_index_scale_lower_arcsec(const quadfile_t* qf);

qfits_header* quadfile_get_header(const quadfile_t* qf);

int quadfile_write_header_to(quadfile_t* qf, FILE* fid);

int quadfile_write_all_quads_to(quadfile_t* qf, FILE* fid);

#endif
