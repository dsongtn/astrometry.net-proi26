/*
 This file is part of the Astrometry.net suite.
 Licensed under a 3-clause BSD style license - see LICENSE.
 */

#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "anqfits.h"
#include "qfits_card.h"
#include "qfits_std.h"

static void set_card(
    unsigned char* block,
    size_t card_index,
    const char* text) {
    unsigned char* card = block + card_index * FITS_LINESZ;
    size_t length = strlen(text);

    if (length > FITS_LINESZ) {
        length = FITS_LINESZ;
    }
    memset(card, ' ', FITS_LINESZ);
    memcpy(card, text, length);
}

static int write_all(int fd, const unsigned char* data, size_t size) {
    while (size) {
        ssize_t written = write(fd, data, size);

        if (written < 0) {
            if (errno == EINTR) {
                continue;
            }
            return -1;
        }
        if (!written) {
            errno = EIO;
            return -1;
        }
        data += (size_t)written;
        size -= (size_t)written;
    }
    return 0;
}

static int open_bytes(
    unsigned char* bytes,
    size_t byte_count,
    anqfits_t** result) {
    char path[] = "/tmp/test-anqfits-header.XXXXXX";
    int fd;

    if (!bytes || !byte_count || !result) {
        return -1;
    }
    *result = NULL;
    fd = mkstemp(path);
    if (fd < 0) {
        return -1;
    }
    if (write_all(fd, bytes, byte_count)) {
        int saved_errno = errno;

        close(fd);
        unlink(path);
        errno = saved_errno;
        return -1;
    }
    if (close(fd)) {
        int saved_errno = errno;

        unlink(path);
        errno = saved_errno;
        return -1;
    }
    *result = anqfits_open(path);
    unlink(path);
    return 0;
}

static void set_empty_image_extension(unsigned char* block) {
    set_card(block, 0U, "XTENSION= 'IMAGE   '");
    set_card(block, 1U, "BITPIX  =                    8");
    set_card(block, 2U, "NAXIS   =                    0");
    set_card(block, 3U, "PCOUNT  =                    0");
    set_card(block, 4U, "GCOUNT  =                    1");
    set_card(block, 5U, "END");
}

static int test_header_card_block_boundary(void) {
    unsigned char bytes[2U * FITS_BLOCK_SIZE];
    anqfits_t* fits = NULL;
    size_t card;
    int failed = 0;

    memset(bytes, ' ', sizeof(bytes));
    set_card(bytes, 0U, "SIMPLE  =                    T");
    set_card(bytes, 1U, "BITPIX  =                    8");
    set_card(bytes, 2U, "NAXIS   =                    0");
    set_card(bytes, 3U, "EXTEND  =                    F");
    for (card = 4U; card + 1U < FITS_NCARDS; card++) {
        set_card(bytes, card, "COMMENT bounded header-card regression");
    }
    set_card(bytes, FITS_NCARDS - 1U,
             "BOUNDARY=                    1");
    set_card(bytes + FITS_BLOCK_SIZE, 0U, "END");

    if (open_bytes(bytes, sizeof(bytes), &fits) || !fits) {
        return 1;
    }
    failed += anqfits_n_ext(fits) != 1;
    failed += anqfits_header_size(fits, 0) !=
        (off_t)(2U * FITS_BLOCK_SIZE);
    failed += anqfits_data_size(fits, 0) != 0;
    anqfits_close(fits);
    return failed;
}

static int test_malformed_keys(void) {
    char no_equals[FITS_LINESZ];
    char leading_equals[FITS_LINESZ];
    char key[FITS_LINESZ + 1U];
    int failed = 0;

    memset(no_equals, 'X', sizeof(no_equals));
    memset(leading_equals, ' ', sizeof(leading_equals));
    leading_equals[0] = '=';
    failed += qfits_getkey_r(no_equals, key) != NULL;
    failed += qfits_getkey_r(leading_equals, key) != NULL;
    return failed;
}

static int test_zero_length_first_axis(void) {
    unsigned char bytes[2U * FITS_BLOCK_SIZE];
    unsigned char* extension = bytes + FITS_BLOCK_SIZE;
    anqfits_t* fits = NULL;
    int failed = 0;

    memset(bytes, ' ', sizeof(bytes));
    set_card(bytes, 0U, "SIMPLE  =                    T");
    set_card(bytes, 1U, "BITPIX  =                    8");
    set_card(bytes, 2U, "NAXIS   =                    2");
    set_card(bytes, 3U, "NAXIS1  =                    0");
    set_card(bytes, 4U, "NAXIS2  =                  100");
    set_card(bytes, 5U, "EXTEND  =                    T");
    set_card(bytes, 6U, "END");
    set_empty_image_extension(extension);

    if (open_bytes(bytes, sizeof(bytes), &fits) || !fits) {
        return 1;
    }
    if (anqfits_n_ext(fits) != 2) {
        failed++;
    } else {
        failed += fits->exts[0].data_size != 0;
        failed += fits->exts[1].hdr_start != 1;
    }
    anqfits_close(fits);
    return failed;
}

static int test_random_groups_data_extent(void) {
    unsigned char bytes[3U * FITS_BLOCK_SIZE];
    unsigned char* extension = bytes + 2U * FITS_BLOCK_SIZE;
    anqfits_t* fits = NULL;
    int failed = 0;

    memset(bytes, ' ', sizeof(bytes));
    set_card(bytes, 0U, "SIMPLE  =                    T");
    set_card(bytes, 1U, "BITPIX  =                    8");
    set_card(bytes, 2U, "NAXIS   =                    2");
    set_card(bytes, 3U, "NAXIS1  =                    0");
    set_card(bytes, 4U, "NAXIS2  =                    2");
    set_card(bytes, 5U, "GROUPS  =                    T");
    set_card(bytes, 6U, "PCOUNT  =                    1");
    set_card(bytes, 7U, "GCOUNT  =                    2");
    set_card(bytes, 8U, "EXTEND  =                    T");
    set_card(bytes, 9U, "END");
    set_empty_image_extension(extension);

    if (open_bytes(bytes, sizeof(bytes), &fits) || !fits) {
        return 1;
    }
    if (anqfits_n_ext(fits) != 2) {
        failed++;
    } else {
        failed += fits->exts[0].data_size != 1;
        failed += fits->exts[1].hdr_start != 2;
    }
    anqfits_close(fits);
    return failed;
}

static int test_extension_growth(void) {
    const size_t extension_count = 1024U;
    const size_t block_count = extension_count + 1U;
    const size_t byte_count = block_count * FITS_BLOCK_SIZE;
    unsigned char* bytes = calloc(1U, byte_count);
    anqfits_t* fits = NULL;
    size_t extension;
    int failed = 0;

    if (!bytes) {
        return 1;
    }
    memset(bytes, ' ', byte_count);
    set_card(bytes, 0U, "SIMPLE  =                    T");
    set_card(bytes, 1U, "BITPIX  =                    8");
    set_card(bytes, 2U, "NAXIS   =                    0");
    set_card(bytes, 3U, "EXTEND  =                    T");
    set_card(bytes, 4U, "END");
    for (extension = 0U; extension < extension_count; extension++) {
        unsigned char* block = bytes +
            (extension + 1U) * FITS_BLOCK_SIZE;

        set_empty_image_extension(block);
    }

    if (open_bytes(bytes, byte_count, &fits) || !fits) {
        free(bytes);
        return 1;
    }
    failed += anqfits_n_ext(fits) != (int)block_count;
    for (extension = 0U;
         extension < (size_t)anqfits_n_ext(fits);
         extension++) {
        failed += fits->exts[extension].header == NULL;
        failed += fits->exts[extension].table != NULL;
        failed += fits->exts[extension].image != NULL;
        failed += fits->exts[extension].hdr_start != (int)extension;
        failed += fits->exts[extension].data_start != (int)extension + 1;
        failed += fits->exts[extension].hdr_size != 1;
        failed += fits->exts[extension].data_size != 0;
    }
    anqfits_close(fits);
    free(bytes);
    return failed;
}

static int malformed_data_size_rejected(const char* const* cards) {
    unsigned char bytes[FITS_BLOCK_SIZE];
    anqfits_t* fits = NULL;
    size_t card = 0U;

    memset(bytes, ' ', sizeof(bytes));
    while (cards[card]) {
        set_card(bytes, card, cards[card]);
        card++;
    }
    set_card(bytes, card, "END");
    if (open_bytes(bytes, sizeof(bytes), &fits)) {
        return 1;
    }
    if (fits) {
        anqfits_close(fits);
        return 1;
    }
    return 0;
}

static int test_malformed_data_sizes(void) {
    static const char* missing_bitpix[] = {
        "SIMPLE  =                    T",
        "NAXIS   =                    0",
        NULL
    };
    static const char* invalid_bitpix[] = {
        "SIMPLE  =                    T",
        "BITPIX  =                   24",
        "NAXIS   =                    0",
        NULL
    };
    static const char* missing_axis[] = {
        "SIMPLE  =                    T",
        "BITPIX  =                    8",
        "NAXIS   =                    2",
        "NAXIS2  =                  100",
        NULL
    };
    static const char* zero_gcount[] = {
        "SIMPLE  =                    T",
        "BITPIX  =                    8",
        "NAXIS   =                    0",
        "GCOUNT  =                    0",
        NULL
    };
    static const char* negative_pcount[] = {
        "SIMPLE  =                    T",
        "BITPIX  =                    8",
        "NAXIS   =                    0",
        "PCOUNT  =                   -1",
        NULL
    };
    static const char* overflowing_axes[] = {
        "SIMPLE  =                    T",
        "BITPIX  =                    8",
        "NAXIS   =                    3",
        "NAXIS1  =           2147483647",
        "NAXIS2  =           2147483647",
        "NAXIS3  =           2147483647",
        NULL
    };
    static const char* unrepresentable_axis[] = {
        "SIMPLE  =                    T",
        "BITPIX  =                    8",
        "NAXIS   =                    1",
        "NAXIS1  = 999999999999999999999999999999",
        NULL
    };
    static const char* truncated_primary[] = {
        "SIMPLE  =                    T",
        "BITPIX  =                    8",
        "NAXIS   =                    1",
        "NAXIS1  =                    1",
        NULL
    };
    static const char* truncated_primary_extend_false[] = {
        "SIMPLE  =                    T",
        "BITPIX  =                    8",
        "NAXIS   =                    1",
        "NAXIS1  =                    1",
        "EXTEND  =                    F",
        NULL
    };

    return malformed_data_size_rejected(missing_bitpix) +
        malformed_data_size_rejected(invalid_bitpix) +
        malformed_data_size_rejected(missing_axis) +
        malformed_data_size_rejected(zero_gcount) +
        malformed_data_size_rejected(negative_pcount) +
        malformed_data_size_rejected(overflowing_axes) +
        malformed_data_size_rejected(unrepresentable_axis) +
        malformed_data_size_rejected(truncated_primary) +
        malformed_data_size_rejected(truncated_primary_extend_false);
}

static int malformed_extension_rejected(const char* const* cards) {
    unsigned char bytes[2U * FITS_BLOCK_SIZE];
    unsigned char* extension = bytes + FITS_BLOCK_SIZE;
    anqfits_t* fits = NULL;
    size_t card = 0U;

    memset(bytes, ' ', sizeof(bytes));
    set_card(bytes, 0U, "SIMPLE  =                    T");
    set_card(bytes, 1U, "BITPIX  =                    8");
    set_card(bytes, 2U, "NAXIS   =                    0");
    set_card(bytes, 3U, "EXTEND  =                    T");
    set_card(bytes, 4U, "END");
    while (cards[card]) {
        set_card(extension, card, cards[card]);
        card++;
    }
    set_card(extension, card, "END");
    if (open_bytes(bytes, sizeof(bytes), &fits)) {
        return 1;
    }
    if (fits) {
        anqfits_close(fits);
        return 1;
    }
    return 0;
}

static int test_malformed_extensions(void) {
    static const char* negative_pcount[] = {
        "XTENSION= 'IMAGE   '",
        "BITPIX  =                    8",
        "NAXIS   =                    0",
        "PCOUNT  =                   -1",
        "GCOUNT  =                    1",
        NULL
    };
    static const char* truncated_data[] = {
        "XTENSION= 'IMAGE   '",
        "BITPIX  =                    8",
        "NAXIS   =                    1",
        "NAXIS1  =                    1",
        "PCOUNT  =                    0",
        "GCOUNT  =                    1",
        NULL
    };
    static const char* extension_groups[] = {
        "XTENSION= 'IMAGE   '",
        "BITPIX  =                    8",
        "NAXIS   =                    2",
        "NAXIS1  =                    0",
        "NAXIS2  =                    1",
        "GROUPS  =                    T",
        "PCOUNT  =                    0",
        "GCOUNT  =                    1",
        NULL
    };

    return malformed_extension_rejected(negative_pcount) +
        malformed_extension_rejected(truncated_data) +
        malformed_extension_rejected(extension_groups);
}

static int test_extension_without_end(void) {
    unsigned char bytes[2U * FITS_BLOCK_SIZE];
    unsigned char* extension = bytes + FITS_BLOCK_SIZE;
    anqfits_t* fits = NULL;

    memset(bytes, ' ', sizeof(bytes));
    set_card(bytes, 0U, "SIMPLE  =                    T");
    set_card(bytes, 1U, "BITPIX  =                    8");
    set_card(bytes, 2U, "NAXIS   =                    0");
    set_card(bytes, 3U, "EXTEND  =                    T");
    set_card(bytes, 4U, "END");
    set_card(extension, 0U, "XTENSION= 'IMAGE   '");
    set_card(extension, 1U, "BITPIX  =                    8");
    set_card(extension, 2U, "NAXIS   =                    0");
    set_card(extension, 3U, "PCOUNT  =                    0");
    set_card(extension, 4U, "GCOUNT  =                    1");
    if (open_bytes(bytes, sizeof(bytes), &fits)) {
        return 1;
    }
    if (fits) {
        anqfits_close(fits);
        return 1;
    }
    return 0;
}

int main(void) {
    int failures = 0;

    failures += test_header_card_block_boundary();
    failures += test_malformed_keys();
    failures += test_zero_length_first_axis();
    failures += test_random_groups_data_extent();
    failures += test_extension_growth();
    failures += test_malformed_data_sizes();
    failures += test_malformed_extensions();
    failures += test_extension_without_end();
    if (failures) {
        fprintf(stderr, "test_anqfits_header: %i failures\n", failures);
        return 1;
    }
    printf("test_anqfits_header: PASS\n");
    return 0;
}
