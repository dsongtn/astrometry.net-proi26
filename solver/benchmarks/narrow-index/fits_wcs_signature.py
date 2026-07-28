#!/usr/bin/env python3
"""Canonical scientific signatures for Astrometry.net FITS WCS headers.

This module intentionally uses only the Python standard library.  It parses the
primary FITS header directly, selects the coordinate transformation cards, and
normalizes their typed values before hashing.  Formatting, card order, comments,
HISTORY, DATE and other non-WCS metadata therefore cannot change the signature.
"""

from __future__ import annotations

from decimal import Decimal, InvalidOperation
import hashlib
import json
from pathlib import Path
import re
from typing import Any, Dict, List, Tuple


SIGNATURE_METHOD = "fits-wcs-v1"
FITS_CARD_BYTES = 80
FITS_BLOCK_BYTES = 2880

_NUMBER_RE = re.compile(
    r"^[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[DE][+-]?[0-9]+)?$",
    re.IGNORECASE,
)
_STANDARD_WCS_RE = re.compile(
    r"^(?:"
    r"WCSAXES[A-Z]?|WCSNAME[A-Z]?|"
    r"CTYPE[0-9]+[A-Z]?|CUNIT[0-9]+[A-Z]?|CNAME[0-9]+[A-Z]?|"
    r"CRPIX[0-9]+[A-Z]?|CRVAL[0-9]+[A-Z]?|CDELT[0-9]+[A-Z]?|"
    r"CRDER[0-9]+[A-Z]?|CSYER[0-9]+[A-Z]?|CROTA[0-9]+[A-Z]?|"
    r"PC[0-9]+_[0-9]+[A-Z]?|CD[0-9]+_[0-9]+[A-Z]?|"
    r"PV[0-9]+_[0-9]+[A-Z]?|PS[0-9]+_[0-9]+[A-Z]?|"
    r"LONPOLE[A-Z]?|LATPOLE[A-Z]?|RADESYS[A-Z]?|RADECSYS|EQUINOX[A-Z]?"
    r")$"
)
_SIP_RE = re.compile(r"^(?:A|B|AP|BP)_(?:ORDER|[0-9]+_[0-9]+)$")
_IMAGE_GEOMETRY = frozenset({"NAXIS1", "NAXIS2", "IMAGEW", "IMAGEH"})


class WcsSignatureError(ValueError):
    """The FITS header cannot provide an unambiguous scientific WCS signature."""


def _header_cards(path: Path) -> List[str]:
    cards: List[str] = []
    try:
        with path.open("rb") as stream:
            while True:
                block = stream.read(FITS_BLOCK_BYTES)
                if len(block) != FITS_BLOCK_BYTES:
                    raise WcsSignatureError(
                        f"{path}: primary FITS header ended before an END card"
                    )
                try:
                    text = block.decode("ascii")
                except UnicodeDecodeError as error:
                    raise WcsSignatureError(
                        f"{path}: primary FITS header is not ASCII"
                    ) from error
                for offset in range(0, FITS_BLOCK_BYTES, FITS_CARD_BYTES):
                    card = text[offset : offset + FITS_CARD_BYTES]
                    keyword = card[:8].strip().upper()
                    if keyword == "END":
                        return cards
                    cards.append(card)
    except OSError as error:
        raise WcsSignatureError(f"cannot read FITS WCS {path}: {error}") from error


def _quoted_string(field: str, path: Path, keyword: str) -> str:
    value: List[str] = []
    position = 1
    while position < len(field):
        character = field[position]
        if character != "'":
            value.append(character)
            position += 1
            continue
        if position + 1 < len(field) and field[position + 1] == "'":
            value.append("'")
            position += 2
            continue
        remainder = field[position + 1 :].lstrip()
        if remainder and not remainder.startswith("/"):
            raise WcsSignatureError(
                f"{path}: malformed trailing text in FITS string card {keyword}"
            )
        return "".join(value).rstrip()
    raise WcsSignatureError(f"{path}: unterminated FITS string card {keyword}")


def _unquoted_value(field: str) -> str:
    return field.split("/", 1)[0].strip()


def _canonical_decimal(value: str, path: Path, keyword: str) -> str:
    try:
        decimal = Decimal(value.replace("D", "E").replace("d", "e"))
    except InvalidOperation as error:
        raise WcsSignatureError(
            f"{path}: invalid numeric value {value!r} in {keyword}"
        ) from error
    if not decimal.is_finite():
        raise WcsSignatureError(f"{path}: non-finite numeric value in {keyword}")
    if decimal.is_zero():
        return "0"
    sign, digits_tuple, exponent = decimal.as_tuple()
    digits = list(digits_tuple)
    while digits and digits[-1] == 0:
        digits.pop()
        exponent += 1
    prefix = "-" if sign else ""
    return f"{prefix}{''.join(str(digit) for digit in digits)}e{exponent}"


def _typed_value(card: str, path: Path, keyword: str) -> Dict[str, Any]:
    if card[8:10] != "= ":
        raise WcsSignatureError(f"{path}: WCS card {keyword} has no FITS value")
    field = card[10:]
    stripped = field.lstrip()
    if not stripped:
        raise WcsSignatureError(f"{path}: WCS card {keyword} has an empty value")
    if stripped.startswith("'"):
        return {"type": "string", "value": _quoted_string(stripped, path, keyword)}
    value = _unquoted_value(stripped)
    if value in ("T", "F"):
        return {"type": "boolean", "value": value == "T"}
    if _NUMBER_RE.fullmatch(value):
        return {
            "type": "number",
            "value": _canonical_decimal(value, path, keyword),
        }
    raise WcsSignatureError(
        f"{path}: unsupported FITS value {value!r} in WCS card {keyword}"
    )


def _is_signature_keyword(keyword: str) -> bool:
    return bool(
        keyword in _IMAGE_GEOMETRY
        or _STANDARD_WCS_RE.fullmatch(keyword)
        or _SIP_RE.fullmatch(keyword)
    )


def canonical_wcs_payload(path: Path) -> Dict[str, Any]:
    """Return normalized WCS cards suitable for stable JSON serialization."""
    selected: List[Tuple[str, Dict[str, Any]]] = []
    seen: set[str] = set()
    for card in _header_cards(path):
        keyword = card[:8].strip().upper()
        if not _is_signature_keyword(keyword):
            continue
        if keyword in seen:
            raise WcsSignatureError(f"{path}: duplicate WCS card {keyword}")
        seen.add(keyword)
        selected.append((keyword, _typed_value(card, path, keyword)))

    required = {"CTYPE1", "CTYPE2", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2"}
    missing = sorted(required - seen)
    if missing:
        raise WcsSignatureError(f"{path}: missing TAN core WCS cards {missing}")
    ctype_values = {
        keyword: value.get("value")
        for keyword, value in selected
        if keyword in ("CTYPE1", "CTYPE2")
    }
    if any("TAN" not in str(ctype_values[keyword]).upper() for keyword in ("CTYPE1", "CTYPE2")):
        raise WcsSignatureError(f"{path}: CTYPE1/CTYPE2 do not describe a TAN projection")
    has_cd = all(keyword in seen for keyword in ("CD1_1", "CD1_2", "CD2_1", "CD2_2"))
    has_cdelt = all(keyword in seen for keyword in ("CDELT1", "CDELT2"))
    if not has_cd and not has_cdelt:
        raise WcsSignatureError(
            f"{path}: WCS has neither a complete 2x2 CD matrix nor CDELT1/CDELT2"
        )
    has_sip = any("SIP" in str(ctype_values[keyword]).upper() for keyword in ("CTYPE1", "CTYPE2"))
    if has_sip and not {"A_ORDER", "B_ORDER"}.issubset(seen):
        raise WcsSignatureError(f"{path}: SIP projection lacks A_ORDER/B_ORDER")
    for family in ("A", "B", "AP", "BP"):
        has_coefficients = any(
            re.fullmatch(rf"{family}_[0-9]+_[0-9]+", keyword) for keyword in seen
        )
        if has_coefficients and f"{family}_ORDER" not in seen:
            raise WcsSignatureError(
                f"{path}: {family} SIP coefficients lack {family}_ORDER"
            )

    return {
        "method": SIGNATURE_METHOD,
        "cards": [
            {"keyword": keyword, **value}
            for keyword, value in sorted(selected, key=lambda item: item[0])
        ],
    }


def canonical_wcs_bytes(path: Path) -> bytes:
    return (
        json.dumps(
            canonical_wcs_payload(path),
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
        )
        + "\n"
    ).encode("ascii")


def wcs_signature(path: Path) -> str:
    return hashlib.sha256(canonical_wcs_bytes(path)).hexdigest()
