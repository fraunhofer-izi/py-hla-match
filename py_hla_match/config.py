from __future__ import annotations

import logging
import re
import threading
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import FrozenSet, Optional, Pattern

from py_hla_match.policy import (
    ExpressionSuffixPolicy
)

logger = logging.getLogger(__name__)


CANONICAL_HLA_LOCI: FrozenSet[str] = frozenset(
    {
        'A', 'B', 'C',  # NOTE: standard HLA Class I loci
        'DRB1',  # NOTE: standard HLA Class II locus
        'DQB1',  # NOTE: standard HLA Class II locus
        'DPB1',  # NOTE: not standard but often discussed
        'DQA1',  # NOTE: interesting candidate
        'DPA1',  # NOTE: interesting candidate
        'DRB3', 'DRB4', 'DRB5',  # NOTE: interesting candidate(s)
        'DRB345',  # NOTE: NOT A LOCUS, but used here as generic DRB3/4/5
        'DRBX',  # NOTE: placeholder to indicate missing DRB3/4/5
        # However, DRB3/4/5 may need special handling

        'MICA', 'MICB',  # NOTE: interesting candidates
        'DMA', 'DMB', 'DOA', 'DOB',  # NOTE: interesting candidates
        # the exhaustive list of all known HLA loci
        # should be considered non-standard or experimental

        'E', 'F', 'G',  # limited data
        'DQA2', 'DQB2',  # limited data

        # currently not considered for matching ---

        # 'MICC', 'MICD', 'MICE'  # non functional related to MICA/MICB
        # 'TAP1', 'TAP2',  # not direct cell-surface targets
        # 'PSMB9', 'PSMB8',  # not direct cell-surface targets
        # 'HFE',  # hemochromatosis locus, not relevant to HLA matching
        # 'DRA',  # largely monomorphic/minimally polymorphic DR alpha chain
        # likely pseudogenes ---
        # 'DPA2', 'DPB2', 'DPA3', 'DQB3',
        # --- likely pseudogenes
        # pseudogenes ---
        # 'DRB2', 'DRB6', 'DRB7', 'DRB8', 'DRB9',
        # 'H', 'J', 'K', 'L', 'N', 'P', 'R', 'S', 'T',
        # 'U', 'V', 'W', 'X', 'Y', 'Z',
        # --- pseudogenes

        # --- currently not considered for matching
    }
)

# we now support processing of DRB345
CANONICAL_DRB345_SUB_LOCI = {
    'DRB3', 'DRB4', 'DRB5',
    'DRBX'  # generic locus indicating "missing"
}


@dataclass(frozen=True)
class HLAMatchConfig:
    # ARD configuration (read once at singleton creation)
    ard_imgt_version: str = "Latest"
    ard_data_dir: Optional[str] = None

    # Extensions only, canonical set remains owned by the library
    extra_valid_loci: FrozenSet[str] = field(default_factory=frozenset)
    strict_loci: bool = True  # must be set to False, by user

    # NA tokens normalized across datasets
    na_tokens: FrozenSet[str] = field(
        default_factory=lambda: frozenset(
            {"NA", "NE", "NEW", "UNKNOWN", "ND", "NULL"}
        )
    )

    # Expression suffix policy
    expression_suffix_policy: ExpressionSuffixPolicy = field(
        default_factory=ExpressionSuffixPolicy
    )

    # Internal caches
    _nomenclature_pattern: Optional[Pattern[str]] = field(
        default=None, init=False, repr=False
    )
    _effective_valid_loci: Optional[FrozenSet[str]] = field(
        default=None, init=False, repr=False
    )

    @property
    def effective_valid_loci(self) -> FrozenSet[str]:
        if self._effective_valid_loci is None:
            object.__setattr__(
                self,
                "_effective_valid_loci",
                frozenset(
                    set(CANONICAL_HLA_LOCI) | set(self.extra_valid_loci)
                ),
            )
        return self._effective_valid_loci

    @property
    def drb345_sub_loci(self) -> FrozenSet[str]:
        """
        Placeholder: canonical DRB345 sub-loci (DRB3/DRB4/DRB5/DRBX).
        Exposed via config for future adjustability without changing hla
        parser logic.
        """
        return CANONICAL_DRB345_SUB_LOCI

    def _compile_nomenclature_pattern(self) -> Pattern[str]:
        # Build regex with configured NA tokens; match IMGT/HLA field, suffix,
        # and group code constraints.
        na_alt = "|".join(sorted(re.escape(tok) for tok in self.na_tokens))
        pat = rf"""
        # optional prefix
        ^ (?:HLA-)?
        # locus
        (?P<locus>[A-Z0-9]+)
        # asterisk
        \*
        # three alternatives
        (
            # 1–4 fields
            (?P<allele_fields>\d{{2,4}}(?::\d{{2,4}}){{0,3}})
            # suffix AND group (error))
            (?![NLSCAQ][GP]|[GP][NLSCAQ])
            # optional expression suffix
            (?P<suffix>[NLSCAQ])?
            # optional group code
            (?P<group_code>[GP])?
            $
        |
            # configured NA tokens
            (?P<nan>{na_alt})
            $
        |
            # anything else for error message
            (?P<remainder>.*)
            $
        )
        """
        return re.compile(pat, re.VERBOSE)

    @property
    def nomenclature_pattern(self) -> Pattern[str]:
        if self._nomenclature_pattern is None:
            self._nomenclature_pattern = self._compile_nomenclature_pattern()
        return self._nomenclature_pattern

    def recompile_patterns(self) -> None:
        object.__setattr__(
            self, "_nomenclature_pattern", self._compile_nomenclature_pattern()
        )
        # Invalidate effective loci cache in case config changed
        object.__setattr__(self, "_effective_valid_loci", None)


# local config storage from global defaults
_THREADLOCAL = threading.local()
_GLOBAL_DEFAULTS: HLAMatchConfig = HLAMatchConfig()
# global version token for cache invalidation
_CONFIG_VERSION: int = 1


def _get_threadlocal_config() -> HLAMatchConfig:
    """
    Get (thread-)local HLAMatchConfig, from _GLOBAL_DEFAULTS on first access
    for thread.
    """
    if not hasattr(_THREADLOCAL, "config"):
        # new instance seeded from defaults
        _THREADLOCAL.config = HLAMatchConfig(
            ard_imgt_version=_GLOBAL_DEFAULTS.ard_imgt_version,
            ard_data_dir=_GLOBAL_DEFAULTS.ard_data_dir,
            extra_valid_loci=_GLOBAL_DEFAULTS.extra_valid_loci,
            strict_loci=_GLOBAL_DEFAULTS.strict_loci,
            na_tokens=_GLOBAL_DEFAULTS.na_tokens,
            expression_suffix_policy=_GLOBAL_DEFAULTS.expression_suffix_policy,
        )
        _THREADLOCAL.config.recompile_patterns()
    return _THREADLOCAL.config


def _set_threadlocal_config(config: HLAMatchConfig) -> None:
    _THREADLOCAL.config = config


def get_config() -> HLAMatchConfig:
    return _get_threadlocal_config()


def get_config_version() -> int:
    return _CONFIG_VERSION


def set_config(config: HLAMatchConfig) -> None:
    global _CONFIG_VERSION
    previous_config = _get_threadlocal_config()
    if config.extra_valid_loci:
        if config.strict_loci:
            raise ValueError(
                "extra_valid_loci not allowed in strict mode: "
                f"{sorted(config.extra_valid_loci)}"
            )
        logger.warning(
            "Using extra_valid_loci (additive to canonical set): %s",
            sorted(config.extra_valid_loci),
        )

    # new (thread-)local config
    _set_threadlocal_config(config)
    config.recompile_patterns()
    # bump version to signal cache invalidation
    _CONFIG_VERSION += 1

    # reset singleton if ARD defaults changed
    try:
        if (
            previous_config.ard_imgt_version != config.ard_imgt_version
            or previous_config.ard_data_dir != config.ard_data_dir
        ):
            # late import to prevent circular import at module import
            from .singleton import _reset_ard_instance
            _reset_ard_instance()
    except Exception:
        # log instead of fail ARD set_config on reset issues
        logger.debug("Optional ARD singleton reset failed.", exc_info=True)


@contextmanager
def config_context(config: HLAMatchConfig):
    """
    Context manager to temporarily override (thread-)local config.
    The previous (potentially default) config is restored on exit or exception.
    """
    previous_config = _get_threadlocal_config()
    set_config(config)
    try:
        yield
    finally:
        set_config(previous_config)
