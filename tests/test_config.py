import unittest
import threading
from dataclasses import FrozenInstanceError

from py_hla_match.policy import (
    ExpressionSuffixPolicy,
    ExpressionSuffixMatchLevel
)
from py_hla_match.config import (
    HLAMatchConfig,
    CANONICAL_HLA_LOCI,
    get_config,
    get_config_version,
    set_config,
    config_context
)


class TestConfigDefaults(unittest.TestCase):
    """Basic defaults and invariants."""

    def tearDown(self):
        # Reset to library defaults to avoid side effects across tests
        set_config(HLAMatchConfig())

    def test_default_config_values(self):
        config = get_config()
        self.assertTrue(config.strict_loci)
        self.assertIsInstance(
            config.expression_suffix_policy, ExpressionSuffixPolicy
        )
        self.assertIsNotNone(config.nomenclature_pattern)
        self.assertIsInstance(get_config_version(), int)

    def test_effective_valid_loci_contains_canonical(self):
        config = get_config()
        effective_loci = config.effective_valid_loci
        self.assertIsInstance(effective_loci, frozenset)
        # canonical loci are included
        for locus in ["A", "B", "C", "DRB1", "DQB1", "DPB1", "DRB345", "DRBX"]:
            self.assertIn(locus, effective_loci)
        # canonical set equals subset of effective
        self.assertTrue(set(CANONICAL_HLA_LOCI).issubset(set(effective_loci)))


class TestConfigSetBehavior(unittest.TestCase):
    """set_config warnings/errors and version bump."""

    def tearDown(self):
        set_config(HLAMatchConfig())

    def test_strict_loci_rejects_extras(self):
        config = HLAMatchConfig(
            extra_valid_loci=frozenset({"DRA"}), strict_loci=True
        )
        with self.assertRaises(ValueError):
            set_config(config)

    def test_additive_extras_warn_and_update_effective(self):
        prev_ver = get_config_version()
        config = HLAMatchConfig(
            extra_valid_loci=frozenset({"DRA"}), strict_loci=False
        )
        with self.assertLogs("py_hla_match.config", level="WARNING") as cm:
            set_config(config)
        self.assertTrue(any("extra_valid_loci" in msg for msg in cm.output))
        self.assertIn("DRA", get_config().effective_valid_loci)
        self.assertEqual(get_config_version(), prev_ver + 1)


class TestExpressionSuffixPolicy(unittest.TestCase):
    """Explicit policy decisions remain configurable."""

    def tearDown(self):
        set_config(HLAMatchConfig())

    def test_policy_defaults(self):
        policy = get_config().expression_suffix_policy
        self.assertEqual(
            policy.equal_risk, ExpressionSuffixMatchLevel.NOT_APPLICABLE
        )
        self.assertEqual(
            policy.risk_vs_none, ExpressionSuffixMatchLevel.ALLELE_MISMATCH
        )
        self.assertEqual(
            policy.risk_vs_different_risk,
            ExpressionSuffixMatchLevel.ALLELE_MISMATCH
        )
        self.assertEqual(
            policy.q_present, ExpressionSuffixMatchLevel.NOT_APPLICABLE
        )

    def test_set_explicit_policy(self):
        policy = ExpressionSuffixPolicy(
            equal_risk=ExpressionSuffixMatchLevel.NOT_APPLICABLE,
            risk_vs_none=ExpressionSuffixMatchLevel.ALLELE_MISMATCH,
            risk_vs_different_risk=ExpressionSuffixMatchLevel.ALLELE_MISMATCH,
            q_present=ExpressionSuffixMatchLevel.NOT_APPLICABLE,
        )
        config = HLAMatchConfig(expression_suffix_policy=policy)
        set_config(config)
        self.assertIs(get_config().expression_suffix_policy, policy)


class TestConfigContextAndThreadLocal(unittest.TestCase):
    def tearDown(self):
        set_config(HLAMatchConfig())

    def test_scoped_override_and_restore(self):
        base = HLAMatchConfig()
        set_config(base)
        self.assertTrue(get_config().strict_loci)
        temp = HLAMatchConfig(strict_loci=False)
        with config_context(temp):
            self.assertFalse(get_config().strict_loci)
        # restored
        self.assertTrue(get_config().strict_loci)

    def test_nested_contexts_restore(self):
        base = HLAMatchConfig(na_tokens=frozenset({"NA", "NE"}))
        set_config(base)
        outer = HLAMatchConfig(na_tokens=frozenset({"NA"}))
        inner = HLAMatchConfig(na_tokens=frozenset({"NA", "NONE"}))
        with config_context(outer):
            self.assertEqual(get_config().na_tokens, frozenset({"NA"}))
            with config_context(inner):
                self.assertEqual(
                    get_config().na_tokens, frozenset({"NA", "NONE"})
                )
            # restored to outer
            self.assertEqual(get_config().na_tokens, frozenset({"NA"}))
        # restored to base
        self.assertEqual(get_config().na_tokens, frozenset({"NA", "NE"}))

    def test_thread_local_isolation(self):
        # main thread remains with strict_loci=True
        set_config(HLAMatchConfig(strict_loci=True))
        result = {"child_strict": None}

        def worker():
            set_config(HLAMatchConfig(strict_loci=False))
            result["child_strict"] = get_config().strict_loci

        t = threading.Thread(target=worker)
        t.start()
        t.join()

        # Child thread had False, main thread stays True
        self.assertFalse(result["child_strict"])
        self.assertTrue(get_config().strict_loci)


class TestConfigImmutability(unittest.TestCase):
    def tearDown(self):
        set_config(HLAMatchConfig())

    def test_get_config_is_not_mutable(self):
        cfg = get_config()
        with self.assertRaises(FrozenInstanceError):
            # direct mutation must be disallowed; use set_config instead
            cfg.strict_loci = False
