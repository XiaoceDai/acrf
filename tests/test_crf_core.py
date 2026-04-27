"""Regression tests for the portable CRF method core (must match figure scripts)."""

from __future__ import annotations

import pytest

from crf import (
    TriangulationInput,
    acai,
    acai_rel,
    bulk_strength_from_labels,
    cai,
    dur_direction,
    fig2c_trend_state,
    interpret_acai_rel,
    ncp_baseline_sign,
    triangulate,
)


def test_acai_rel_identity_on_healthy_reference() -> None:
    assert acai_rel(1.5, 1.5) == pytest.approx(1.0)
    assert interpret_acai_rel(1.1) == "amplified"
    assert interpret_acai_rel(0.5) == "attenuated"
    assert interpret_acai_rel(-0.1) == "reversed"


def test_ncp_sign_and_acai() -> None:
    assert ncp_baseline_sign("LA-biased") == 1
    assert ncp_baseline_sign("RA-biased") == -1
    raw = cai(0.3, 0.1)
    assert raw == pytest.approx(0.2)
    assert acai(raw, 1) == pytest.approx(0.2)
    assert acai(raw, -1) == pytest.approx(-0.2)


def test_bulk_strength_matches_fig5a_table() -> None:
    assert bulk_strength_from_labels("Reversed", "")[0] == -2
    assert bulk_strength_from_labels("Amplified", "")[0] == 2
    assert bulk_strength_from_labels("", "Preserved")[1] == "Preserved"


def test_dur_direction_threshold() -> None:
    assert dur_direction(0.29, thr=0.30) == 0
    assert dur_direction(0.30, thr=0.30) == 1
    assert dur_direction(-0.31, thr=0.30) == -1


def test_fig2c_trend_state_edges() -> None:
    assert fig2c_trend_state(-0.03, 1.0) == "Reversed"
    assert fig2c_trend_state(0.5, 1.0) == "Attenuated"  # ratio 0.5 < 0.67
    assert fig2c_trend_state(0.8, 1.0) == "Preserved"
    assert fig2c_trend_state(1.4, 1.0) == "Amplified"


def test_triangulation_metabolism_led() -> None:
    row = TriangulationInput(
        bulk_strength=-1,
        dur_dir=-1,
        baseline_direction=1,
        af_major_carrier_ct="Unresolved",
    )
    assert triangulate(row) == "Metabolism-led"


def test_triangulation_structure_led_fibroblast() -> None:
    row = TriangulationInput(
        bulk_strength=-2,
        dur_dir=0,
        baseline_direction=-1,
        af_major_carrier_ct="fibroblast",
    )
    assert triangulate(row) == "Structure-led"


def test_triangulation_convergent_struct_arm() -> None:
    row = TriangulationInput(
        bulk_strength=-1,
        dur_dir=1,
        baseline_direction=-1,
        af_major_carrier_ct="Unresolved",
    )
    assert triangulate(row) == "Convergent"
