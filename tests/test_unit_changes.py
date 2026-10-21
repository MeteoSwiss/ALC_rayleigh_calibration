"""Instrument-unit change detection for the dashboard (monitoring/unit_changes.py)."""
import json

import pandas as pd

from monitoring.unit_changes import current_identity, unit_changes


def _hk(tmp_path, key, rows):
    d = tmp_path / key
    d.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(d / f"{key}_hk.csv", index=False)


def test_module_and_firmware_changes_detected(tmp_path):
    key = "0-20000-0-00001_A"
    _hk(tmp_path, key, [
        {"date": "20250101", "laser": 100, "unit_id": "TUB200009", "firmware": "1.13"},
        {"date": "20250102", "laser": 100, "unit_id": "", "firmware": ""},          # gap: ignored
        {"date": "20250402", "laser": 100, "unit_id": "TUB140016", "firmware": "1.13"},
        {"date": "20260509", "laser": 100, "unit_id": "TUB140016", "firmware": "1.15"},
    ])
    ch = unit_changes(tmp_path, key, events_file=tmp_path / "none.json")
    kinds = [(c["date"].strftime("%Y%m%d"), c["kind"]) for c in ch]
    assert kinds == [("20250402", "unit"), ("20260509", "firmware")]
    assert "TUB200009" in ch[0]["label"] and "TUB140016" in ch[0]["label"]
    ident = current_identity(tmp_path, key)
    assert ident["unit_id"] == "TUB140016" and ident["firmware"] == "1.15"


def test_placeholder_serial_is_not_a_change_and_events_are_listed(tmp_path):
    key = "0-20000-0-00002_B"
    _hk(tmp_path, key, [
        {"date": "20260601", "laser": 99, "unit_id": "CL0", "firmware": "205"},
        {"date": "20260708", "laser": 99, "unit_id": "CL0", "firmware": "205"},
    ])
    ev = tmp_path / "events.json"
    ev.write_text(json.dumps({key: [{"date": "2026-07-07", "label": "optical block replaced"}]}),
                  encoding="utf-8")
    ch = unit_changes(tmp_path, key, events_file=ev)
    assert [c["kind"] for c in ch] == ["event"]
    ident = current_identity(tmp_path, key)
    assert ident["placeholder"] is True and ident["unit_id"] == "CL0"


def test_tub_prefix_and_bare_number_are_the_same_unit(tmp_path):
    key = "0-20000-0-00003_A"
    _hk(tmp_path, key, [
        {"date": "20250101", "laser": 100, "unit_id": "TUB140016", "firmware": "1.13"},
        {"date": "20250105", "laser": 100, "unit_id": "140016", "firmware": "1.13"},
    ])
    assert unit_changes(tmp_path, key, events_file=tmp_path / "none.json") == []
