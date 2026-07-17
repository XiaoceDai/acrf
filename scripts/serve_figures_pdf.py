#!/usr/bin/env python3
"""Serve figures/ over HTTP so Cursor Simple Browser can render PDFs remotely.

Chromium (Cursor built-in / Glass browser) can display application/pdf natively.
Custom PDF editor extensions often fail on Cursor 3 + Remote-SSH; this is the
reliable remote preview path.

Usage:
  python3 scripts/serve_figures_pdf.py          # start on :8765
  python3 scripts/serve_figures_pdf.py --port 8765 --check
"""
from __future__ import annotations

import argparse
import json
import socket
import time
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "figures"
DEBUG_LOG = ROOT / ".cursor" / "debug-5793e0.log"
DEFAULT_PORT = 8765


def _log(hypothesis_id: str, message: str, data: dict) -> None:
    # #region agent log
    DEBUG_LOG.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "sessionId": "5793e0",
        "runId": "post-fix",
        "hypothesisId": hypothesis_id,
        "location": "serve_figures_pdf.py",
        "message": message,
        "data": data,
        "timestamp": int(time.time() * 1000),
    }
    with DEBUG_LOG.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(payload) + "\n")
    # #endregion


class FiguresHandler(SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(FIG), **kwargs)

    def end_headers(self):
        # Allow Simple Browser / Glass to embed PDFs.
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Cache-Control", "no-cache")
        super().end_headers()

    def log_message(self, fmt: str, *args) -> None:
        print(f"[figures-pdf] {self.address_string()} - {fmt % args}")


def port_open(port: int) -> bool:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.settimeout(0.3)
        return s.connect_ex(("127.0.0.1", port)) == 0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", type=int, default=DEFAULT_PORT)
    ap.add_argument("--check", action="store_true", help="Only probe port/PDF headers")
    args = ap.parse_args()

    sample = FIG / "fig6a_anchor_coloc_calibration.pdf"
    _log(
        "F",
        "serve_figures_probe",
        {
            "port": args.port,
            "port_open": port_open(args.port),
            "sample_exists": sample.is_file(),
            "figures_dir": str(FIG),
        },
    )

    if args.check:
        print(f"port_open={port_open(args.port)} sample={sample.is_file()}")
        print(f"open in Cursor browser: http://127.0.0.1:{args.port}/{sample.name}")
        return

    if port_open(args.port):
        print(f"already listening on :{args.port}")
        print(f"open: http://127.0.0.1:{args.port}/fig6a_anchor_coloc_calibration.pdf")
        _log("G", "server_already_running", {"port": args.port})
        return

    httpd = ThreadingHTTPServer(("127.0.0.1", args.port), FiguresHandler)
    print(f"Serving {FIG} at http://127.0.0.1:{args.port}/")
    print("Open a PDF in Cursor Simple Browser / Glass, e.g.:")
    print(f"  http://127.0.0.1:{args.port}/fig6a_anchor_coloc_calibration.pdf")
    _log("G", "server_started", {"port": args.port, "dir": str(FIG)})
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nstopped")


if __name__ == "__main__":
    main()
