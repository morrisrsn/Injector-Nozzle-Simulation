from dataclasses import asdict, fields
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse
import argparse
import json
import subprocess
import sys
import threading
import warnings
import webbrowser

from Thermodynamics import NozzleInput, calculate_nozzle


ROOT = Path(__file__).resolve().parent
FRONTEND_FILE = ROOT / "frontend" / "index.html"
INPUT_DEFAULTS = asdict(NozzleInput())
INPUT_FIELDS = {field.name: field for field in fields(NozzleInput)}


def _coerce_config(payload):
    config = {}

    for name, default in INPUT_DEFAULTS.items():
        if name not in payload:
            continue

        value = payload[name]
        if isinstance(default, float):
            config[name] = float(value)
        elif isinstance(default, int):
            config[name] = int(value)
        else:
            config[name] = str(value)

    return config


def _send_bytes(handler, body, content_type, status=200):
    handler.send_response(status)
    handler.send_header("Content-Type", content_type)
    handler.send_header("Content-Length", str(len(body)))
    handler.end_headers()
    handler.wfile.write(body)


def _send_json(handler, payload, status=200):
    body = json.dumps(payload, ensure_ascii=False).encode("utf-8")
    _send_bytes(handler, body, "application/json; charset=utf-8", status)


def _open_browser(url):
    opened = webbrowser.open(url)
    if not opened and sys.platform == "darwin":
        subprocess.Popen(["open", url])


class NozzleRequestHandler(BaseHTTPRequestHandler):
    server_version = "NozzleFrontend/1.0"

    def log_message(self, format, *args):
        print("%s - - %s" % (self.address_string(), format % args))

    def do_GET(self):
        path = urlparse(self.path).path

        if path == "/" or path == "/index.html":
            body = FRONTEND_FILE.read_bytes()
            _send_bytes(self, body, "text/html; charset=utf-8")
            return

        if path == "/api/defaults":
            _send_json(
                self,
                {
                    "defaults": INPUT_DEFAULTS,
                    "fields": list(INPUT_FIELDS.keys()),
                    "options": {
                        "mode": ["combustion", "cold_gas"],
                        "chamber_method": ["A", "B"],
                    },
                },
            )
            return

        _send_json(self, {"error": "Nicht gefunden"}, status=404)

    def do_POST(self):
        path = urlparse(self.path).path
        if path != "/api/calculate":
            _send_json(self, {"error": "Nicht gefunden"}, status=404)
            return

        try:
            length = int(self.headers.get("Content-Length", "0"))
            raw_body = self.rfile.read(length).decode("utf-8") if length else "{}"
            payload = json.loads(raw_body)
            config = _coerce_config(payload)
            num_points = int(payload.get("num_points", 180))

            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                result = calculate_nozzle(config, num_points=num_points)

            result["warnings"] = [str(item.message) for item in caught]
            _send_json(self, result)
        except Exception as exc:
            _send_json(self, {"error": str(exc)}, status=400)


def main():
    parser = argparse.ArgumentParser(description="Lokales HTML-Frontend für Thermodynamics.py")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", default=8000, type=int)
    parser.add_argument("--no-browser", action="store_true", help="Browser beim Start nicht automatisch öffnen")
    args = parser.parse_args()

    server = ThreadingHTTPServer((args.host, args.port), NozzleRequestHandler)
    browser_host = "127.0.0.1" if args.host in ("", "0.0.0.0") else args.host
    url = f"http://{browser_host}:{args.port}"
    print(f"Frontend läuft unter {url}")
    if not args.no_browser:
        threading.Timer(0.5, lambda: _open_browser(url)).start()
    server.serve_forever()


if __name__ == "__main__":
    main()
