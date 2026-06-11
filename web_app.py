from dataclasses import asdict, fields
from datetime import datetime, timezone
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import unquote, urlparse
import argparse
import json
import re
import subprocess
import sys
import threading
import warnings
import webbrowser

from Thermodynamics import NozzleInput, calculate_nozzle


ROOT = Path(__file__).resolve().parent
FRONTEND_FILE = ROOT / "frontend" / "index.html"
CONFIG_DIR = ROOT / "configs"
CONFIG_SCHEMA_VERSION = 1
CONFIG_NAME_PATTERN = re.compile(r"[^A-Za-z0-9_.-]+")
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


def _safe_config_id(name):
    cleaned = CONFIG_NAME_PATTERN.sub("_", str(name).strip())
    cleaned = cleaned.strip("._-")[:64]
    if not cleaned:
        raise ValueError("Bitte einen Namen für die Konfiguration eingeben.")
    return cleaned


def _config_path(config_id):
    safe_id = _safe_config_id(config_id)
    return CONFIG_DIR / f"{safe_id}.json"


def _now_iso():
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _coerce_saved_config(payload):
    config = _coerce_config(payload)
    if "num_points" in payload:
        num_points = int(payload["num_points"])
        config["num_points"] = max(30, min(1200, num_points))
    return config


def _read_config_file(path):
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict) or "config" not in data:
        raise ValueError(f"Ungültige Konfigurationsdatei: {path.name}")
    return data


def _list_configs():
    CONFIG_DIR.mkdir(exist_ok=True)
    configs = []
    for path in sorted(CONFIG_DIR.glob("*.json")):
        try:
            data = _read_config_file(path)
        except Exception as exc:
            configs.append({"id": path.stem, "name": path.stem, "error": str(exc)})
            continue
        configs.append({
            "id": path.stem,
            "name": data.get("name", path.stem),
            "description": data.get("description", ""),
            "updated_at": data.get("updated_at", ""),
            "created_at": data.get("created_at", ""),
            "dashboard": data.get("config", {}).get("chamber_source", "manual"),
        })
    return configs


def _save_config(payload):
    name = str(payload.get("name", "")).strip()
    description = str(payload.get("description", "")).strip()
    raw_config = payload.get("config", {})
    if not isinstance(raw_config, dict):
        raise ValueError("Konfiguration muss ein JSON-Objekt sein.")

    config_id = _safe_config_id(name)
    CONFIG_DIR.mkdir(exist_ok=True)
    path = _config_path(config_id)
    created_at = _now_iso()
    if path.exists():
        try:
            created_at = _read_config_file(path).get("created_at", created_at)
        except Exception:
            pass

    document = {
        "schema_version": CONFIG_SCHEMA_VERSION,
        "name": name,
        "description": description,
        "created_at": created_at,
        "updated_at": _now_iso(),
        "config": _coerce_saved_config(raw_config),
    }
    tmp_path = path.with_suffix(".json.tmp")
    tmp_path.write_text(json.dumps(document, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    tmp_path.replace(path)
    return {"id": config_id, **document}


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

        if path == "/api/configs":
            _send_json(self, {"configs": _list_configs()})
            return

        if path.startswith("/api/configs/"):
            config_id = unquote(path.removeprefix("/api/configs/"))
            config_path = _config_path(config_id)
            if not config_path.exists():
                _send_json(self, {"error": "Konfiguration nicht gefunden"}, status=404)
                return
            try:
                data = _read_config_file(config_path)
                _send_json(self, {"id": config_path.stem, **data})
            except Exception as exc:
                _send_json(self, {"error": str(exc)}, status=400)
            return

        _send_json(self, {"error": "Nicht gefunden"}, status=404)

    def do_POST(self):
        path = urlparse(self.path).path
        if path == "/api/configs":
            try:
                length = int(self.headers.get("Content-Length", "0"))
                raw_body = self.rfile.read(length).decode("utf-8") if length else "{}"
                payload = json.loads(raw_body)
                saved = _save_config(payload)
                _send_json(self, saved)
            except Exception as exc:
                _send_json(self, {"error": str(exc)}, status=400)
            return

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
