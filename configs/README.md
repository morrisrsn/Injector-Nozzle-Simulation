# Gespeicherte Konfigurationen

Dieser Ordner enthaelt JSON-Konfigurationen, die ueber das Web-Frontend gespeichert werden.

- Jede Konfiguration liegt als einzelne `.json`-Datei vor.
- Dateinamen werden vom Server bereinigt, damit keine Pfade ausserhalb dieses Ordners beschrieben werden.
- Gespeichert werden nur bekannte Eingabeparameter aus `NozzleInput` plus `num_points`.
- Die Dateien koennen normal mit dem Repository versioniert und geteilt werden.
