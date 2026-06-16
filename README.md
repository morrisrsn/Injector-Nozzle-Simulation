# Düsenrechner

Ein lokaler Düsen- und Injektorrechner mit Python-Backend und HTML-Frontend. Das Projekt berechnet thermodynamische Zustände, Düsenströmung, Schub, Geometrie, Kontur- und Strömungsprofile für Kaltgas- und Verbrennungsanwendungen.

## Kurzbeschreibung

Dieses Repository stellt ein lokales Webtool zur Auslegung und Analyse von Düsen bereit. Über eine browserbasierte Oberfläche können Betriebsbedingungen, Medien, Geometrie- und Injektorparameter eingegeben werden. Die eigentliche Berechnung läuft in Python und liefert Kennwerte wie Kammerzustand, Austrittsmachzahl, Austrittsgeschwindigkeit, Schub, Düsenabmessungen, Konturpunkte und Strömungsprofile.

## Hauptfunktionen

- Lokales HTML-Frontend im Browser
- Python-Backend über einen einfachen HTTP-Server
- Berechnung von Kaltgas- und Verbrennungsfällen
- Auslegung einer Düse aus gegebenem Massenstrom
- Auswertung einer vorgegebenen Düsengeometrie
- Berechnung von:
  - Kammerdruck und Kammerzustand
  - Temperatur, Dichte, spezifischer Gaskonstante und Gamma
  - Hals-, Austritts- und Kammergeometrie
  - Austritts-Machzahl
  - Austrittsdruck
  - Austrittstemperatur
  - Austrittsgeschwindigkeit
  - Druckschub und Gesamtschub
- Optionale Kopplung von Injektor und Düse
- Darstellung von:
  - Kennwerten
  - Thermodynamikdaten
  - Performance-Daten
  - Geometrie
  - Injektordaten
  - Düsenkontur
  - Strömungsprofilen
- Export der Konturdaten als CSV

## Projektstruktur

```text
.
├── web_app.py
├── Thermodynamics.py
└── frontend
    └── index.html
```

### `web_app.py`

`web_app.py` ist der lokale Server des Projekts. Die Datei startet einen HTTP-Server, liefert das Frontend aus und stellt zwei API-Endpunkte bereit.

Wichtige Aufgaben:

- Startet den lokalen Webserver
- Öffnet optional automatisch den Browser
- Liefert `frontend/index.html` aus
- Sendet Standardwerte an das Frontend
- Empfängt Eingaben vom Frontend
- Ruft die Berechnungsfunktion `calculate_nozzle(...)` auf
- Gibt die Ergebnisse als JSON an den Browser zurück

Wichtige Routen:

```text
GET  /
GET  /index.html
GET  /api/defaults
POST /api/calculate
```

### `Thermodynamics.py`

`Thermodynamics.py` enthält die eigentliche Rechenlogik. Dort sind die Eingabeparameter, thermodynamischen Modelle, Geometrieberechnung, Injektorberechnung und Düsenströmung implementiert.

Wichtige Bestandteile:

- `NozzleInput`: zentrale Dataclass für alle Eingabeparameter
- `GasState`: thermodynamischer Zustand des Gases
- `NozzleGeometry`: geometrische Beschreibung der Düse
- Thermodynamikfunktionen für Kaltgas und Verbrennung
- vereinfachte Gleichgewichts- bzw. vollständige Verbrennungsrechnung
- Injektormodell für Kraftstoff- und Oxidatorseite
- Düsenströmungsmodell
- Konturgenerierung
- Profilberechnung entlang der Düse
- zentrale Funktion `calculate_nozzle(...)`

Die Funktion `calculate_nozzle(...)` ist der wichtigste Einstiegspunkt für die Berechnung. Sie nimmt ein Konfigurations-Dictionary entgegen und gibt ein Ergebnis-Dictionary zurück.

Typische Ergebnisbereiche:

```text
inputs
state
performance
geometry
contour
flowfield
injector
files
warnings
```

### `frontend/index.html`

`index.html` ist die Benutzeroberfläche. Sie enthält HTML, CSS und JavaScript in einer Datei.

Wichtige Aufgaben:

- Aufbau der Eingabemaske
- Umschaltung zwischen Dashboards
- Umschaltung zwischen Auslegungsmodus und Geometriemodus
- Senden der Eingaben an `/api/calculate`
- Anzeigen der Ergebnisdaten
- Tabellen für Thermodynamik, Performance, Geometrie und Injektor
- Visualisierung der Düsenkontur
- Visualisierung von Profilen
- CSV-Download

## Installation

### 1. Repository klonen

```bash
git clone <repository-url>
cd <repository-name>
```

### 2. Virtuelle Umgebung erstellen

Windows:

```bash
python -m venv .venv
.venv\Scripts\activate
```

macOS/Linux:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

### 3. Abhängigkeiten installieren

Mindestens benötigt:

```bash
pip install numpy pandas matplotlib scipy
```

Optional:

```bash
pip install CoolProp cantera
```

Die optionalen Bibliotheken verbessern die Stoffdaten- bzw. Verbrennungsmodellierung. Wenn sie nicht installiert sind, nutzt das Programm vereinfachte Fallback-Modelle.

## Programm starten

```bash
python web_app.py
```

Standardmäßig läuft das Frontend danach unter:

```text
http://127.0.0.1:8000
```

Falls der Browser nicht automatisch geöffnet werden soll:

```bash
python web_app.py --no-browser
```

Mit anderem Port:

```bash
python web_app.py --port 8080
```

Mit anderem Host:

```bash
python web_app.py --host 0.0.0.0 --port 8000
```

## Bedienung

### Dashboard „Düse“

Dieses Dashboard ist für klassische Düsenrechnungen mit manuell vorgegebenen Kammerbedingungen gedacht.

Typische Eingaben:

- Betriebsmodus: Verbrennung oder Kaltgas
- Kammerdruck
- Umgebungsdruck
- Eintrittstemperaturen
- Massenstrom
- Düsenwirkungsgrad
- Entladungskoeffizient
- Brennstoff und Oxidator
- Äquivalenzverhältnis
- Geometrie- bzw. Auslegungsparameter

### Dashboard „Injektor + Düse“

Dieses Dashboard koppelt ein vereinfachtes Injektormodell mit der Düsenrechnung.

Typische Eingaben:

- Kraftstofftankdruck
- Kraftstofftanktemperatur
- Kraftstofftankvolumen
- Kraftstoff-Injektordurchmesser
- Anzahl Kraftstoffbohrungen
- Kraftstoff-Injektor-Cd
- Oxidatortankdruck
- Oxidatortanktemperatur
- Oxidatortankvolumen
- Oxidator-Injektordurchmesser
- Anzahl Oxidatorbohrungen
- Oxidator-Injektor-Cd
- Düsengeometrie hinter dem Injektor

Das Modell bestimmt daraus unter anderem Massenströme, Mischungsverhältnis, Kammerdruck und anschließend die Düsenperformance.

## Rechenmodi

### Kaltgasmodus

Im Kaltgasmodus wird ein Gas ohne Verbrennungsrechnung betrachtet. Der Zustand wird aus den vorgegebenen Druck- und Temperaturwerten sowie Stoffdaten berechnet.

Geeignet für:

- einfache Druckgasdüsen
- Sauerstoff, Luft oder andere kalte Gase
- schnelle Abschätzungen ohne Reaktion

### Verbrennungsmodus

Im Verbrennungsmodus wird ein Reaktionsgas aus Brennstoff und Oxidator berechnet. Unterstützt werden aktuell Ethanol, Methan, n-Butan/C4H10 und Wasserstoff.

Unterstützte bzw. vorgesehene Oxidatoren:

- Luft
- Sauerstoff
- Lachgas/N2O

Es gibt zwei Modellvarianten:

- `equilibrium`: vereinfachte Gibbs-Gleichgewichtsrechnung mit NASA-7-Daten
- `complete`: vollständige Verbrennung ohne detaillierte Dissoziation

## Auslegungsmodus und Geometriemodus

### Düse auslegen

Im Auslegungsmodus wird die Düse aus dem gewünschten Massenstrom und dem Kammerzustand berechnet. Dabei werden Halsfläche, Austrittsfläche und weitere Geometriegrößen bestimmt.

Typische Ausgabegrößen:

- Halsdurchmesser
- Austrittsdurchmesser
- Kammerdurchmesser
- konvergente Länge
- divergente Länge
- Gesamtlänge der Düse

### Austrittsbedingungen aus Vorgabe

Im Geometriemodus wird eine vorgegebene Düse ausgewertet. Der Nutzer gibt Halsdurchmesser, Austrittsdurchmesser, Kammerdurchmesser und Kammerlänge vor. Das Programm berechnet daraus die Strömungs- und Performancewerte.

## Ergebnisbereiche

### Kennwerte

Die wichtigsten Ergebnisse werden direkt oben im Dashboard angezeigt:

- Schub
- Austritts-Machzahl
- Halsdurchmesser
- Austrittsdurchmesser

### Thermodynamik

Dieser Bereich zeigt den berechneten Gaszustand:

- Kammerdruck
- Kammer-/Gastemperatur
- Dichte
- spezifische Gaskonstante
- molare Masse
- Isentropenexponent Gamma
- Wärmekapazität
- Schallgeschwindigkeit
- Gaszusammensetzung

### Performance

Dieser Bereich zeigt die düsenbezogenen Leistungswerte:

- Choking-Zustand
- effektiver Massenstrom
- Flächenverhältnis
- Austritts-Machzahl
- Austrittsdruck
- ideale und reale Austrittsgeschwindigkeit
- Austrittstemperatur
- Druckschub
- Gesamtschub

### Geometrie

Dieser Bereich zeigt die berechneten oder vorgegebenen geometrischen Größen:

- Kammerfläche
- Halsfläche
- Austrittsfläche
- Kammerdurchmesser
- Halsdurchmesser
- Austrittsdurchmesser
- Kammerlänge
- konvergente Düsenlänge
- divergente Düsenlänge
- Gesamt-Düsenlänge
- Konvergenzwinkel
- Divergenzwinkel

### Injektor

Dieser Bereich zeigt die Ergebnisse der gekoppelten Injektor-Düsen-Rechnung, falls das Injektor-Dashboard aktiv ist.

Typische Werte:

- berechneter Kammerdruck
- Gesamtmassenstrom
- Düsenkapazität
- reales O/F-Verhältnis
- reales Äquivalenzverhältnis
- Kraftstoffmassenstrom
- Oxidatormassenstrom
- Austrittsgeschwindigkeiten der Injektorströme
- Phaseninformationen
- grobe Tankmassen
- geschätzte Brenndauer

### Kontur

Die Konturtabelle enthält die Koordinaten der Düsenkontur entlang der x-Achse. Typische Größen:

- axiale Position `x`
- Radius `r`

### Profile

Die Profiltabelle enthält Strömungsgrößen entlang der Düse:

- axiale Position `x`
- Machzahl
- Druck
- Temperatur
- Geschwindigkeit

## Visualisierung

Das Frontend stellt verschiedene Plotmodi bereit:

- 2D-Kontur
- alle Verläufe
- 3D-Ansicht
- Schnittansicht

Zusätzlich kann ein Feldoverlay gewählt werden, zum Beispiel:

- Machzahl
- Druck
- Temperatur
- Geschwindigkeit

## CSV-Export

Über den CSV-Button können die berechneten Konturpunkte exportiert werden. Der Export enthält aktuell die Koordinaten:

```text
x [m], r [m]
```

## Hinweise zu Dateinamen

In `web_app.py` wird aktuell importiert mit:

```python
from Thermodynamics import NozzleInput, calculate_nozzle
```

Daher muss die Rechendatei exakt `Thermodynamics.py` heißen. Auf manchen Systemen ist die Groß-/Kleinschreibung wichtig. Falls die Datei `thermodynamics.py` heißt, sollte entweder die Datei umbenannt oder der Import angepasst werden:

```python
from thermodynamics import NozzleInput, calculate_nozzle
```

## Bekannte Einschränkungen

- Die Verbrennungsrechnung ist vereinfacht und nutzt reduzierte Produktchemie mit groben Stoffdaten-Fallbacks.
- Ohne CoolProp werden Stoffdaten über vereinfachte Fallback-Modelle berechnet.
- Ohne Cantera wird keine detaillierte externe Verbrennungschemie genutzt.
- Die Injektorberechnung ist ein vereinfachtes 0D-Modell.
- N2O-Zweiphasenströmung wird nur grob angenähert.
- Die Warnung zur Enthalpie-Nullstelle kann auftreten, wenn im aktuellen Suchbereich keine saubere Vorzeichenänderung der Enthalpiebilanz gefunden wird.
- Die Ergebnisse sind als Engineering-Abschätzung zu verstehen und nicht als validierte Auslegung für sicherheitskritische Hardware.

## Typischer Workflow

1. Web-App starten.
2. Dashboard auswählen:
   - `Düse` für manuelle Kammerbedingungen
   - `Injektor + Düse` für gekoppelte Injektor-Düsen-Rechnung
3. Rechenmodus wählen:
   - `Düse auslegen`
   - `Austrittsbedingungen`
4. Betriebsmodus wählen:
   - `Verbrennung`
   - `Kaltgas`
5. Parameter eintragen.
6. Berechnung starten.
7. Kennwerte, Tabellen und Plots prüfen.
8. Kontur bei Bedarf als CSV exportieren.

## Beispielstart

```bash
python web_app.py
```

Dann im Browser öffnen:

```text
http://127.0.0.1:8000
```

## Beispielhafte Eingaben

### Verbrennung mit Ethanol und Luft

```text
mode = combustion
fuel = C2H5OH
oxidizer = Air
phi = 1.0
p_c = 2e6 Pa
p_amb = 101325 Pa
mdot = 0.1 kg/s
eta_nozzle = 0.95
```

### Kaltgas mit Sauerstoff

```text
mode = cold_gas
gas_medium = O2
p_c = 3e6 Pa
T_in = 300 K
p_amb = 101325 Pa
mdot = 0.1 kg/s
```
