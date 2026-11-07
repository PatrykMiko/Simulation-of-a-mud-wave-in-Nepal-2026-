# Jak zainstalować i uruchomić model (krok po kroku)

Instrukcja jest napisana tak, żeby dało się wszystko zrobić bez niczyjej pomocy. Wszystkie polecenia
sprawdzono na Ubuntu 24.04 z GRASS GIS 8.3.2 i r.avaflow 4.0G Revision 7. Polecenia wpisujesz w zwykłym
terminalu, chyba że napisano inaczej.

## 0. Co jest potrzebne

- **Linux** (r.avaflow 4.0G działa na systemach Linux i opiera się na GRASS GIS). Na Windowsie użyj WSL2, maszyny
  wirtualnej albo Dockera.
- GRASS GIS (nasz test: 8.3.2 z pakietów Ubuntu; autorzy r.avaflow dokumentują wersję 8.4.1, ale na 8.3.2 też
  działa), R z kilkoma pakietami (wbudowane wykresy r.avaflow), Python 3 z pakietami z `requirements.txt`.
- Około 2–3 GB wolnego miejsca (głównie pakiety R i dane przykładowe).

## 1. Instalacja

### 1a. Pakiety systemowe i Python

```bash
sudo apt update
sudo apt install grass grass-dev grass-doc gdal-bin build-essential \
                 python3-venv python3-pip r-base r-base-dev
cd <katalog-projektu>
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

### 1b. Pakiety R (kompilacja `terra` potrafi trwać kilkanaście minut)

```bash
sudo Rscript -e 'install.packages(c("Rcpp","terra","lattice","sp","raster","ROCR"), repos="https://cloud.r-project.org")'
Rscript -e 'for (p in c("Rcpp","terra","lattice","sp","raster","ROCR")) cat(p, requireNamespace(p, quietly=TRUE), "\n")'
```

Ostatnia komenda ma wypisać `TRUE` przy każdym pakiecie.

### 1c. Kod r.avaflow 4.0G

Pobierz źródła r.avaflow 4.0G ze strony projektu (`https://www.landslidemodels.org/r.avaflow/`, strona
r.avaflow.direct z instrukcją i linkiem do pobrania) i rozpakuj do katalogu `r.avaflow.40G/` w głównym
katalogu projektu. Wewnątrz mają leżeć pliki źródłowe, a nie zagnieżdżony drugi folder `r.avaflow.40G/`.
Ten katalog nie jest w Gicie (patrz `.gitignore`), więc każdy pobiera go sam. Zapisz wersję: 4.0G Revision 7
(2 września 2025).

### 1d. Lokalizacja GRASS i instalacja rozszerzenia

Jeśli używasz condy, wyłącz ją na czas pracy z GRASS-em (biblioteki z condy mogą kolidować z systemowym GRASS-em):

```bash
conda deactivate
mkdir -p ~/grassdata
grass -c XY ~/grassdata/xy_loc
```

Komunikat o braku `wxPython` i ostrzeżenie `SyntaxWarning` są nieszkodliwe (nie potrzebujemy graficznego
interfejsu). Wyjdź z sesji poleceniem `exit`, a potem:

```bash
grass --text ~/grassdata/xy_loc/PERMANENT
```

Prompt zmieni się na `GRASS xy_loc/PERMANENT:~ >`. Tam wpisz (podstaw swoją ścieżkę):

```bash
g.extension extension=r.avaflow.40G url=/ścieżka/do/projektu/r.avaflow.40G/
g.extension -a
exit
```

Lista z `g.extension -a` powinna zawierać `r.avaflow.40G`, `r.avaflow.main`, `r.avaflow.mult`,
`r.avaflow.background`, `r.avaflow.paraview`, `r.lakefill`, `r.scarp`.

## 2. Testy bez GRASS-a

```bash
python -m pytest -q
```

Powinno przejść 11 testów (metryki, walidator, generator, scenariusze, parser wyników). Nie wymagają GRASS-a.

## 3. Pierwszy przebieg: przykład Cascadia (test całego łańcucha)

```bash
mkdir -p data/cascadia && cd data/cascadia
wget https://www.landslidemodels.org/r.avaflow/data/Cascadia_r.avaflow40G.zip
unzip Cascadia_r.avaflow40G.zip
cd ../..

python -m ravaflow_tools.runner configs/ca_exp2.yaml     # około 3 minut
python -m ravaflow_tools.score  configs/ca_exp2.yaml     # kilka sekund
```

Używaj `python` z naszego środowiska wirtualnego (`venv/bin/python`) i bez aktywnej condy.

**Oczekiwany wynik** (referencja z naszych uruchomień): `runner` wypisze czas modelu około 100 s, całkowity około
3 minut (resztę zajmuje rysowanie wykresów) i CSI r.avaflow równe **0,549**. `score` wypisze własny indeks
Jaccarda około **0,556** (różnica wynika z przepróbowania referencji z 5 m na siatkę modelu 20 m), stosunek pól
około 1,51. Jeśli wyniki są wyraźnie inne, coś jest skonfigurowane inaczej.

## 4. Co dostajemy z przebiegu

Dla przebiegu o prefiksie `ca_exp2` w katalogu roboczym (`workdir` z konfiguracji):

| Ścieżka | Zawartość |
|---|---|
| `ca_exp2.run.sh`, `ca_exp2.run.log` | wygenerowany skrypt startowy i pełny log GRASS-a |
| `ca_exp2_results/ca_exp2_ascii/` | rastry ESRI ASCII: `*_elev.asc`, `*_hflowNNNN.asc` (grubość przepływu w kolejnych krokach czasu), `*_hflow_max.asc` (maksimum w czasie), analogicznie `tflow`, `pflow` (frakcje), `treach` (czas dojścia czoła) |
| `ca_exp2_results/ca_exp2_files/` | pliki tekstowe, m.in. `*_time.txt` (czas obliczeń w s) i `*_evaluation.txt` (ocena względem `impactarea`: CSI, FoC, długości spływu itd.) |
| `ca_exp2_results/ca_exp2_plots/` | wykresy i animacje generowane przez r.avaflow (R) |
| `ca_exp2_results/ca_exp2_vr/` | katalog dodatkowy r.avaflow (nie analizowaliśmy jego zawartości) |
| `ca_exp2_score/` | **nasze** wyniki oceny: `metrics.json`, `extent_comparison.png` (zgodne / nadmiar / niedobór), `hmax.png` |

Zasięg modelu to komórki, w których `*_hflow_max` ≥ 0,1 m (pole `threshold` w konfiguracji).

## 5. Własny przebieg: plik konfiguracji

Konfiguracja to plik YAML (wzór: `configs/ca_exp2.yaml`):

| Klucz | Znaczenie |
|---|---|
| `workdir` | katalog roboczy; ścieżki w `inputs` są względem niego, tu powstają wyniki |
| `grass_mapset` | `~/grassdata/xy_loc/PERMANENT` (lokalizacja z kroku 1d) |
| `region_from` | nazwa rastra GRASS, z którego ustawiamy region obliczeń |
| `inputs` | `nazwa_rastra_w_GRASS: plik.tif` (importowane przez `r.in.gdal -o`) |
| `pre_commands` | opcjonalnie: polecenia GRASS przed modelem (np. `r.lakefill`, `r.mapcalc`) |
| `flags` | opcjonalne flagi `r.avaflow.40G` (np. `-a`, `-m`) |
| `params` | parametry modelu `nazwa: wartość`; listy YAML zamieniane na `a,b,c`; **`prefix` obowiązkowy** |
| `reference_mask` | plik z zasięgiem rzeczywistym do oceny (`score`) |
| `threshold` | próg grubości [m] zasięgu modelu (domyślnie 0,1; **musi być taki sam dla obu modeli**) |
| `timeout_seconds` | opcjonalnie: przerwij przebieg po tylu sekundach |

Nazwy i znaczenie parametrów modelu (`friction`, `density`, `hrelease`, `time` itd.) pochodzą z r.avaflow, a
nie z naszego kodu. Pełny opis: manual na stronie projektu; krótka pomoc: w sesji GRASS `r.avaflow.40G --help`.
Dobrym wzorcem są skrypty z paczki Cascadia (`data/cascadia/ca.avaflow.start.sh`, 9 przykładowych eksperymentów).

## 6. Dane syntetyczne i dane z Nepalu

```bash
# dane syntetyczne w formacie kontraktu danych (GeoTIFF)
python -m ravaflow_tools.synthetic data/synthetic
python -m ravaflow_tools.check_inputs data/synthetic --volume 2e6 --cellsize 30 --epsg 32645
python -m ravaflow_tools.runner configs/synthetic_test.yaml
python -m ravaflow_tools.score  configs/synthetic_test.yaml
```

Zasięg „rzeczywisty” w danych syntetycznych jest wymyślony, więc Jaccard z tego przebiegu **nic nie mówi o jakości
modelu**; test sprawdza tylko, że cały łańcuch przechodzi na danych w tym formacie.

Dla prawdziwych danych: rozpakuj je do `data/nepal/`, uruchom walidator (`data/README.md`), skopiuj
`configs/synthetic_test.yaml` do `configs/langtang_base.yaml`, zmień nazwy plików, `workdir`, `cellsize`, czas
symulacji i parametry, a potem uruchom `runner` i `score`.

## 7. Scenariusze (S1–S3) i serie przebiegów

Plik `configs/scenarios_synthetic.yaml` pokazuje format: konfiguracja bazowa plus lista wariantów.

```bash
python -m ravaflow_tools.scenarios configs/scenarios_synthetic.yaml          # tylko przygotuj rastry i konfiguracje
python -m ravaflow_tools.scenarios configs/scenarios_synthetic.yaml --run    # przygotuj i uruchom po kolei
```

Warianty: `scale_release` (mnożnik grubości masy, np. 0,25, 0,5, 2), `barrier` (tama: odcinek we współrzędnych
mapy, wysokość i szerokość w metrach, dodawane do NMT) i `params` (nadpisanie dowolnych parametrów, np. udziału
faz). Po uruchomieniu powstaje `scenarios_summary.csv` (Jaccard, czas itd. dla każdego wariantu). Wszystkie
scenariusze są eksploracyjne: parametry pochodzą z kalibracji na zdarzeniu bazowym, więc wyniki porównujemy
głównie między modelami i względem bazy.

## 8. Ważne zasady i pułapki

- **Jeden mapset GRASS to jeden przebieg naraz.** Przebiegi uruchamiaj po kolei. Równoległe wymagają osobnych
  lokalizacji (jeszcze nie zaimplementowane).
- **Nie naciskaj `Ctrl+Z` podczas przebiegu** (zatrzymuje proces, który potem blokuje mapset). Do przerwania służy
  `Ctrl+C`.
- Rysowanie wykresów przez r.avaflow zajmuje ponad połowę czasu przebiegu. Do pomiaru kosztu obliczeń używaj
  czasu modelu (`*_time.txt`, pole `model_seconds`), nie czasu całkowitego.
- Czas mierz zawsze na tym samym sprzęcie i zapisuj liczbę rdzeni (wymóg porównania kosztów w RQ5).
- Raster `elevation` to powierzchnia **podstawy** masy, a nie teren sprzed zdarzenia (zob. `data/README.md`).
- Rastry `release_height` poza obszarem odrywania mają mieć dokładnie 0, nie wartość brakującą.

## 9. Typowe problemy

| Objaw | Przyczyna i rozwiązanie |
|---|---|
| `g.extension: nie znaleziono polecenia` | polecenie wpisane poza sesją GRASS; użyj `grass --text ~/grassdata/xy_loc/PERMANENT` |
| `ERROR: ... currently running GRASS in selected mapset ... .gislock` | inny przebieg trwa albo został przerwany; `pgrep -a grass`; jeśli pusto: `rm ~/grassdata/xy_loc/PERMANENT/.gislock` |
| `runner` nic nie wypisuje przez kilka minut | to normalne, przebieg trwa; wypisuje się dopiero po zakończeniu (log na bieżąco w `*.run.log`) |
| `top` pokazuje stan `T` procesu `grass` | proces zatrzymany; `pkill -9 -P <PID>; kill -9 <PID>`, potem usuń `.gislock` |
| `ERROR: wxGUI requires wxPython` | nieszkodliwe, GUI nie jest potrzebne |
| `[uwaga] siatka referencji (5.0 m) != siatka modelu (20.0 m)` | informacja, nie błąd; referencja jest przepróbowana najbliższym sąsiadem |
| wykresy się nie generują | sprawdź pakiety R (krok 1b) |
