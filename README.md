# Symulacja fali błotnej w Nepalu (2026): SCIDDICA-SS3 vs r.avaflow v4

Odtwarzamy numerycznie zdarzenie z 26 sierpnia 2026
(zawalenie skał i lodowca w Parku Narodowym Langtang, spływ gruzowy i powódź wzdłuż rzek Lende Khola i
Trishuli) dwoma modelami o skrajnie różnej złożoności: automatem komórkowym **SCIDDICA-SS3** i wielofazowym
modelem ciągłym **r.avaflow v4**. Porównujemy je pod względem dokładności zasięgu, wrażliwości na dane i
kosztu obliczeń, a potem sprawdzamy scenariusze „co by było, gdyby” (tama na trasie, inna masa lodowca, inny
udział lodu i wody).



## Gdzie co jest

```
.
├── README.md                  ten plik
├── requirements.txt           pakiety Pythona
├── configs/                   konfiguracje przebiegów (YAML)
│   ├── ca_exp2.yaml             test na przykładzie Cascadia (eksperyment 2, ze znanym zasięgiem)
│   ├── synthetic_test.yaml      test na danych syntetycznych w formacie kontraktu danych
│   └── scenarios_synthetic.yaml przykład serii scenariuszy S1–S3
├── ravaflow_tools/            nasz kod
│   ├── runner.py                konfiguracja -> skrypt startowy -> GRASS -> r.avaflow, pomiar czasu
│   ├── score.py                 ocena przebiegu względem zasięgu rzeczywistego, mapy i metryki
│   ├── metrics.py               metryki zgodności zasięgu (niezależne od modelu), wczytywanie rastrów
│   ├── viz.py                   mapy porównania zasięgu, maksymalna grubość, animacje
│   ├── check_inputs.py          walidator danych wejściowych (kontrakt danych)
│   ├── synthetic.py             generator danych syntetycznych
│   └── scenarios.py             scenariusze (masa, tama, parametry) z pliku YAML
├── tests/test_all.py          testy automatyczne (nie wymagają GRASS-a)
├── docs/
│   ├── uruchamianie.md          jak zainstalować i uruchomić (krok po kroku)
│   ├── model_naukowo.md         jak działa model, jak liczymy metryki, ograniczenia
│   └── stan_prac.md             co zrobiono, jakie eksperymenty, co oznaczają wyniki
├── data/                      dane (poza Gitem; zob. data/README.md)
└── r.avaflow.40G/             kod r.avaflow (poza Gitem; pobierany osobno)
```

## Szybki start

1. Zainstaluj środowisko (GRASS GIS, r.avaflow 4.0G, pakiety Pythona): `docs/uruchamianie.md`, część 1.
2. Sprawdź, że wszystko działa: `python -m pytest -q` (testy bez GRASS-a), a potem przebieg na przykładzie
   Cascadia: `docs/uruchamianie.md`, część 3.
3. Dla własnych danych: `python -m ravaflow_tools.check_inputs data/nepal --volume ... --cellsize 30`, potem
   konfiguracja w `configs/` i `python -m ravaflow_tools.runner configs/<plik>.yaml`.

## Protokół porównania obu modeli

Wyniki mają być porównywalne, więc oba modele korzystają z tego samego:

- **danych wejściowych**: to samo NMT (siatka 30 m), ten sam obszar odrywania i objętość, ten sam zasięg
  rzeczywisty (kontrakt danych: `data/README.md`),
- **definicji zasięgu modelu**: komórka należy do zasięgu, jeśli maksymalna (w czasie) grubość przepływu jest
  równa lub większa od progu **0,1 m**,
- **metryk**: indeks Jaccarda e = |R ∩ S| / |R ∪ S| (iloraz, bez pierwiastka), precision, recall, F1, stosunek
  pól, odległość Hausdorffa. Kod w `ravaflow_tools/metrics.py`.

## Dla osoby od SS3: jak ocenić wynik swojego modelu naszymi metrykami

Potrzebujesz rastra maksymalnej (w czasie) grubości złoża z SS3, zapisanego w tej samej siatce co dane
(GeoTIFF albo ESRI ASCII). Jeśli siatki się różnią, funkcja `sample_to_grid` przepróbuje referencję.

```python
from ravaflow_tools import metrics as M

hmax, meta = M.read_raster("wynik_ss3_maks_grubosc.asc")
ref, ref_meta = M.read_raster("data/nepal/reference_extent.tif")
if not M.same_grid(meta, ref_meta):
    print(M.grid_differences(meta, ref_meta))   # co dokładnie się różni
    ref = M.sample_to_grid(ref, ref_meta, meta)

obs, valid = M.reference_masks(ref, {**meta, "nodata": ref_meta["nodata"]})
model = M.model_extent(hmax, meta, threshold=0.1)
print(M.score(model, obs, valid, meta["cellsize"]))   # jaccard, precision, recall, f1, ...
```
Najprościej zapisać wynik SS3 dokładnie na siatce reference_extent.tif.
## Dla osób od danych

Format, jakiego potrzebuje model: `data/README.md`. Zanim wyślesz pliki, sprawdź je walidatorem:

```bash
python -m ravaflow_tools.check_inputs data/nepal --volume 200e6 --cellsize 30 --epsg 32645
```

Najczęstsze błędy: raster `release_height` z wartością brakującą zamiast zera poza obszarem odrywania, inna
siatka lub rozmiar komórki w poszczególnych plikach, układ geograficzny (stopnie) zamiast rzutowanego (metry),
`elevation` będący powierzchnią sprzed zdarzenia zamiast powierzchni podstawy masy.

## Stan prac

Krótko: pipeline r.avaflow jest uruchomiony i sprawdzony na przykładzie Cascadia, a walidator, dane syntetyczne
i scenariusze są gotowe i przetestowane bez GRASS-a. Szczegóły, wyniki i co z nich wynika: `docs/stan_prac.md`.


## Źródła i licencje

- r.avaflow 4.0G (Mergili i in. 2025, GMD 18: 9879–9896) jest oprogramowaniem na licencji GPL; pobieramy je
  osobno ze strony projektu (`https://www.landslidemodels.org/r.avaflow/`) i nie wersjonujemy w tym repozytorium.
- Dane Cascadia pochodzą z paczki przykładowej r.avaflow; dane z Nepalu mają źródła opisane w notkach do danych.

