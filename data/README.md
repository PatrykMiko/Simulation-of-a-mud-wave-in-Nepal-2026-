# Katalog `data/`

Zawartość tego katalogu **nie trafia do Gita** (patrz `.gitignore`): rastry są duże i mają własne
licencje. Do repozytorium trafia tylko ten plik.

## Układ

| Podkatalog | Co zawiera | Skąd |
|---|---|---|
| `data/cascadia/` | przykład szkoleniowy r.avaflow 4.0G (dane i skrypt startowy) | `https://www.landslidemodels.org/r.avaflow/data/Cascadia_r.avaflow40G.zip` |
| `data/synthetic/` | dane syntetyczne w formacie kontraktu danych | `python -m ravaflow_tools.synthetic data/synthetic` |
| `data/nepal/` | **prawdziwe dane (Langtang 2026)**, dostarczane przez osoby od danych | do uzupełnienia |

## Kontrakt danych dla `data/nepal/`

Wszystkie rastry w formacie GeoTIFF, z **identyczną siatką**: ten sam zasięg, ten sam rozmiar komórki
(plan: 30 m) i ten sam układ rzutowany w metrach (propozycja: UTM 45N, EPSG:32645).

| Plik | Znaczenie | Wymagany |
|---|---|---|
| `elevation.tif` | powierzchnia **podstawy** masy [m n.p.m.], czyli teren minus masa osuwiska; `elevation + release_height` = teren sprzed zdarzenia | tak |
| `release_height.tif` | grubość masy w obszarze odrywania [m]; poza nim **dokładnie 0** (nie nodata) | tak |
| `reference_extent.tif` | zasięg rzeczywisty: 1 = zasięg, 0 = brak zasięgu, wartość ujemna lub nodata = brak informacji | tak |
| `entrainment_height.tif` | grubość materiału do unoszenia [m]; poza obszarem 0 | nie |

Do każdego pliku krótka notka: źródło, data, sposób powstania, założenia (np. `data/nepal/NOTES.md`).

## Walidacja przed użyciem

```bash
python -m ravaflow_tools.check_inputs data/nepal --volume 200e6 --cellsize 30 --epsg 32645
```

Skrypt wypisuje listę błędów i ostrzeżeń; kod wyjścia 1 oznacza, że dane nie spełniają kontraktu.
`--volume` to objętość, którą założyliście przy tworzeniu `release_height.tif` (w m³).

