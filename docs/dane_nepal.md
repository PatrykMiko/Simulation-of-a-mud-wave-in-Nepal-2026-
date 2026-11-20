# Dane wejściowe Nepal: opis dla osób od modeli

Dotyczy wersji danych **v2** (`data/nepal/derived/v2/`). Dane obejmują odcinek od krótkiej blokady
rzeki Lhende Khola przy zbiegu z doliną lodowcową w dół doliny, do krawędzi obszaru roboczego na zachodzie.
Etap od oderwania lodowca do blokady NIE jest w tych danych.

Źródła: DEM-y Copernicus GLO-30 (COP30), JAXA ALOS AW3D30 i SRTM GL1, przeliczone interpolacją dwuliniową
na siatkę 30 m w UTM 45N. Oś rzeki, zasięg obserwowany i polygon źródła narysowano ręcznie w QGIS na zdjęciach
satelitarnych; bufor 2 km wokół osi rzeki wyznacza obszar roboczy.

## 1. Szybki start

1. Kontrola danych (z głównego folderu repo):
   ```
   python -m ravaflow_tools.verify_derived --version v2 --volume 140e6
   python -m ravaflow_tools.check_inputs data/nepal/derived/v2/copernicus --volume 140e6 --cellsize 30 --epsg 32645
   ```
   `check_inputs` powtórzcie dla `alos` i `srtm`. Oczekiwane: `błędów 0` (3 uwagi o grubości masy są znane,
   patrz punkt 6) i `WYNIK: OK`.
2. Każdy model dostaje **folder z czterema plikami**, tak jak dotąd: `data/nepal/derived/v2/<dem>`,
   gdzie `<dem>` to `copernicus`, `alos` albo `srtm`.
3. Uruchamiajcie **wszystkie trzy DEM-y** z tymi samymi parametrami modelu. Rozrzut wyników między
   DEM-ami to miara niepewności danych terenowych.

Sposób uruchamiania modeli znacie. Nie edytujcie plików `.tif` ręcznie.

## 2. Pliki (w każdym folderze DEM-u)

Wszystkie: GeoTIFF, EPSG:32645 (UTM 45N), komórka 30 m, **identyczna siatka** (1971 kolumn x 2096 wierszy), float32.

| Plik | Znaczenie | Wartości specjalne |
|---|---|---|
| `elevation.tif` | teren sprzed zdarzenia [m n.p.m.] | poza obszarem roboczym dokładnie 0 |
| `release_height.tif` | grubość masy w t = 0 [m], leży NA `elevation` | poza polygonem źródła dokładnie 0 |
| `entrainment_height.tif` | maks. grubość materiału do porwania [m]: 0,5 / 1,5 / 3,0 | poza obszarem roboczym 0 |
| `reference_extent.tif` | zasięg obserwowany ze zdjęć: 1 = zasięg, 0 = brak zasięgu | **-1 = brak informacji (poza obszarem roboczym)** |
| `manifest.json` | parametry tej wersji (objętość, poziom Z, hash polygonu) | |

Opis powstania plików i statystyki (liczba dziur, załatanych komórek itd.) jest w `data/nepal/derived/v2/NOTES.md`
(jeden plik na całą wersję).

Poza folderami DEM-ów: `data/nepal/domain_mask.tif` (1 = obszar roboczy), `data/nepal/<dem>/dem_original.tif`
(DEM po obróbce, -9999 poza obszarem; wejście do ponownego generowania).

**Uwaga:** folder `data/nepal/<dem>/` zawiera też starszy zestaw danych (źródło w bliźnie oderwania, inna siatka:
1969 x 2094). Do etapu od blokady używajcie wyłącznie `data/nepal/derived/v2/<dem>`.

## 3. Najważniejsze fakty o danych

- **`elevation` NIE ma odjętej masy.** Źródło leży w dnie doliny (masa nadeszła z góry), więc
  `release_height` jest nałożone na niezmieniony teren. To odstępstwo od opisu w kontrakcie danych
  (tam `elevation + release_height` = teren sprzed zdarzenia, co dotyczy źródła w bliźnie oderwania).
  W r.avaflow wysokość uwolnienia jest nakładana na model terenu, więc to właściwy układ.
- **Masa źródłowa (v2):** 140 mln m³ w polygonie 1,589 km² (1766 komórek), środek ok. x = 351 204, y = 3 135 060
  (UTM 45N). Składa się z warstwy bazowej 44 m w całym polygonie (69,9 mln m³) oraz poziomej tafli
  (70,1 mln m³) wypełniającej obniżenia terenu do poziomu Z (3071,8 m dla Copernicusa; 3072,1 ALOS;
  3076,3 SRTM). Maksymalna grubość 263-274 m (w najgłębszym miejscu wąwozu).
- **Podział na fazy (stała / płynna) nie jest zapisany w rastrach.** `release_height` to całkowita grubość.
  Udział fazy płynnej ustawiacie w konfiguracji modelu (parametr scenariusza). r.avaflow nie ma wbudowanego
  topnienia lodu (według dokumentacji, sprawdźcie w swojej wersji), więc udział płynu na starcie to decyzja
  scenariuszowa. Stosujcie ten sam udział w obu modelach.
- **Porywanie:** `entrainment_height` to górny limit grubości, nie wartość porywana. Cały obszar roboczy ma
  limit (łącznie ok. 585 mln m³ dostępnych). Wartość docelowa porwania z dokumentu kalibracyjnego to
  58,6 mln m³, więc limit nie jest tu ograniczeniem. Zapisujcie w wynikach, ile model faktycznie porwał.
- **Obszar roboczy** = bufor 2 km wokół osi rzeki przycięty do pokrycia wszystkich trzech DEM-ów
  (434 238 komórek). Fala wychodzi z obszaru na zachodniej krawędzi (595 komórek obserwowanego zasięgu
  leży już poza obszarem). Poza obszarem `elevation` = 0, czyli tworzy się „klif”: masa, która tam dojdzie,
  spada i znika. Warto liczyć, ile masy opuszcza obszar.

## 4. Uwagi do modeli

**r.avaflow**
- Wejście: `elevation`, `release_height` (u Was w jednym rastrze, podział na fazy w parametrach),
  `entrainment_height` (rastr maksymalnej porywanej grubości).
- Początkowa kolumna masy jest gruba (do ok. 270 m). Spodziewajcie się gwałtownego początku, więc
  zwróćcie uwagę na krok czasowy i stabilność w pierwszych sekundach.

**SCIDDICA SS3**
- Konwersja rastrów na siatkę heksagonalną: `ravaflow_tools/to_ss3.py`. Powtarzajcie ją po każdej nowej
  wersji danych. Siatka w modelu SS3 musi uwzględniać, że poza obszarem roboczym wysokość to 0 (klif).
- Sprawdźcie, że `release_height` jest w SS3 **dodawane** do `elevation` (a nie odejmowane ani traktowane jak
  wysokość bezwzględna).
- `reference_extent` zawiera -1: nie wczytujcie go jako „0 = brak zasięgu”.

## 5. Ocena wyników

- Zasięg: porównujcie z `reference_extent` **tylko w komórkach ≥ 0**. Komórki -1 wyłączcie z liczenia.
  Rozważcie wyłączenie pasa kilku km przy zachodniej krawędzi, gdzie obserwowany zasięg jest ucięty
  granicą obszaru.
- Markery kalibracyjne (czas dotarcia do Gyirong, prędkości, nabieganie): `data/calibration/dane_do_kalibracji.md`.
  Przed użyciem w pracy sprawdźcie źródła liczb w tym pliku.
- Kalibrujcie na części markerów (np. zasięg), a resztę zostawcie do walidacji. Dajcie obu modelom ten sam
  budżet strojenia.

## 6. Znane ograniczenia i założenia (opisać w pracy)

1. Objętość i kształt masy źródłowej to **założenie** (literatura podaje rozbieżne szacunki całkowitej objętości).
   Maks. grubość 263-274 m wynika głównie z kształtu wąwozu (wąski polygon w stromej dolinie). Skrypt zgłasza to jako UWAGA.
2. DEM-y są **sprzed zdarzenia** i mają 30 m, więc wąski wąwóz jest słabo rozdzielony.
3. **Łatanie artefaktów:** w ALOS załatano 1873 komórki, w SRTM 1014 komórek wartościami z Copernicusa
   (tam, gdzie DEM odbiegał od Copernicusa o ponad 200 m, a drugi DEM się z nim zgadzał). Na torze fali nie
   zostały komórki z różnicą ponad 200 m. W łatanych miejscach trzy DEM-y nie są niezależne.
4. Grubości porywania (0,5 / 1,5 / 3,0 m wg nachylenia) to założenie empiryczne.
5. Zasięg obserwowany (`reference_extent`) pochodzi z ręcznego obrysu na zdjęciach satelitarnych.

## 7. Generowanie innej wersji danych

Skrypty: `ravaflow_tools/prepare_nepal_data.py`, kontrola: `ravaflow_tools/verify_derived.py`.
Wejścia: `data/nepal/vector/nepal_vectors.gpkg` (warstwy `rzeka`, `bufor_rzeki`, `zasieg`),
`data/nepal/source/source_<wersja>.gpkg` (warstwa `zrodlo`), surowe DEM-y (ścieżki na górze skryptu).

Przebudowa DEM-ów od zera (bez `--reuse-dem`) wymaga surowych plików, a ścieżki do nich w `RAW` na górze skryptu
są lokalne (komputer autora). Pozostałym wystarczy `--reuse-dem`, bo `dem_original.tif` i `domain_mask.tif` są
w repozytorium.

Inna objętość lub grubość bazy (bez przebudowy DEM-ów), przykład:
```
python -m ravaflow_tools.prepare_nepal_data --version v2_110 --volume 110e6 --base-m 34.6 --reuse-dem --source-gpkg data/nepal/source/source_v1.gpkg
python -m ravaflow_tools.verify_derived --version v2_110 --volume 110e6
```
- `--base-m` = warstwa bazowa w całym polygonie [m]; „pół na pół” to `volume / (2 x pole polygonu)`.
  Skrypt zatrzymuje się, gdy sama baza jest większa niż objętość.
- `--release level` = sama pozioma tafla (bez bazy).
- Nowa wersja zawsze w nowym folderze (`--version`), skrypt nie nadpisuje istniejących.
- Inny polygon źródła: nowy `source_<wersja>.gpkg` z tabelą `zrodlo`, w obrębie obszaru roboczego.