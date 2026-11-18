# Dane wejściowe, wersja źródła `v2`

Wygenerowano: 2026-11-18 skryptem `ravaflow_tools/prepare_nepal_data.py`.

## Źródła
- DEM-y (surowe pliki): copernicus: `output_hh.tif`, alos: `output_AW3D30.tif`, srtm: `output_SRTMGL1.tif`
- Polygon źródła: `source_v1.gpkg`, warstwa `zrodlo` (SHA-256 2d8e50692e08...)
- Zasięg rzeczywisty: warstwa `zasieg`; strefa erozji: cały obszar roboczy (plik `nepal_vectors.gpkg`)
- Obszar roboczy: gotowa warstwa 'bufor_rzeki'

## Siatka
GeoTIFF, EPSG:32645, komórka 30 m, 1971x2096, identyczna dla wszystkich plików i DEM-ów.
Poza obszarem roboczym: `elevation` = 0, `release_height` = 0, `entrainment_height` = 0,
`reference_extent` = -1 (brak informacji).

## Jak powstały pliki
- `elevation.tif`: DEM sprzed zdarzenia, przeliczony na siatkę (interpolacja dwuliniowa), dziury
  (wartości < 100 m) wypełnione interpolacją (do 30 komórek), a gdy nie wystarczyła,
  wartościami z Copernicusa. Artefakty ALOS/SRTM (różnica od Copernicusa > 200 m, gdy drugi DEM zgadza się z Copernicusem w 100 m) zastąpiono Copernicusem, łata powiększona o 3 komórki.
- `release_height.tif`: warstwa bazowa 44 m w całym polygonie źródła + pozioma tafla (poziom Z) wypełniająca obniżenia terenu (DEM + baza) resztą objętości. Objętość = 1.4e+08 m3. Poza polygonem dokładnie 0.
- `reference_extent.tif`: 1 w zasięgu, 0 w obszarze roboczym poza zasięgiem, -1 poza obszarem roboczym.
- `entrainment_height.tif`: z nachylenia Copernicusa (metoda Horna): >35 st. 0.5 m, 20-35 st. 1.5 m,
  <20 st. 3.0 m; tylko w strefie erozji, poza nią 0. Wartości miąższości są założeniem empirycznym.

## Odstępstwo od opisu kontraktu
Kontrakt mówi, że `elevation + release_height` = teren sprzed zdarzenia (źródło w bliźnie oderwania).
Tu źródło jest w dnie doliny (masa nadeszła z góry), więc `elevation` = teren sprzed zdarzenia bez
odejmowania masy, a `release_height` leży na nim. Odejmowanie wykopałoby w dolinie rów.

## Statystyki
| DEM | dziur | interpolacja | z Copernicusa | załatane artefakty | h_max [m] | objętość [mln m3] |
|---|---|---|---|---|---|---|
| copernicus | 0 | 0 | 0 | 0 | 270.8 | 140.000 |
| alos | 0 | 0 | 0 | 1873 | 273.5 | 140.000 |
| srtm | 0 | 0 | 0 | 1014 | 262.6 | 140.000 |
