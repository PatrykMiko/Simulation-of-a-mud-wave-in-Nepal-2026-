# Stan prac: co zrobiono, jakie eksperymenty i co one oznaczają

Wyniki w częściach 2 i 3 pochodzą z przykładu szkoleniowego Cascadia i z danych syntetycznych, a nie z Nepalu.

## 1. Zrobione

| Obszar | Stan |
|---|---|
| Instalacja | GRASS GIS 8.3.2, r.avaflow 4.0G Revision 7 (2 IX 2025), R i jego pakiety, lokalizacja GRASS `xy_loc` (Ubuntu 24.04) |
| Test dymny r.avaflow | przykład Cascadia, eksperymenty 1 i 2: przeszły, wbudowane wykresy działają |
| Wrapper Pythona | `runner.py`: konfiguracja YAML → skrypt startowy → GRASS (`--exec`) → wyniki, czas, parser oceny; sprawdzony na Cascadia exp. 2 oraz danych syntetycznych|
| Metryki własne | `metrics.py`: Jaccard, precision, recall, F1, stosunek pól, Hausdorff; zgodne z CSI z r.avaflow (0,556 vs 0,549) (0,535 vs 0,525 na danych syntetycznych)|
| Wizualizacja | `viz.py`/`score.py`: mapa zgodności zasięgu, mapa maksymalnej grubości, animacja GIF |
| Walidator danych | `check_inputs.py`: sprawdza kontrakt danych; przetestowany na danych poprawnych i celowo zepsutych |
| Dane syntetyczne | `synthetic.py`: zestaw w formacie kontraktu (GeoTIFF lub ASCII) do testu łańcucha, łańcuch przeszedł end-to-end na GeoTIFF w EPSG:32645 |
| Scenariusze | `scenarios.py`: S2 (skalowanie masy), S1 (tama w NMT), nadpisanie parametrów (S3); przygotowanie rastrów przetestowane, uruchomienie serii na r.avaflow: następny etap|
| Testy automatyczne | 13 testów (`python -m pytest -q`), bez GRASS-a |
| Dokumentacja | README, instrukcja uruchomienia, opis naukowy, ten plik |

## 2. Eksperymenty na przykładzie Cascadia

### 2.1 Eksperyment 1: lawina skalna, jedna faza

Konfiguracja z paczki Cascadia (gęstość 2800 kg/m³, tarcie `35,15,-1000`, spójność 3,0, czas 200 s z zapisem co 5 s,
siatka 20 m, `-a`, wycinek obszaru `aoicoords`). Brak zasięgu rzeczywistego, więc **bez oceny zgodności**.

- Czas obliczeń modelu: **ok. 86 s**; z rysowaniem wykresów całość 2 min 19 s.
- Z logu: maksymalna prędkość rośnie do ok. 118 m/s w 35 s, energia kinetyczna osiąga maksimum ok. 4,5·10⁷ w
  35 s i po ok. 85 s spada o ponad dwa rzędy wielkości (do ok. 10⁵); później zostaje cienka warstwa o
  prędkości rzędu 16–20 m/s.
- Objętość w logu maleje z 7957,9 do 7866,0 (ok. −1,2%). Nie sprawdzaliśmy przyczyny (np. odpływ poza obszar
  obliczeń), **do wyjaśnienia**, zanim zaczniemy wnioskować o bilansie masy.
- Wnioski: model uruchamia się i daje fizycznie wiarygodny przebieg (rozpęd, hamowanie, zatrzymanie), a wykresy
  (R) generują się poprawnie.

### 2.2 Eksperyment 2: spływ gruzowy, trzy fazy, ze zasięgiem rzeczywistym

Konfiguracja z paczki Cascadia (3 fazy, `rhrelease1 = 0,5`, tarcie `35,0,0,35,0,0,0,0,0.05`, czas 200 s, siatka 20 m,
`impactarea` na siatce 5 m). **Parametry domyślne, niekalibrowane.**

| Miara | r.avaflow (`*_evaluation.txt`) | Nasz kod (`score`) |
|---|---|---|
| Jaccard e (CSI) | 0,549 | 0,556 |
| stosunek pól model/obserwacja (FoC) | 1,503 | 1,51 |
| precision / recall | rTP/rOP = 0,887 (recall) | 0,593 / 0,898 |
| F1 | | 0,715 |
| długość spływu: model / obserwacja | 2280 m / 2220 m (różnica 60 m, 3%) | |
| odległość Hausdorffa | | 424 m |
| TP / FP / FN | | 1158 / 794 / 131 komórek |
| czas modelu | 96 s (pierwszy przebieg), 99,6 s (przez wrapper) | |
| czas całkowity z wykresami | | 178,5 s |

Inne wskaźniki z pliku oceny r.avaflow: HSS 0,642, AUROC 0,882, odległość do doskonałej klasyfikacji 0,167,
SPI 0,701.

**Co to znaczy:**

1. **Łańcuch działa.** Nasz niezależny kod daje prawie ten sam Jaccard co r.avaflow (0,556 vs 0,549), a stosunek pól
   się zgadza. Próg 0,1 m i definicja zasięgu są więc ustawione tak jak w r.avaflow. Różnica 0,007 to efekt
   przepróbowania referencji z 5 m na 20 m. Ten sam kod może oceniać SS3.
2. **Model zawyża szerokość rozlewu, ale dobrze trafia długość.** Precision 0,59 przy recall 0,90 oraz błąd długości
   spływu zaledwie 3% oznaczają, że masa idzie wzdłuż koryta mniej więcej tak jak w rzeczywistości, lecz rozlewa się
   zbyt szeroko, zwłaszcza na płaskim terenie u podstawy stoku (widać to na `extent_comparison.png`).
3. **Wynik jest poniżej progu e ≥ 0,7 z RQ1, ale to niekalibrowane parametry domyślne.** To jest dokładnie punkt
   wyjścia pod kalibrację: poprawa Jaccarda ma przyjść ze zmiany parametrów, nie z kodu oceny.
4. **Nie wolno z tego wnioskować o jakości r.avaflow na zdarzeniu z Nepalu.** To przykład szkoleniowy, a sama
   odległość Hausdorffa (424 m) jest wrażliwa na pojedyncze odstające komórki.

## 3. Test łańcucha na danych syntetycznych (format kontraktu danych)

Dane z `synthetic.py` (GeoTIFF, EPSG:32645, siatka 200 × 300 komórek po 30 m), parametry z Cascadia exp. 2
(niekalibrowane), czas symulacji 400 s. Zasięg „rzeczywisty" jest wymyślony.

| Miara | r.avaflow | Nasz kod (`score`) |
|---|---|---|
| Jaccard e (CSI) | 0,525 | 0,535 |
| stosunek pól (FoC) | 1,388 | 1,392 |
| precision / recall | | 0,599 / 0,833 |
| F1 | | 0,697 |
| TP / FP / FN | | 2054 / 1377 / 411 komórek |
| odległość Hausdorffa | | 1984 m |
| czas modelu | 472 s | |

**Co to znaczy:**

1. **Łańcuch działa na danych w formacie kontraktu**: walidator zaakceptował GeoTIFF i EPSG:32645,
   r.avaflow zaimportował rastry z dużymi współrzędnymi UTM, `score` policzył metryki.
2. **Jaccard nic nie mówi o jakości modelu**, bo zasięg jest wymyślony.
3. **Różnica 0,525 vs 0,535.** Siatka wyjściowa r.avaflow ma 201 × 301 komórek, referencja 200 × 300, a róg
   modelu jest przesunięty o 5 m. Z ułamków w `*_evaluation.txt` wynika różnica ok. 30 komórek na brzegu
   zasięgu (ok. 1,3% jego pola), co jest spójne z przesunięciem siatek.
4. **Metryki zależne od TN (HSS, AUROC, D2PC, SPI) z r.avaflow nie są porównywalne z naszymi**: wynika z nich,
   że r.avaflow liczy je na mniejszej domenie (ok. 15 tys. komórek zamiast ok. 60 tys. u nas). Jaccard,
   precision, recall i F1 nie zależą od TN. Do porównania modeli używamy tylko metryk z `metrics.py`.
5. **Koszt:** 472 s dla siatki 200 × 300 i 400 s symulacji, ok. 5 razy więcej niż Cascadia (ok. 100 s).
   Obszar Langtang–Trishuli będzie wielokrotnie większy, więc trzeba ustalić obszar i rozdzielczość
   kalibracji (RQ5). `model_seconds` zapisujemy zawsze.   

