# Stan prac: co zrobiono, jakie eksperymenty i co one oznaczają

Stan na 21 października 2026. Dokument jest dziennikiem tego, co zostało sprawdzone, z konkretnymi liczbami,
oraz listą rzeczy, które jeszcze trzeba zrobić. Wszystkie wyniki w części 2 pochodzą z **przykładu szkoleniowego
Cascadia dołączonego do r.avaflow**, a nie z Nepalu.

## 1. Zrobione

| Obszar | Stan |
|---|---|
| Instalacja | GRASS GIS 8.3.2, r.avaflow 4.0G Revision 7 (2 IX 2025), R i jego pakiety, lokalizacja GRASS `xy_loc` (Ubuntu 24.04) |
| Test dymny r.avaflow | przykład Cascadia, eksperymenty 1 i 2: przeszły, wbudowane wykresy działają |
| Wrapper Pythona | `runner.py`: konfiguracja YAML → skrypt startowy → GRASS (`--exec`) → wyniki, czas, parser oceny; sprawdzony na Cascadia exp. 2 |
| Metryki własne | `metrics.py`: Jaccard, precision, recall, F1, stosunek pól, Hausdorff; zgodne z CSI z r.avaflow (0,556 vs 0,549) |
| Wizualizacja | `viz.py`/`score.py`: mapa zgodności zasięgu, mapa maksymalnej grubości, animacja GIF |
| Walidator danych | `check_inputs.py`: sprawdza kontrakt danych; przetestowany na danych poprawnych i celowo zepsutych |
| Dane syntetyczne | `synthetic.py`: zestaw w formacie kontraktu (GeoTIFF lub ASCII) do testu łańcucha |
| Scenariusze | `scenarios.py`: S2 (skalowanie masy), S1 (tama w NMT), nadpisanie parametrów (S3); przygotowanie rastrów przetestowane |
| Testy automatyczne | 11 testów (`python -m pytest -q`), bez GRASS-a |
| Dokumentacja | README, instrukcja uruchomienia, opis naukowy, ten plik |

**Nie zostało jeszcze zweryfikowane:** przebieg na danych syntetycznych w GeoTIFF przez r.avaflow, wczytywanie
GeoTIFF-ów w naszym kodzie na prawdziwym sprzęcie (testy używają ASCII), uruchomienie scenariuszy przez GRASS-a.

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
