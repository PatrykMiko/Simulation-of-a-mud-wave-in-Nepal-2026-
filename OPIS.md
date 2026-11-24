# Rekonstrukcja fali poniżej blokady w dolinie Lhende Khola

Badamy, która z 20 konfiguracji r.avaflow 4.0G Rev.7 na siatce 30 m najlepiej odtwarza obserwowany zasięg fali od blokady do zachodniego końca kartowanego odcinka. Każda próba jest pełnym rozwiązaniem równań. Zwycięzca jest gotową symulacją. Instrukcję uruchomienia zawiera [SYMULACJA.md](SYMULACJA.md).

## Jak wybieramy zwycięzcę

Optuna maksymalizuje tę samą liczbę, którą sortuje ranking: **końcowe IoU całego ocenianego obszaru**. W rastrze modelowym komórka jest zalana, jeśli maksymalna w całym przebiegu łączna głębokość mieszaniny osiągnęła 0,1 m. Na natywnej siatce porównujemy ten zasięg z mapą obserwowaną. Bierzemy wszystkie komórki ze znanym stanem 0 lub 1 poza narzuconym poligonem źródła. IoU to TP / (TP + FP + FN). Nadmierny rozlew zwiększa FP, a pominięty fragment doliny zwiększa FN. Przy identycznej wartości wygrywa niższy numer próby.

Raportujemy osobno dojście do końca, czas zakończenia, precision, recall, nadmiar i niedobór pola oraz diagnostykę obu części doliny. Próba z najwyższym IoU zachowuje pierwsze miejsce także wtedy, gdy jej czoło zatrzymało się wcześniej. Zachodnia część doliny uczestniczy w optymalizacji. Wynik jest dopasowaniem do mapy tego zdarzenia. Niezależna ocena prędkości, głębokości i czasu dojścia wymaga osobnych obserwacji. Poprzedni częściowy przebieg przy 6000 s miał IoU całości **0,240441**, kalibracyjne **0,301716**. Porównanie z końcowym wynikiem nowej serii obejmie również zmianę horyzontu.

## Dane i założenia o źródle

Używamy czterech rastrów v2 Copernicusa w EPSG:32645, 1971 × 2096 komórek po 30 m. Teren poza domeną dostaje numeryczne przedłużenie najbliższą wysokością. Takie komórki są pomijane przy ocenie. Źródło przy blokadzie obejmuje 1766 komórek. Suma z rastra wynosi w przybliżeniu 140 000 061 m³. We wszystkich próbach zadajemy 140 mln m³ mieszaniny. To założenie dla tego etapu kaskady. Rekonstrukcja lawiny skalno-lodowej i dalszego przepływu uzasadnia sprawdzenie różnych czasów uwalniania, lecz nie wyznacza hydrogramu wypływu z naszego poligonu. [Badanie zdarzenia](https://jms.imde.ac.cn/en/article/doi/10.1007/s11629-026-1482-y?viewType=HTML).

Model wielofazowy używa fazy stałej 1 i wodnej 3, faza 2 ma początkowo udział zerowy. Gęstości faz 1 i 3 wynoszą 2700 i 1000 kg/m³. Zmiana proporcji dotyczy objętości. Wariant 50/50 oznacza po 70 mln m³. Zakres od 50% do bazowych 63,73316% objętości fazy stałej pozwala sprawdzić wrażliwość mieszaniny na koncentrację. Jest to przedział eksperymentalny, a skład faktycznej fali pozostaje niepewny. [Iverson, fizyka spływów gruzowych](https://pubs.usgs.gov/publication/70193210).

W wariancie natychmiastowym cała wysokość źródła jest dostępna od początku. W stopniowym solver stosuje własny mechanizm ekstruzji od t=0, ze stałym tempem w każdej komórce i wspólnym czasem zakończenia. Dzielenie grubości przez czas powiększony o 1 s uwzględnia początkową porcję w implementacji Rev.7. Sprawdzamy zadaną całkę i objętość początkową zapisaną przez solver. Mapa dostępnej miąższości erozyjnej oraz udział stały 0,65 w pobieranym materiale są wspólne dla prób. Stroimy współczynnik erozji. Eksperyment nie identyfikuje osobnych właściwości każdego odcinka dna.

## Dziesięć konfiguracji i dziesięć adaptacyjnych

Pierwsza konfiguracja powtarza parametry najlepszej wcześniejszej próby kalibracyjnej 120 m, użyte potem na 30 m. W tabeli podano zmiany względem próby 1. Powtarzamy wielkości fizyczne tej próby, stosując nową wspólną regułę końca.

| Nr | Parametry lub zmiana | Co sprawdzamy |
|---|---|---|
| 1 | φ=35°, δ=8,1915967643°, udział stały 0,6373316048, Manning 0,03033694194, log₁₀ erozji −8,491546965, uwalnianie 200,606630818 s | punkt odniesienia |
| 2 | uwolnienie natychmiastowe | sam przebieg czasowy zasilania |
| 3 | udział stały 0,50 | większy udział wody |
| 4 | natychmiastowo, 50% fazy stałej / 50% wody, opory i erozja jak w próbie 1 | wpływ impulsu przy składzie 50/50 (porównanie z 3), wpływ Manninga przy tym samym impulsie (porównanie z 9) |
| 5 | δ=6° | mniejsze tarcie podłoża |
| 6 | Manning 0,022 | mniejszy opór fazy płynnej |
| 7 | log₁₀ erozji −8,0 | szybsze pobieranie materiału |
| 8 | log₁₀ erozji −9,0 | wolniejsze pobieranie materiału |
| 9 | natychmiastowo, udział stały 0,50, Manning 0,025 | impuls bardziej wodnej mieszaniny |
| 10 | uwalnianie 120 s, δ=7°, φ=33°, log₁₀ erozji −8,7 | szybszy dopływ, mniejsze tarcie i słabsza erozja |

Po komplecie dziesięciu wyników Optuna TPE proponuje 11–15. Po komplecie piętnastu proponuje 16–20. W obrębie grupy propozycje oczekujące są oznaczone jako trwające (`constant_liar`), a plan jest zapisany przed startem pracowników. Zakresy adaptacyjne to φ 30–38°, δ 5–12°, udział stały od 0,50 do wartości bazowej 0,6373316047825536, Manning 0,020–0,045 w skali logarytmicznej, log₁₀ współczynnika erozji od −9 do −7,5, tryb natychmiastowy lub stopniowy, a dla stopniowego czas 60–300 s w skali logarytmicznej. Poprzednie wyniki pomogły wybrać punkt startowy. Uczenie TPE korzysta wyłącznie z prób ocenianych na 30 m tym samym protokołem.

## Kiedy kończy się symulacja

Co 120 s czasu fizycznego zapisujemy klatkę. Wszystkie konfiguracje mają te same warunki końca. Do decyzji program używa geometrii domeny, dynamiki i zasięgu modelowanego. Zgodność z mapą obserwowaną jest liczona po zakończeniu.

- **Dojście i stabilizacja** (`outlet_stable`). Zachodni koniec narysowanej osi rzeki (E 300 051,57 m, N 3 076 565,33 m) ma zero obserwowanych zalanych komórek w promieniu 300 m. Punkt kontrolny umieszczamy w środku najbliższej obserwowanej zalanej komórki, około 355 m od końca osi, na E 300 225 m, N 3 076 875 m (UTM 45N). W promieniu 300 m wokół punktu jest 294 ocenianych komórek, w tym 89 obserwowanych zalanych. Dojście rejestrujemy, gdy w jednej klatce przynajmniej trzy oceniane komórki w promieniu 300 m mają co najmniej 0,1 m głębokości chwilowej. Symulujemy jeszcze minimum 600 s od pierwszego wykrycia i wymagamy pięciu kolejnych przyrostów skumulowanego zalanego pola po najwyżej 0,5% na interwał. Czas dojścia ma rozdzielczość klatek.
- **Zatrzymanie przed końcem doliny** (`stalled`). Od 7200 s sprawdzamy, czy przez pełną godzinę symulowaną skumulowany zasięg poza źródłem zwiększył się najwyżej o 0,5%, a najmniejsza odległość zalanej komórki od punktu końcowego poprawiła się o mniej niż 30 m. Większa zmiana zeruje okno obserwacji. Najwcześniejszy taki koniec wypada około 10 800 s. To geometryczna heurystyka. Na zakolach odległość w linii prostej może słabo opisywać ruch wzdłuż koryta.
- **Ustanie przepływu** (`flow_stopped`). Akceptujemy fizyczne zakończenie zgłoszone przez solver.
- **Limit obliczeniowy** (`safety_cap`). Ruchomy wariant może liczyć się do 43 200 s (12 godzin zdarzenia, z dokładnością do kroku numerycznego). To jawny limit kosztu, dobrany do eksperymentu. Jest sześć razy dłuższy od poprzedniego horyzontu. Zakończenie na limicie trzeba tak oznaczyć w publikacji. Jeśli żaden wariant nie dotrze do końca, raport pokaże to wprost.

Prywatna modyfikacja Rev.7 odczytuje pola po kroku fizycznym, stosuje tę regułę końcową i ogranicza pośredni eksport na dysk. Równania przepływu i sposób obliczania końcowych pól pochodzą z Rev.7. Test porównuje `hflow_max` z oryginalnym rdzeniem skompilowanym na tej samej maszynie. [Źródła r.avaflow](https://www.landslidemodels.org/r.avaflow/software/r.avaflow.40G.zip).

## Jak odczytać wyniki

`ranking.csv` podaje IoU, czas zdarzenia, dojście i powód zakończenia każdej próby. `records/NN.json` dodatkowo zachowuje dokładne parametry, bilans i tablice metryk. `protocol.json` zawiera skróty danych i wersję środowiska. Duże pliki przegranych są usuwane po przyjęciu raportu. Po zakończeniu wszystkich prób program renderuje zwycięzcę z zachowanych klatek, bez ponownego liczenia równań.

Wyniki zwycięzcy obejmują mapę zgodności TP/FP/FN, maksymalną głębokość, film chwilowego zalania, rastry GeoTIFF `hflow_max`, `treach`, `vfront`, `vhmax`, `basechange_fin` oraz szeregi czasowe i parametry zastosowane przez solver. Pliki `frames/NNNNN.bin` zawierają rekordy little-endian: indeks komórki int32 w porządku wierszowym od północy, głębokość float32 [m], prędkości faz 1 i 3 float32 [m/s]. Zapis obejmuje komórki z głębokością ≥0,1 m we wnętrzu siatki obliczeniowej. `timeline.csv` ma kolumny bez nagłówka: numer klatki, czas [s], liczba zalanych kiedykolwiek ocenianych komórek, liczba aktualnie zalanych ocenianych komórek przy ujściu, czas pierwszego wykrycia (−1 do czasu dojścia), stan. Czas trwania 20 prób na komputerze kolegi można oszacować po pierwszej piątce. Przy pięciu równoległych procesach czas całej serii jest w przybliżeniu sumą czasów najwolniejszej próby z każdej z czterech grup.
