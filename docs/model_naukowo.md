# Jak działa model: strona naukowa

## 1. Problem fizyczny

Zdarzenie z 26 sierpnia 2026 w Langtang: zawalenie skał i lodowca przekształciło się w spływ gruzowy i powódź
na ok. 100 km wzdłuż rzek Lende Khola i Trishuli. Mamy kaskadę procesów: lawina skalno-lodowa, uwodnienie
i unoszenie materiału po drodze, przejście w spływ gruzowy, powódź. Wstępne szacunki objętości źródła mocno się
różnią (od ok. 55 mln ton lodu i skał do 200 mln m³ z inwersji sejsmicznej), więc objętość traktujemy jako
niepewną (plan przewiduje przedział ok. 20–200 mln m³).

## 2. Dwie klasy modeli w projekcie

| | SCIDDICA-SS3 | r.avaflow v4 |
|---|---|---|
| Typ | automat komórkowy, jedna faza efektywna | model ciągły, wielofazowy, open source (GRASS GIS) |
| Rola | model prosty, szybki, przejrzysty | model odniesienia fizycznego |
| Woda i lód | brak osobnych faz (jeden efektywny materiał) | rozdzielone fazy, unoszenie, zatrzymanie i depozycja |
| Ograniczenia | walidowany na mniejszych zdarzeniach | wymaga parametrów reologicznych, cięższy obliczeniowo |

(Opis wg planu projektu.) Porównanie obu modeli na jednym zdarzeniu, z identycznym protokołem, jest wartością
naukową projektu.

## 3. r.avaflow

### 3.1 Co robi

r.avaflow to narzędzie do symulacji spływów masowych (lawiny, spływy gruzowe, powodzie, kaskady procesów),
działające jako zestaw modułów GRASS GIS. Według planu i artykułu o wersji 4 (Mergili i in. 2025, GMD 18:
9879–9896) rozwiązuje równania zachowania masy i pędu uśrednione po głębokości, z unoszeniem i depozycją
materiału. Oferuje model mieszaniny typu Voellmy oraz model wielofazowy, a wersja 4 dodaje model warstwowy,
który ma ułatwiać symulację oddziaływania osuwiska ze zbiornikiem (przydatne w scenariuszu z tamą). 

### 3.2 Idea modelu uśrednionego po głębokości (ogólna wiedza, *do weryfikacji* w artykułach)

Zamiast rozwiązywać przepływ w trzech wymiarach, model opisuje masę przez grubość przepływu `h(x, y, t)` i
prędkość uśrednioną po grubości `u(x, y, t)`. Schematycznie dla jednej fazy:

- **bilans masy:** zmiana `h` w czasie = napływ i odpływ masy do sąsiednich komórek + źródło (unoszenie
  materiału z podłoża, minus depozycja),
- **bilans pędu:** zmiana pędu `h·u` = składowa grawitacji wzdłuż stoku − opór podłoża (tarcie, opór
  turbulentny/lepki) − gradient ciśnienia + człon związany z unoszeniem.

Model wielofazowy prowadzi takie bilanse osobno dla faz (w eksperymencie 9 z paczki Cascadia fazy 1, 2 i 3 mają
gęstości 2800, 1800 i 1000 kg/m³, co odpowiada skale, drobnej frakcji i wodzie) i uwzględnia oddziaływanie
między fazami. Szczegóły sformułowania i znaczenie parametrów (kąty tarcia, współczynniki oporu, lepkość,
parametry unoszenia) są w artykułach Mergili i in. 2017 (GMD 10: 553–569) oraz 2025 i w manualu.
Schemat numeryczny jest wysokiej rozdzielczości. Krok czasowy dobierany jest
adaptacyjnie wg warunku CFL; w naszych logach liczba CFL oscylowała wokół 0,37–0,40.

### 3.3 Dane wejściowe (jak je rozumie r.avaflow)

| Wejście | Parametr r.avaflow | W naszym kontrakcie |
|---|---|---|
| model wysokościowy podłoża | `elevation` | `elevation.tif` (powierzchnia **podstawy** masy) |
| masa w obszarze odrywania | `hrelease`, `hrelease1..3` | `release_height.tif` |
| materiał do unoszenia | `hentrmax`, `hentrmax1` | `entrainment_height.tif` (opcjonalnie) |
| zasięg rzeczywisty do oceny | `impactarea` | `reference_extent.tif` |
| gęstość, tarcie, lepkość itd. | `density`, `friction`, `viscosity`, ... | parametry w `params` w YAML |
| czas | `time` (krok zapisu, czas całkowity) | `params.time` |

Konwencja wysokości: suma `elevation` i `release_height` ma dawać teren sprzed zdarzenia. Poza obszarem odrywania
`release_height` ma być dokładnie 0.

### 3.4 Wyniki

Dla każdego kroku czasu: grubość przepływu (`hflow`) i jego frakcji, a także maksima w czasie (`*_max`), czas
dojścia czoła (`treach`), wykresy i animacje, plik z oceną względem `impactarea` i czas obliczeń. Układ plików
opisuje `docs/uruchamianie.md`, część 4. W logu przebiegu każdy krok wypisuje m.in. maksymalną grubość, maksymalną
prędkość, objętość i energię kinetyczną, co pozwala ocenić, kiedy masa się zatrzymuje.

## 4. Jak oceniamy wyniki (to jest nasz kod, `ravaflow_tools/metrics.py`)

### 4.1 Zasięg modelu i zasięg rzeczywisty

- **Zasięg modelu S:** komórki, w których maksymalna w czasie grubość przepływu `h_max ≥ 0,1 m` (próg ustawiany w
  konfiguracji, domyślnie 0,1 m, tak jak pierwsza wartość parametru `thresholds` w r.avaflow). Próg zmienia wynik,
  więc **musi być identyczny dla r.avaflow i SS3**.
- **Zasięg rzeczywisty R:** raster z konwencją: > 0 zasięg, 0 brak, ujemne lub brak danych = nieznane. Komórki
  nieznane są wyłączone z liczenia metryk.
- Jeśli siatka referencji różni się od siatki modelu, referencja jest próbkowana na siatkę modelu metodą
  najbliższego sąsiada (to wprowadza niewielką rozbieżność; w teście Cascadia: 0,549 wg r.avaflow vs 0,556 wg nas).

### 4.2 Metryki

Niech TP, FP, FN, TN to liczby komórek (zgodne, nadmiar modelu, niedobór modelu, zgodny brak zasięgu).

| Metryka | Wzór | Znaczenie |
|---|---|---|
| **Jaccard e** (główna) | TP / (TP + FP + FN) = \|R ∩ S\| / \|R ∪ S\| | zgodność zasięgu; nasz próg z RQ1: e ≥ 0,7. W r.avaflow to **CSI** |
| precision | TP / (TP + FP) | jaka część zasięgu modelu jest prawdziwa (niska = model zawyża) |
| recall | TP / (TP + FN) | jaka część zasięgu prawdziwego jest odtworzona (niska = model zaniża) |
| F1 | 2·TP / (2·TP + FP + FN) | średnia harmoniczna precision i recall |
| stosunek pól | (TP + FP) / (TP + FN) | w r.avaflow FoC; > 1: model zawyża powierzchnię |
| odległość Hausdorffa | największa odległość [m] między zbiorami | wrażliwa na pojedyncze odstające komórki, traktować pomocniczo |

**Uwaga o definicjach:** w pracy o SS2 (2009) „e” to pierwiastek z powyższego ilorazu, w SS3 sam iloraz. My
używamy ilorazu (zgodnie z SS3 i z CSI w r.avaflow). Nie mylić z SPI (synthetic performance index) z pliku
`*_evaluation.txt`, który jest inną miarą.

### 4.3 Czego te metryki nie mówią

- Zgodność zasięgu nie oznacza zgodności dynamiki (grubość złoża, czas dojścia czoła): to RQ2. Dobry obrys może
  wynikać z przypadkowej kombinacji parametrów (ekwifinalność, RQ4).
- Dane o dynamice (czasy, przepływy) są wstępne i medialne, więc do kontroli jakościowej, nie do kalibracji.

## 5. Założenia i ograniczenia

1. **Parametry reologiczne są słabo ograniczone.** Autorzy r.avaflow wprost zaznaczają, że potrzebne są dalsze
   badania parametrów, zanim symulacje będą wiarygodne w prognozach (Mergili i in. 2017), a w pracy o Piz Cengalo
   dwa scenariusze źródła wody nie dały się rozstrzygnąć (Mergili i in. 2020). Nasze wyniki to odtworzenie
   zdarzenia i analiza wrażliwości, nie prognoza.
2. **Dane wejściowe są niepewne:** błędy NMT w Himalajach, teren zmienił się po zdarzeniu, objętość źródła
   różni się kilkukrotnie między szacunkami.
3. **Scenariusze „co by było, gdyby” są eksploracyjne**: nie da się ich zwalidować na danych, parametry zostają
   z kalibracji na zdarzeniu bazowym, więc dochodzi niepewność ekstrapolacji.
4. **Rozdzielczość:** siatka 30 m jest gruba względem wąskich dolin; zawyżony rozlew na płaskim terenie widać już w
   teście Cascadia.
5. **Różnica modeli:** SS3 nie ma wody ani lodu, więc w scenariuszach zmieniających udział lodu i wody jego
   reakcję trzeba modelować przez parametry mobilności, co jest arbitralne i musi być jawnie zaznaczone.

## 6. Słownik

| Termin | Znaczenie |
|---|---|
| obszar odrywania | miejsce, z którego startuje masa (`release_height > 0`) |
| `hflow`, `hflow_max` | grubość przepływu w kroku czasu, jej maksimum w czasie |
| `treach` | czas dojścia czoła do komórki |
| `impactarea` | raster zasięgu rzeczywistego dla oceny w r.avaflow |
| unoszenie (entrainment) | wciąganie materiału z podłoża do przepływu |
| ekwifinalność | różne zestawy parametrów dające ten sam wynik |
| CSI | critical success index, w r.avaflow to samo co Jaccard |
| mapset / lokalizacja (GRASS) | katalog roboczy GRASS-a; nasza lokalizacja `xy_loc` ma generyczny układ XY |

