# Jak uruchomić 20 symulacji Nepalu na siatce 30 m

Po zakończeniu wyślij plik **`wyniki/paczka.zip`**. Zawiera ranking według IoU całego ocenianego obszaru, wyniki wszystkich 20 prób oraz mapy, animację i dane zwycięzcy.

## Pierwsze uruchomienie po sklonowaniu repozytorium

Instrukcja zakłada świeży klon gałęzi `main` na natywnym Ubuntu 24.04 LTS (x86-64). Potrzebujesz konta z `sudo`, internetu na czas instalacji i lokalnego dysku z systemem plików Linux. Solver używa CPU.

1. Otwórz terminal **w katalogu sklonowanego repozytorium**, w którym widzisz ten plik i folder `scripts`. W tej samej lokalizacji będziesz uruchamiać wszystkie poniższe polecenia. Ze względu na wrapper GRASS wybierz ścieżkę bez spacji, złożoną z liter ASCII, cyfr, myślników, kropek, ukośników lub podkreśleń. Uruchom jednorazową instalację.

   ```bash
   ./scripts/eksperyment-install.sh
   ```

2. Instalator używa APT do instalacji **GRASS GIS 8** i nagłówków `grass-dev`, narzędzi **GDAL**, kompilatora C/C++, narzędzi budowy, systemowego Pythona, Pillow i **FFmpeg**. Następnie tworzy lokalne środowisko `.venv` i instaluje w nim m.in. NumPy, SciPy, Rasterio, Matplotlib i Optunę. Odczytuje z repo źródła **r.avaflow 4.0G Rev.7**, sprawdza SHA-256 i kompiluje rdzeń z monitorem oraz oryginalny rdzeń porównawczy w `.local-tools/nepal20`. Buduje również **DMTCP 4.2.0** z wtyczką do zapisu otwartych plików. Dane 30 m są już w repo. PyCharm, CUDA i konto Optuny są zbędne.
3. Instalator uruchamia krótki test obu trybów uwalniania i porównuje syntetyczny przebieg z oryginalnym rdzeniem. Następnie zapisuje dwie trwające jednocześnie małe próby, odtwarza je i porównuje wyniki z przebiegiem ciągłym. Poczekaj na komunikat o zgodności rastrów, klatek, metryk i bilansu. Sam test można powtórzyć poleceniem `./scripts/eksperyment.sh verify`.

## Sprawdź zasoby przed pierwszą próbą

```bash
./scripts/eksperyment.sh zasoby
```

Polecenie wyświetla bieżącą **dostępną** pamięć RAM, wolne miejsce na partycji projektu i liczbę dostępnych procesorów logicznych. Przy pierwszym starcie program ponawia odczyt i automatycznie dobiera **od jednej do pięciu** równoległych prób. Orientacyjnie rezerwuje 3 GiB RAM i 8 GiB dysku dla wspólnych operacji, następnie liczy po 10 GiB dostępnej RAM, 18 GiB wolnego dysku i jednym logicznym CPU na próbę. Dla pięciu daje to około **53 GiB dostępnej RAM, 98 GiB wolnego dysku i 5 CPU**. Podstawą szacunku dysku jest poprzedni pojedynczy finał z chwilowym budżetem checkpointu około 15,5 GiB. Rzeczywiste wykorzystanie zależy od zasięgu fali, fazy eksportu i innych aplikacji. Przy oszacowaniu zero automatyczny start wstrzymuje uruchomienie.

**Jeśli pięć naraz to za dużo**, możesz od początku ustawić np. dwie, nawet gdy automat proponuje więcej. Wybierz tylko jedno z poleceń startowych poniżej. Po uruchomieniu zajrzyj do `./scripts/eksperyment.sh status`, `free -h` i `df -h .`, szczególnie po załadowaniu pierwszych solverów. W razie potrzeby zapisz grupę przed wyłączeniem. Wznowienie z checkpointu odtwarza tę samą grupę procesów i jej dotychczasową równoległość.

## Start jednym poleceniem

W katalogu projektu wykonaj poniższe polecenie.

```bash
./scripts/eksperyment.sh start
```

Program działa w tle. Terminal można zamknąć. Aby od początku ograniczyć obciążenie do dwóch prób, jako polecenie startowe wybierz `NEPAL_WORKERS=2 ./scripts/eksperyment.sh start`. Ręczny limit ustal na poziomie sugerowanym przez `zasoby` lub niższym. Próby i momenty aktualizacji Optuny pozostają takie same. Pierwsze dziesięć konfiguracji zajmuje dwie serie, następnie Optuna proponuje pięć konfiguracji na podstawie dziesięciu wyników i pięć na podstawie piętnastu.

## Ile czasu zdarzenia liczy jedna próba?

Każda próba sama dobiera moment zakończenia według tych samych reguł. Co 120 s symulowanego zdarzenia program zapisuje klatkę i sprawdza ruch fali.

- **Po dotarciu do końca badanego odcinka** program czeka jeszcze co najmniej 600 s czasu zdarzenia. Kończy próbę, gdy przez pięć kolejnych interwałów skumulowane pole zalania przyrasta za każdym razem najwyżej o 0,5%. Dojście oznacza co najmniej trzy oceniane komórki z chwilową głębokością ≥0,1 m w promieniu 300 m od punktu końcowego.
- **Gdy fala jeszcze nie dotarła**, od 7200 s program sprawdza, czy przez pełną symulowaną godzinę pole zalane kiedykolwiek poza źródłem zwiększyło się najwyżej o 0,5%, a zalanie zbliżyło się do punktu końcowego o mniej niż 30 m. Jeśli oba warunki się utrzymują, kończy próbę z powodem `stalled`. Najwcześniej wypada to około 10 800 s (3 h zdarzenia).
- **Gdy solver zgłosi ustanie przepływu**, próba może skończyć się wcześniej. Jeśli fala wciąż się przemieszcza, górny limit wynosi 43 200 s, czyli 12 godzin **modelowanego zdarzenia**.

Czas pracy komputera zależy od sprzętu, liczby równoległych prób i kroku numerycznego. Po pierwszej grupie masz pierwszy pomiar dla swojej maszyny. Końcowe IoU całego ocenianego obszaru wyznacza miejsce w rankingu. Raport podaje osobno czas symulowany, powód zakończenia i informację o dojściu.

## Postęp, planowana przerwa i wznowienie

W nowym terminalu otwartym w katalogu repozytorium sprawdzisz postęp poleceniem

```bash
./scripts/eksperyment.sh status
```

Status pokazuje liczbę zakończonych prób, czas zdarzenia ostatniej zapisanej klatki aktywnych prób, stan dysku i lidera. Ranking zaczyna się pojawiać po pierwszej zakończonej próbie. Szczegóły uruchomienia są w `wyniki/supervisor.log` i `wyniki/checkpoint/job.log`.

Przed wyłączeniem komputera uruchom zapis.

```bash
./scripts/eksperyment.sh save-stop
```

**Poczekaj na komunikat „ZAPIS UKOŃCZONY”**. Zapis obejmuje menedżera i aktywne procesy, także w środku przebiegu. Po ponownym uruchomieniu komputera otwórz terminal w tym samym katalogu repozytorium i wznów serię.

```bash
./scripts/eksperyment.sh start
```

Zachowaj ten sam katalog, kernel, biblioteki, programy, dane i środowisko. Nadzór sprawdza wolne miejsce i szacuje koszt checkpointu na podstawie pamięci i otwartych plików. Przy presji dysku podejmuje próbę zapisu. Przy komunikacie „Pauza w RAM” **zostaw komputer włączony**, zwolnij miejsce poza katalogiem obliczeń i ponów `save-stop`. Nagła utrata zasilania podczas pracy może kosztować obliczenia od początku aktywnej próby. Po takim błędzie sprawdź logi. Przy wygaszonych procesach `./scripts/eksperyment.sh retry` zachowuje przyjęte raporty i powtarza niedokończone próby z zapisanymi parametrami.

## Co wysłać

Po komunikacie „Zakończono” wyślij **`wyniki/paczka.zip`**. Program sam tworzy archiwum i sprawdza jego CRC. `ranking.csv` jest posortowany malejąco według końcowego IoU całości. `best.json` wskazuje katalog zwycięzcy. W nim są m.in. `zasieg.png`, `glebokosc.png`, `fala.mp4`, georeferencyjne rastry, bilans i klatki z głębokością i prędkościami faz. Z każdej przegranej próby pozostaje raport JSON. Szczegółowe założenia i definicje zawiera [OPIS.md](OPIS.md).
