# Jak uruchomić 20 symulacji Nepalu na siatce 30 m

Po zakończeniu wyślij plik **`wyniki/paczka.zip`**. Zawiera ranking według IoU całego ocenianego obszaru, wyniki wszystkich 20 prób oraz mapy, animację i dane zwycięzcy.

## Udostępnienie przez GitHub

**Patryk** wykonuje jednorazowo poniższe polecenia w swojej kopii projektu. Kod eksperymentu, dokumentacja, źródła r.avaflow oraz dane są już w lokalnej historii `main`. Wystarczy wypchnąć tę gałąź.

```bash
cd "$HOME/PycharmProjects/Simulation-of-a-mud-wave-in-Nepal-2026-"
git push origin main
```

**Osoba wykonująca obliczenia** po otrzymaniu dostępu przechodzi do instalacji poniżej. Jeśli ma już klon w `$HOME/nepal30`, przed pierwszą instalacją wystarczy wykonać `cd "$HOME/nepal30"`, `git switch main` i `git pull --ff-only origin main`, a następnie uruchomić instalator. Aktualizacje repo wykonuj przed rozpoczęciem serii. Wznowienie zapisanego przebiegu wymaga tej samej wersji kodu.

## Instalacja na nowym komputerze

Instrukcja dotyczy natywnego Ubuntu 24.04 LTS, x86-64. Potrzebne są konto z `sudo`, internet do instalacji oraz lokalny dysk z systemem plików Linux. Solver używa CPU. Wszystkie polecenia po sklonowaniu wykonuj z katalogu repozytorium.

1. Na czystym Ubuntu zainstaluj Git i sklonuj repo do krótkiej ścieżki na lokalnym dysku. Jeśli repo jest prywatne, najpierw przyjmij zaproszenie właściciela na GitHubie i skorzystaj ze swojego uwierzytelnienia Git. Ścieżka projektu powinna zawierać wyłącznie litery ASCII, cyfry, myślniki, kropki, ukośniki lub podkreślenia.

   ```bash
   sudo apt-get update && sudo apt-get install -y git
   git clone --branch main https://github.com/PatrykMiko/Simulation-of-a-mud-wave-in-Nepal-2026-.git "$HOME/nepal30"
   cd "$HOME/nepal30"
   ./scripts/eksperyment-install.sh
   ```

2. Instalator używa APT do instalacji **GRASS GIS 8** i nagłówków `grass-dev`, narzędzi **GDAL**, kompilatora C/C++, narzędzi budowy, systemowego Pythona, Pillow i **FFmpeg**. Następnie tworzy lokalne środowisko `.venv` i instaluje w nim m.in. NumPy, SciPy, Rasterio, Matplotlib i Optunę. Odczytuje z repo źródła **r.avaflow 4.0G Rev.7**, sprawdza SHA-256 i kompiluje rdzeń z monitorem oraz oryginalny rdzeń porównawczy w `.local-tools/nepal20`. Buduje również **DMTCP 4.2.0** z wtyczką do zapisu otwartych plików. Dane 30 m są już w repo. PyCharm, CUDA i konto Optuny są zbędne.
3. Instalator uruchamia krótki test obu trybów uwalniania i porównuje syntetyczny przebieg z oryginalnym rdzeniem. Następnie zapisuje dwie trwające jednocześnie małe próby, odtwarza je i porównuje wyniki z przebiegiem ciągłym. Poczekaj na komunikat o zgodności rastrów, klatek, metryk i bilansu. Sam test można powtórzyć poleceniem `./scripts/eksperyment.sh verify`.

## Sprawdź zasoby przed pierwszą próbą

```bash
cd "$HOME/nepal30"
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

## Postęp, planowana przerwa i wznowienie

W nowym terminalu:

```bash
cd "$HOME/nepal30"
./scripts/eksperyment.sh status
```

Status pokazuje liczbę zakończonych prób, czas zdarzenia ostatniej zapisanej klatki aktywnych prób, stan dysku i lidera. Ranking zaczyna się pojawiać po pierwszej zakończonej próbie. Szczegóły uruchomienia są w `wyniki/supervisor.log` i `wyniki/checkpoint/job.log`.

Przed wyłączeniem komputera:

```bash
./scripts/eksperyment.sh save-stop
```

**Poczekaj na komunikat „ZAPIS UKOŃCZONY”**. Zapis obejmuje menedżera i aktywne procesy, także w środku przebiegu. Po ponownym uruchomieniu komputera:

```bash
cd "$HOME/nepal30"
./scripts/eksperyment.sh start
```

Zachowaj ten sam katalog, kernel, biblioteki, programy, dane i środowisko. Nadzór sprawdza wolne miejsce i szacuje koszt checkpointu na podstawie pamięci i otwartych plików. Przy presji dysku podejmuje próbę zapisu. Przy komunikacie „Pauza w RAM” **zostaw komputer włączony**, zwolnij miejsce poza katalogiem obliczeń i ponów `save-stop`. Nagła utrata zasilania podczas pracy może kosztować obliczenia od początku aktywnej próby. Po takim błędzie sprawdź logi. Przy wygaszonych procesach `./scripts/eksperyment.sh retry` zachowuje przyjęte raporty i powtarza niedokończone próby z zapisanymi parametrami.

## Co wysłać

Po komunikacie „Zakończono” wyślij **`wyniki/paczka.zip`**. Program sam tworzy archiwum i sprawdza jego CRC. `ranking.csv` jest posortowany malejąco według końcowego IoU całości. `best.json` wskazuje katalog zwycięzcy. W nim są m.in. `zasieg.png`, `glebokosc.png`, `fala.mp4`, georeferencyjne rastry, bilans i klatki z głębokością i prędkościami faz. Z każdej przegranej próby pozostaje raport JSON. Szczegółowe założenia i definicje zawiera [OPIS.md](OPIS.md).
