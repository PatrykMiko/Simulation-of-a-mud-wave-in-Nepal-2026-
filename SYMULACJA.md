# Jak uruchomić 20 symulacji Nepalu na siatce 30 m

Po sklonowaniu gałęzi `main` otwórz terminal w folderze z tym plikiem. Dane Nepalu i źródła r.avaflow są w repo. Wszystkie próby działają na **CPU**. Gotowym plikiem do wysłania będzie **`wyniki/paczka.zip`**.

## Wybierz sposób uruchomienia

| System komputera | Środowisko obliczeń | Terminal do poniższych poleceń |
|---|---|---|
| Windows z obsługą wirtualizacji | Ubuntu 24.04 w Canonical Multipass | **Git Bash** z Git for Windows, w katalogu klonu |
| macOS z obsługą Multipass | Ubuntu 24.04 w Multipass | Terminal, w katalogu klonu |
| Inny Linux z obsługą wirtualizacji | Ubuntu 24.04 w Multipass | Bash, w katalogu klonu |
| Natywne Ubuntu 24.04 | bezpośrednio na komputerze | Bash, w katalogu klonu |

Na Windowsie, macOS i innym Linuksie główny system zostaje na miejscu. Maszyna wirtualna ma własny dysk i zachowuje zainstalowane programy oraz dane między uruchomieniami. Zgodność architektury komputera, wersji jądra i zapisu pamięci sprawdzi test instalacyjny. W tym projekcie wykonano go na Linuksie x86-64. Przed startem na innym sprzęcie, również na Apple Silicon, poczekaj na zaliczenie testu na tym komputerze.

### Windows, macOS lub inny Linux

1. Zainstaluj **Canonical Multipass** dla swojego systemu. Na Windowsie instalator wymaga uprawnień administratora i obsługi wirtualizacji. Na macOS użyj instalatora `.pkg`. Na Linuksie z obsługą Snap można wykonać `sudo snap install multipass`. Następnie otwórz wskazany w tabeli terminal i sprawdź polecenie `multipass version`. [Instalacja i wymagania Multipass dla danego systemu](https://documentation.ubuntu.com/multipass/latest/how-to-guides/install-multipass/).
2. Sprawdź **wolną pamięć oraz dysk na fizycznym komputerze**. Przykład dla jednej próby równocześnie to maszyna z 2 CPU, 16 GiB RAM i dyskiem 80 GiB. Użyj go, jeśli możesz udostępnić tyle RAM i miejsca, zostawiając zasoby dla systemu. Dysk wirtualny może rosnąć w trakcie pracy. Wybierz rozmiar, który zmieści się na dysku fizycznym także po wypełnieniu.
3. W katalogu klonu uruchom jednorazową instalację. Przykład poniżej tworzy maszynę `nepal30`, przenosi z lokalnego commita pliki z linuksowymi uprawnieniami i instaluje zależności wewnątrz maszyny.

   ```bash
   bash scripts/eksperyment-vm.sh install 2 16G 80G
   ```

   Na stacji z większą pamięcią można przed utworzeniem maszyny wybrać np. `install 6 60G 200G` dla maksymalnie pięciu prób, jeśli gospodarz ma odpowiedni zapas. Parametry określają rozmiar tworzonej maszyny. Powtórzenie instalacji w istniejącej maszynie wykorzystuje jej dotychczasowe zasoby. Maszynę i checkpointy zachowaj do odebrania wyniku.

### Natywne Ubuntu 24.04

Otwórz Bash w katalogu repo na lokalnym dysku linuksowym. Ze względu na wrapper GRASS użyj ścieżki bez spacji i polskich znaków. Potrzebujesz konta z `sudo` oraz internetu do instalacji.

```bash
./scripts/eksperyment-install.sh
```

### Co sprawdza instalator

Obie ścieżki instalują GRASS GIS 8 z nagłówkami, GDAL, kompilator C/C++, Python z lokalnym środowiskiem, NumPy, SciPy, Rasterio, Matplotlib, Optunę i FFmpeg. Instalator kontroluje SHA-256 kodu r.avaflow 4.0G Rev.7, kompiluje solver i buduje DMTCP 4.2.0 do zapisu procesów. Instalowanie sterowników GPU lub CUDA jest zbędne.

Na małym terenie instalator sprawdza dwa sposoby uwalniania oraz zgodność z oryginalnym rdzeniem. Dwie równoległe próby są zapisywane i odtwarzane z ręcznego zapisu, a potem z osobnej kopii po dalszej zmianie plików. Poczekaj na komunikat o zgodności rastrów, klatek, metryk i bilansu oraz odczyt zasobów. Test obejmuje odtworzenie procesów w tej samej instalacji.

## Uruchomienie, postęp i wyniki

Wybierz kolumnę właściwą dla sposobu instalacji. Polecenia wykonuj w katalogu klonu na swoim komputerze.

| Czynność | Windows / macOS / Linux z Multipass | Natywne Ubuntu |
|---|---|---|
| Sprawdź zasoby przed pierwszym startem | `bash scripts/eksperyment-vm.sh zasoby` | `./scripts/eksperyment.sh zasoby` |
| Rozpocznij serię w tle | `bash scripts/eksperyment-vm.sh start` | `./scripts/eksperyment.sh start` |
| Opcjonalnie wybierz dwie próby przy pierwszym starcie | `bash scripts/eksperyment-vm.sh start 2` | `NEPAL_WORKERS=2 ./scripts/eksperyment.sh start` |
| Sprawdź postęp i ostatnią potwierdzoną kopię | `bash scripts/eksperyment-vm.sh status` | `./scripts/eksperyment.sh status` |
| **Zapisz i zatrzymaj przed wyłączeniem laptopa** | `bash scripts/eksperyment-vm.sh save-stop` | `./scripts/eksperyment.sh save-stop` |
| **Wznów po planowanym wyłączeniu** | `bash scripts/eksperyment-vm.sh start` | `./scripts/eksperyment.sh start` |
| Po awarii wróć do potwierdzonej automatycznej kopii | `bash scripts/eksperyment-vm.sh recover` | `./scripts/eksperyment.sh recover` |
| Odbierz sprawdzony ZIP po zakończeniu | `bash scripts/eksperyment-vm.sh export` | gotowy w `wyniki/paczka.zip` |

Automat dobiera równoległość do **dostępnej** RAM, wolnego dysku i CPU w środowisku obliczeń. Orientacyjny budżet to 3 GiB RAM i 10 GiB dysku wspólnie oraz 10 GiB RAM, 32 GiB dysku i jeden logiczny CPU na próbę. Dla pięciu prób daje to około **53 GiB dostępnej RAM i 170 GiB wolnego dysku**. Szacunek uwzględnia dodatkowe miejsce na automatyczną kopię. W czasie jej wymiany chwilowo istnieją dwa zapisy obok bieżących obliczeń. Rzeczywisty rozmiar zależy od liczby aktywnych prób, plików solvera i systemu plików. W maszynie wirtualnej polecenie `zasoby` widzi przestrzeń wewnątrz maszyny. Sprawdzaj również rzeczywiste miejsce po stronie gospodarza.

Przy pierwszym starcie możesz wybrać mniej pracowników. Po załadowaniu solverów sprawdź zużycie zasobów. Wznowienie odtwarza grupę z dotychczasową równoległością. Pierwsze dziesięć prób ma ustalone parametry, a następne dziesięć dobiera Optuna w dwóch grupach po pięć.

## Zapis w trakcie próby i wyłączenie laptopa

- **Planowane wyłączenie.** W dowolnym momencie aktywnych obliczeń wywołaj `save-stop` z tabeli i **poczekaj na końcowe potwierdzenie**. DMTCP zapisuje stan menedżera, aktywnych solverów i otwartych plików. Polecenie dla Multipass zatrzymuje także maszynę wirtualną po udanym zapisie. Po ponownym włączeniu `start` odtwarza ten stan. Zachowaj tę samą maszynę, ścieżkę, jądro, kod, dane i biblioteki.
- **Automatyczny zapis.** Co godzinę rzeczywistej pracy między wznowieniami oraz po ukończeniu 5., 10. i 15. próby nadzór zatrzymuje grupę, utrwala jej pamięć **razem z odpowiadającymi jej plikami roboczymi** i wznawia obliczenia, o ile jest miejsce. Kopia jest pod `wyniki/.autosave` wewnątrz środowiska obliczeń. Wcześniejsza potwierdzona kopia pozostaje dostępna do publikacji kolejnej. Po każdej ukończonej próbie osobno zapisywane są raport i ranking. Klatki co 120 s zdarzenia służą do animacji.
- **Po nagłej utracie zasilania.** `start` potrafi wznowić potwierdzony zapis ręczny. Jeśli wyświetli komunikat o kopii automatycznej, użyj `recover`. To polecenie cofa pliki i procesy do ostatniego potwierdzonego punktu, a potem podejmuje obliczenia. Praca wykonana później może wymagać powtórzenia. Przed pierwszą potwierdzoną kopią pozostają ukończone raporty i możliwość powtórzenia niedokończonych prób poleceniem `retry`.
- **Przy małej ilości miejsca.** Nadzór podejmuje próbę zapisu i zatrzymania. Stan `saved` z powodem `low-disk-after-backup` lub `backup-failed` oznacza zachowany bieżący checkpoint. Zwolnij miejsce i użyj `start`. Przy komunikacie **„Pauza w RAM”** pozostaw komputer włączony, zwolnij miejsce poza plikami serii i ponów `save-stop`.

## Ile czasu zdarzenia liczy jedna próba?

Każda próba sama dobiera moment zakończenia według tych samych reguł. Co 120 s symulowanego zdarzenia program zapisuje klatkę i sprawdza ruch fali.

- **Po dotarciu do końca badanego odcinka** program czeka jeszcze co najmniej 600 s czasu zdarzenia. Kończy próbę, gdy przez pięć kolejnych interwałów skumulowane pole zalania przyrasta za każdym razem najwyżej o 0,5%. Dojście oznacza co najmniej trzy oceniane komórki z chwilową głębokością ≥0,1 m w promieniu 300 m od punktu końcowego.
- **Gdy fala jeszcze nie dotarła**, od 7200 s program sprawdza, czy przez pełną symulowaną godzinę pole zalane kiedykolwiek poza źródłem zwiększyło się najwyżej o 0,5%, a zalanie zbliżyło się do punktu końcowego o mniej niż 30 m. Jeśli oba warunki się utrzymują, kończy próbę z powodem `stalled`. Najwcześniej wypada to około 10 800 s (3 h zdarzenia).
- **Gdy solver zgłosi ustanie przepływu**, próba może skończyć się wcześniej. Jeśli fala wciąż się przemieszcza, górny limit wynosi 43 200 s, czyli 12 godzin **modelowanego zdarzenia**.

Czas pracy komputera zależy od sprzętu, liczby równoległych prób i kroku numerycznego. Po pierwszej grupie masz pierwszy pomiar dla swojej maszyny. Końcowe IoU całego ocenianego obszaru wyznacza miejsce w rankingu. Raport podaje osobno czas symulowany, powód zakończenia i informację o dojściu.

## Co wysłać

Po komunikacie „Zakończono” wyślij **`wyniki/paczka.zip`**. Program sam tworzy archiwum i sprawdza jego CRC. `ranking.csv` jest posortowany malejąco według końcowego IoU całości. `best.json` wskazuje katalog zwycięzcy. W nim są m.in. `zasieg.png`, `glebokosc.png`, `fala.mp4`, georeferencyjne rastry, bilans i klatki z głębokością i prędkościami faz. Z każdej przegranej próby pozostaje raport JSON. Szczegółowe założenia i definicje zawiera [OPIS.md](OPIS.md).

## Jeśli instalacja lub start się nie uda

Otwórz repozytorium w **Claude Code** i poproś Claude’a o zdiagnozowanie oraz naprawienie błędu. Podaj system, używaną ścieżkę (Multipass lub natywne Ubuntu), dokładne polecenie i ostatnie linie komunikatu błędu. W przypadku błędu w tle dołącz wynik `status` oraz końcówkę `wyniki/supervisor.log` i `wyniki/checkpoint/job.log`, jeśli powstały. W Multipass logi są wewnątrz maszyny. Możesz je przeczytać np. poleceniem `multipass exec nepal30 -- tail -n 50 /home/ubuntu/nepal30/wyniki/supervisor.log`.

Przykładowa prośba do Claude’a:

> Przeczytaj SYMULACJA.md i OPIS.md. Na moim komputerze polecenie […] zakończyło się komunikatem […]. Znajdź przyczynę, popraw instalację lub uruchomienie i sprawdź naprawę. Zachowaj parametry eksperymentu, metrykę IoU, wyniki i checkpointy. Jeśli istnieje zapis aktywnej serii, przed zmianą jej środowiska sprawdź, czy nadal będzie możliwe wznowienie. Przed pierwszą pełną serią potwierdź wynik syntetycznego testu zapisu i odtworzenia.
