# HUBSON
**Parametry globalne i przygotowanie środowiska**

Wszystkie warstwy rastrowe zostały ujednolicone do wspólnego układu współrzędnych EPSG:32645. Zastosowano stałą rozdzielczość przestrzenną 30 metrów oraz rygorystycznie wymuszono identyczną siatkę dla każdego pliku.

**Masa obszaru odrywania (release\_height.tif)**

Obszar źródłowy został wyznaczony na podstawie wektorowego poligonu, któremu przypisano atrybut grubości generujący docelową objętość 72.45 milionów metrów sześciennych (to zostanie zmienione - zarówno miejsce jak i grubość być może też). Rasteryzację przeprowadzono z wykorzystaniem zasięgu bazowego modelu terenu. Tło poza strefą obrywu zostało uzupełnione twardymi zerami, aby nie było nieporządanych wartości NoData i NaN.

**Model terenu powierzchni poślizgu (elevation.tif)**

Geometrię terenu oparto na trzech niezależnych modelach numerycznych pobranych i dociętych do obszaru roboczego: Copernicus COP30, JAXA ALOS (AW3D30) oraz SRTMGL1. Aby uzyskać fizyczną powierzchnię podstawy masy (teren po oderwaniu się lawiny), od oryginalnych modeli odjęto warstwę masy obrywu przy użyciu Kalkulatora rastrów. Operacja ta została powtórzona dla każdego z trzech modeli z osobna, co wygenerowało trzy idealnie spasowane zestawy topografii.

**Rzeczywisty zasięg spływu (reference\_extent.tif)**

W celu dokładnego odwzorowania drogi spływu w dolinie rzeki, wyznaczono poligon zasięgu na podkładzie zdjęć satelitarnych. Pozwoliło to na precyzyjne wyodrębnienie nowego koryta rzeki po przejściu fali. Utworzony poligon został zrasteryzowany z przypisaniem wartości 1 dla koryta oraz 0 dla tła. Zasięg siatki podczas konwersji wymuszono bezpośrednio z pliku elevation, co zapewniło idealne nałożenie się pikseli na siebie.

**Materiał do unoszenia (entrainment\_height.tif)**

Warstwa erodowalnego podłoża została opracowana w oparciu o mapę spadków (Slope) wygenerowaną z modelu Copernicus. W Kalkulatorze rastrów QGIS wprowadzono formułę logiczną zagnieżdżającą przypisane empirycznie miąższości: 3 metry dla nachyleń do 20 stopni (aluwia w dnach dolin), 1.5 metra dla nachyleń 20-35 stopni oraz 0.5 metra dla stromizn powyżej 35 stopni. Wynikowy raster został ograniczony wyłącznie do strefy erozji wyznaczonej poligonem.

**Walidacja krzyżowa**

Ostatnim etapem było sprawdzenie zestawów danych w skrypcie walidacyjnym w Pythonie. Każdy z folderów bazy (Copernicus, ALOS, SRTM) przeszedł testy obejmujące zgodność układu, brak wartości pustych pod i nad obszarem odrywania, wymiary siatki oraz weryfikację obliczonej objętości, osiągając pełną zgodność z narzuconym kontraktem danych.

# MARSON


# Dane do kalibracji i parametry referencyjne

Dokument zawiera zbiór wartości weryfikacyjnych i parametrów wejściowych dla modelu spływu gruzowo-lodowego w dolinie Gyirong (Donglin Zangbo), na podstawie dostępnej literatury oraz analizy geomorfologicznej.

## 1. Zdarzenie początkowe i bilans masy

*   **Całkowita objętość obrywu:** ok. $1.1 \times 10^8 \text{ m}^3$
*   **Udział lodu w masie początkowej:** $10\% - 15\%$
*   **Materiał skalny:** Głównie wysokostopniowe skały metamorficzne (gnejsy).
*   **Mechanizm transformacji:** Potężne uderzenie spowodowało natychmiastowe sproszkowanie materiału skalnego, stopienie lodu i przekształcenie całości w wysoce upłynnioną powódź błotno-gruzową.
*   **Łączna objętość materiału zerwanego (entrainment):** $5.86 \times 10^7 \text{ m}^3$

## 2. Pomiary terenowe (Markery do kalibracji)

*   **Czas i dystans:** Zjawisko pokonało dystans $23 \text{ km}$ i spadek ponad $3300 \text{ m}$ w zaledwie około $7 \text{ minut}$. Maksymalny dystans oddziaływania wyniósł ponad $100 \text{ km}$.
*   **Wysokość nabiegania (Superelevation):** Na ostrym zakręcie rzeki Donglin Zangbo siła odśrodkowa wyrzuciła płynącą masę na wysokość aż $180 \text{ m}$ na wklęsłym brzegu.
*   **Zarejestrowane prędkości fali:**
    *   ok. $67 \text{ m/s}$ tuż po uderzeniu w główną dolinę (Punkt kontrolny 1)
    *   ok. $56 \text{ m/s}$ na dużym zakręcie rzeki (strefa superelewacji) (Punkt kontrolny 2)
    *   ok. $22 \text{ m/s}$ przy uderzeniu w port Gyirong (Punkt kontrolny 3)

> **Źródło referencyjne:** *The Gyirong debris flow disaster triggered by the ice-rock avalanche at Langtang Lirung Peak, Nepal, on August 26, 2026*

---

## 3. Entrainment material (Masa porywana)

Z uwagi na brak sztywno przyjętych norm dotyczących zależności między nachyleniem terenu a ilością zalegającego materiału, miąższość warstwy erodowanej ($h_{ent}$) została wyznaczona empirycznie za pomocą poniższej formuły warunkowej w kalkulatorze rastrów (dla warstwy spadków terenowych w QGIS):

```text
(("spadki@1" > 35) * 0.5) + (("spadki@1" >= 20 AND "spadki@1" <= 35) * 1.5) + (("spadki@1" < 20) * 3.0)
```

### Tabela zbiorcza charakteryzacji gleb i głębokości zrywania ($h_{ent}$)

| Wysokość n.p.m. | Nachylenie terenu | Typ gleby (US Soil Taxonomy) | Rodzaj osadu i charakterystyka | Przyjęta głębokość zrywania ($h_{ent}$) |
| :--- | :--- | :--- | :--- | :--- |
| $>4500 \text{ m}$ | Strome ($>35^\circ$) lub strefy zerozowane | Entisols | Płytkie, wysoce szkieletowe gleby na depozytach glacjalnych (skała macierzysta: paragnejsy). Znikoma zawartość iłu, słabe wietrzenie chemiczne. | $0.2 - 0.5 \text{ m}$ |
| $3000 - 4500 \text{ m}$ | Średnie do stromych ($20^\circ - 35^\circ$) | Spodosols / Inceptisols | Strefa subalpejska i zarośli alpejskich. Osady glacjalne, fluwioglacjalne oraz rumosz koluwialny na podnóżach stoków. | $1.0 - 1.5 \text{ m}$ |
| $2000 - 3000 \text{ m}$ | Średnie ($20^\circ - 35^\circ$) | Inceptisols (miejscami Alfisols) | Pasmo średniogórskie. Wyższa intensywność wietrzenia chemicznego, większy udział frakcji ilastej, osady koluwialne i deluwialne. | $1.5 - 2.0 \text{ m}$ |
| Dna dolin / koryta | Łagodne ($<20^\circ$) | Entisols (Fluvents) / Inceptisols | Nieskonsolidowane osady aluwialne, głazy, żwiry i piaski w korycie cieku oraz na tarasach zalewowych. Bardzo wysoki potencjał erozji wgłębnej. | $2.0 - 4.0 \text{ m}$ (śr. $3.0 \text{ m}$) |

### Uzasadnienie geomorfologiczne i gleboznawcze

*   **Piętro wysokogórskie ($>4500 \text{ m n.p.m.}$ oraz zbocza $>35^\circ$):**
    Ze względu na intensywne procesy grawitacyjne oraz gwałtowną erozję, na stokach o nachyleniu przekraczającym $35^\circ$ nie dochodzi do akumulacji dojrzałej pokrywy glebowej (Heimsath et al., 2012). Dominują tam płytkie gleby inicjalne typu Entisols, rozwinięte bezpośrednio na morenach i skałach paragnejsowych. Liniowy spadek intensywności wietrzenia chemicznego wraz z wysokością skutkuje znikomą zawartością frakcji ilastej, co przekłada się na brak spójności gruntu ($c \approx 0 \text{ kPa}$) oraz wysoki kąt tarcia wewnętrznego ($\phi = 34^\circ - 36^\circ$).
*   **Piętro średniogórskie ($2000 - 4500 \text{ m n.p.m.}$ oraz zbocza $20^\circ - 35^\circ$):**
    W strefie subalpejskiej i leśnej dominuje wietrzenie umiarkowane, prowadzące do wykształcenia gleb typu Spodosols oraz Inceptisols. Nachylenia terenu w tym przedziale stanowią strefę transportu i częściowej akumulacji osadów koluwialnych, gdzie miąższość erodowalnego materiału osiąga wartości $1.0 - 2.0 \text{ m}$.
*   **Strefa dna dolin i koryta rzecznego ($<20^\circ$):**
    Dna dolin w Himalajach stanowią główny rezerwuar nieskonsolidowanych materiałów aluwialnych (Fluvents). Podczas przejścia fali spływu gruzowego (debris flow) ponad $70\% - 80\%$ mobilizowanej masy pochodzi z erozji koryta i podcinania brzegów (Cook et al., 2022). Warstwy aluwialne osiągają miąższość kilku metrów, co uzasadnia przyjęcie maksymalnego limitu erozji $h_{ent} = 3.0 \text{ m}$.

---

## 4. Parametry fizykochemiczne podłoża do modeli symulacyjnych

Wartości referencyjne przeznaczone do oprogramowania symulującego dynamikę lawin i spływów (np. r.avaflow, RAMMS, SS3):

*   **Gęstość fazy stałej ($\rho_s$):** $2650 - 2700 \text{ kg/m}^3$ (wynikająca z mineralogii paragnejsów – kwarc, skalenie, miki).
*   **Gęstość objętościowa podłoża / entrainment ($\rho_{ent}$):** $1400 - 1700 \text{ kg/m}^3$ (przeliczona z baz danych SoilGrids dla osadów szkieletowych i aluwialnych).
*   **Gęstość początkowa masy obrywanej ($\rho_{release}$):** $2000 \text{ kg/m}^3$ (wartość szacunkowa do testów).
*   **Kąt tarcia wewnętrznego ($\phi$):** $32^\circ - 35^\circ$.
*   **Kąt tarcia przy dnie / basal friction ($\delta$):** $15^\circ - 20^\circ$ (uwzględniający wysokie ciśnienie porowe wody podczas powodzi/spływu).

> **Dodatkowe źródła:** 
> *   *Debris-flow Hazards and Related Phenomena*
> *   *Soil formation*
> *   *Soils of the high mountain region of Eastern Nepal: classification, distribution and soil forming processes*