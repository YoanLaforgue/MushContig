# MushContig

> `MushContig` est une méthodologie conçue pour surmonter les défis de l'identification fongique dans les échantillons complexes. En se concentrant sur le long fragment ribosomique 18S-ITS-LSU, ce pipeline offre une résolution taxonomique supérieure à celle des approches standards basées uniquement sur le gène 18S.

Développé dans un cadre clinique, il vise à fournir un diagnostic rapide et précis, offrant une alternative puissante à la culture fongique traditionnelle, notamment pour les échantillons poly-fongiques.

---

## Contexte

Le règne fongique reste l'un des règnes du vivant les moins caractérisés sur le plan génomique. Les bases de données publiques, bien que vastes, manquent souvent de génomes complets, se limitant principalement à des marqueurs courts comme la région 18S. Cette limitation pose un problème majeur : la faible distance génétique entre certaines espèces fongiques proches rend leur distinction difficile sur la base de ce seul marqueur.

Pour pallier ces limites, `MushContig` cible une région plus informative :
*   **18S rRNA** : Gène conservé, utile pour une classification à des niveaux taxonomiques élevés.
*   **ITS (Internal Transcribed Spacer)** : Région hypervariable, considérée comme le "code-barres" standard pour l'identification des espèces fongiques en raison de son taux de mutation élevé.
*   **LSU (Large Subunit) rRNA** : Offre des informations phylogénétiques complémentaires et robustes.

L'amplification de ce fragment long (environ 2700 pb) est rendue possible par la technologie de séquençage *long-read* d'**Oxford Nanopore**, qui est au cœur de ce pipeline.

---

## Prérequis

*   **Système d'exploitation** : Un système basé sur Linux/Unix avec un shell `Bash`.
*   **Logiciels** :
    *   `kraken2` (v2.1.2)
    *   `nanoplot` (v1.42.0)
    *   `porechop` (v0.2.3, avec SeqAn 2.1.1)
    *   `nanofilt` (v2.7.1)
    *   `bbmap` (v39.00)
    *   `medaka` (v1.11.3)
    *   `minimap2` (v2.28)
    *   `samtools` (v1.15)
    *   `amplicon_sorter` (v2025/05/28)
*   **Python** :
    *   `python` (v3.9)
    *   Packages : `edlib`, `biopython`, `matplotlib`

---

## Tutoriel

Ce tutoriel décrit les différentes étapes de la méthodologie `MushContig`. Les chemins (`$path/...`) et variables (`$nb_threads`, `$numBarcode`) doivent être adaptés à votre environnement.

### Étape 1 : Contrôle Qualité (QC) des *Reads*

Cette première étape évalue la qualité globale du *run* de séquençage.
```bash
NanoPlot -t "$nb_threads" \
        --fastq "$path/to/your/fastq/$numBarcode.fastq" \
        --title "${dateSeq}_${numBarcode}_QC" \
        --outdir "$path/to/output/NANOPLOT_REPORT" \
        --maxlength 10000 \
        --plots dot
```

### STEP 2 : Human Depletion

Selon la nature de certains prélèvements cliniques, on observe parfois une population de reads humaines trop importante. C’est pourquoi il est nécessaire de réaliser une déplétion humaine informatisée.

```bash
kraken2 \
     --threads "$nb_threads" \
     --db "$Human_data_base" \
     --confidence 0.1 \
     --report "$.../Report_Human_Depletion.txt" \
     --use-names \
     --output "$.../Output.txt" \
     --classified-out "$.../Classified.fasta" \
     --unclassified-out "$.../Unclassified.fastq" \
     "$.../$numBarcode.fastq"
```
### STEP 3 : ONT Adapters Trim

```bash
porechop -i "$.../Unclassified.fastq" \
    -o "$.../Unclassified_adapter_trim.fastq"
```
### STEP 4 : Quality & Length Trim

La région 18S-ITS-LSU fait environ ± 2700 pb.  
Lors de cette étape, nous allons extraire deux populations via `NanoFilt`:

- Population > Q15 pour l’assemblage des contigs.
- Population > Q25 pour le polissage des contigs.

![reads_trim](https://github.com/user-attachments/assets/fd0487bd-2fc0-4436-a5de-2c9ccf655db0)

```bash
for q in 15 25; do
    NanoFilt "$.../Unclassified_adapter_trim.fastq" \
        -q $q \
        --headcrop 10 \
        --tailcrop 10 \
        --length 2300 \
        --maxlength 3500 \
        --logfile "$.../nanofilt_Q$q.log" \
        > "$.../TRIM/Q$q/$numBarcode.Q$q.fastq"
done
```
### STEP 5 : Assembly

Pour réaliser l’assemblage, nous utilisons `Amplicon_sorter`, un outil développé pour trier les séquences selon leur similarité et leur longueur, puis construire des séquences consensus robustes.

```bash
python3 "$.../amplicon_sorter.py" \
    -i "$.../TRIM/Q15/$numBarcode.Q15.fastq" \
    -ar -aln -maxr 30000 \  
    -o "$.../Amplicon_sorter"
```
`maxr` est une variable correspondant au nombre maximal de reads en input.

## STEP 6 : Polishing

Le polissage avec `Medaka` a pour objectif de corriger un maximum d’erreurs restantes dans l’assemblage en utilisant uniquement les reads ONT.  
Le but est d’obtenir la meilleure assembly possible basée uniquement sur les données ONT.

```bash
medaka_consensus -i "$.../TRIM/Q25/$numBarcode.Q25.fastq" \
    -d "$.../Amplicon_sorter/contigs.$numBarcode.fasta" \
    -o "$.../POLISHING" \
    -m r1041_e82_400bps_sup_g615 \
    -t "$nb_threads"
```
### STEP 7 : Identification

L'identification des espèces demeure l'un des défis majeurs en métagénomique fongique. Pour surmonter les limites des bases de données (références manquantes, taxonomie incomplète), nous adoptons une **approche d'identification multiple** en comparant nos consensus finaux à plusieurs bases de données de référence :
  - `NCBI nt`
  - `MycoBank`
  - `SILVA LSU 99%`
  - `Unite 99% (v.7.2)`
  - `Unite dynamic (v.9.0)`

Cette stratégie, couplée aux **données cliniques** du patient, permet d'affiner le diagnostic et de se rapprocher au mieux de la réalité biologique de l'échantillon. Les *contigs* sont finalement renommés avec l'identification taxonomique la plus probable.

### STEP 8 : Tableau d'Abondance

Pour quantifier la proportion de chaque espèce identifiée, un alignement des *reads* initiaux (non-humains) sur les consensus finaux est réalisé avec `Minimap2`. Les résultats de cet alignement permettent de générer un tableau (`.tsv`) qui illustre l'**abondance relative** de chaque espèce, offrant une vue quantitative de la composition fongique de l'échantillon.

---

## Perspectives

Il est fort à parier que le principal défi pour la mycologie au cours de la prochaine décennie sera le séquençage d’un large panel de champignons, dans le but de constituer des bases de données plus robustes et représentatives de la diversité fongique.
Pour l’heure, une approche dite de novo, sans biais dans la reconstruction génomique, couplée à une comparaison des séquences sur différentes bases de données, constitue une alternative pertinente. Bien qu’imparfaite, cette méthode présente une fiabilité supérieure à celle reposant uniquement sur l’identification de la région 18S, laquelle ne permet généralement pas une identification au niveau de l’espèce.

`MushContig` s’inscrit dans une niche en analysant la région 18S–ITS–LSU, répondant ainsi à une demande clinique spécifique. Néanmoins, ce pipeline reste perfectible et n’est actuellement pas conçu pour le traitement à haut débit. Il est amené à évoluer avec son temps, en parallèle des avancées techniques et technologiques.

La `communauté Nanopore` est particulièrement active et propose régulièrement des outils bio-informatiques innovants. Par ailleurs, les progrès rapides du basecalling dans ce domaine, ainsi que l’émergence d’outils de correction de lectures (read correction), tendent à rapprocher la qualité des données `ONT` de celle obtenue avec la technologie `Illumina`.

---

## Poster

Concours au 20ᵉ congrès national de la Société Française de Microbiologie (SFM) qui se tiendra au Palais des Congrès de Bordeaux du 24 au 26 septembre 2025.

<img width="350" height="500" alt="Capture-decran-2025-01-17-a-10 46 49-722x1024" src="https://github.com/user-attachments/assets/0331cde4-25a5-466d-b3c8-c1f00504c6c5" />
