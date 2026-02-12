# MushContig

> `MushContig` est une méthodologie conçue pour surmonter les défis de l'identification fongique dans les échantillons complexes. En se concentrant sur le long fragment 18S-ITS-LSU, ce pipeline offre une résolution taxonomique supérieure à celle des approches standards basées uniquement sur le gène 18S.

Développé dans un cadre clinique, il vise à fournir un diagnostic rapide et précis, offrant une alternative à la culture fongique traditionnelle, notamment pour les échantillons poly-fongiques.

---

## Contexte

Le règne fongique reste l'un des règnes du vivant les moins caractérisés sur le plan génomique. Les bases de données publiques, bien que vastes, manquent souvent de génomes complets, se limitant principalement à des marqueurs courts comme la région 18S. Cette limite pose un problème majeur : la faible distance génétique entre certaines espèces fongiques proches rend leur distinction difficile sur la base de ce seul marqueur.

Pour pallier ces limites, `MushContig` exploite la technologie long-read pour séquencer sans fragmentation un fragment incluant:
*   **18S rRNA**
*   **ITS (Internal Transcribed Spacer)** 
*   **LSU (Large Subunit) rRNA** 

---

## Prérequis

*   **Système d'exploitation** : Linux/Unix (shell `Bash`).
*   **Outils principaux** :
    *   `kraken2` (v2.1.2)
    *   `nanoplot` (v1.42.0)
    *   `porechop` (v0.2.3, avec SeqAn 2.1.1)
    *   `nanofilt` (v2.7.1)
    *   `bbmap` (v39.00)
    *   `medaka` (v1.11.3)
    *   `minimap2` (v2.28)
    *   `samtools` (v1.15)
    *   `amplicon_sorter` (v2025/05/28)
    *   `python` (v3.9)
    *   Packages : `edlib`, `biopython`, `matplotlib`
*   **Scripts .py** :
    *   `rename_fasta_with_taxid.py` : Renomme les contigs FASTA selon leur taxonomie.
    *   `base_accuracy_calculator.py` : Calcule le taux d'erreur par base à partir de la qualité médiane des reads.
    *   `plot.py` : Génère des graphiques à barres à partir des fichiers d'alignement.

---

## Tutoriel

Note : Les variables (`$path/...` , `$nb_threads`, `$numBarcode`) doivent être adaptés à votre configuration locale.

### Étape 1 : Contrôle Qualité Initial (QC) des *Reads*

Évaluation de la qualité globale du run de séquençage.

```bash
NanoPlot -t "$nb_threads" \
        --fastq "$path/to/your/fastq/$numBarcode.fastq" \
        --title "${dateSeq}_${numBarcode}_QC" \
        --outdir "$path/to/output/NANOPLOT_REPORT" \
        --maxlength 10000 \
        --plots dot
```

### Étape 2 : Human Depletion

Selon la nature de certains prélèvements cliniques, on observe parfois une population de reads humaines trop importante. C’est pourquoi il est nécessaire de réaliser une déplétion humaine informatisée.

Base de données recommandée : [Human database](https://zenodo.org/records/8339700)

```bash
kraken2 --threads "$nb_threads" --db "$Human_data_base" --confidence 0.1 \
     --report "Report_Human_Depletion.txt" --use-names --output "Output.txt" \
     --classified-out "Human_Reads.fasta" --unclassified-out "Unclassified.fastq" \
     "$numBarcode.fastq"
```
Les *reads* non classifiés (`Unclassified_non_human.fastq`) correspondent aux séquences non humaines qui seront utilisées pour la suite de l'analyse.

### Étape 3 : Suppression des Adaptateurs ONT

Utilisation de `porechop` pour supprimer les séquences d'adaptateurs résiduelles.

```bash
porechop -i "Unclassified.fastq" -o "Unclassified_adapter_trim.fastq"
```

### Étape 4 : Filtrage par Qualité et Longueur

La région 18S-ITS-LSU fait ± 2700 pb.  

Population > Q15 pour l’assemblage des contigs.

<img width="1315" height="495" alt="Capture d’écran 2026-02-11 201808" src="https://github.com/user-attachments/assets/41fb3b50-9e17-4b31-995e-071557099940" />


```bash
NanoFilt "Unclassified_adapter_trim.fastq" -q 15 --headcrop 10 --tailcrop 10 \
         --length 2300 --maxlength 3000 > "$numBarcode.Q15.fastq"
```
### Étape 5 : QC Post-filtrage

Un second passage via `NanoPlot` permet de vérifier le nombre de reads restants, le N50 et la Median read quality.

### Étape 6 : Assemblage des Contigs

Pour réaliser l’assemblage, nous utilisons `Amplicon_sorter`, un outil développé pour trier les séquences selon leur similarité et leur longueur pour générer des contigs robustes.

```bash
python3 amplicon_sorter.py -i "$numBarcode.Q15.fastq" -maxr 30000 -ldc 20 \
        -sc $Base_Accuracy -o "Amplicon_sorter_Output"
```

`-maxr` ou `--maxreads` :
Définit le nombre maximal de reads utilisés en entrée.
Amplicon_sorter étant relativement coûteux en temps de calcul, cette limitation permet de réduire la durée d’exécution. La pertinence de cette valeur peut être validée à l’aide d’une courbe de raréfaction.

'-sc' ou `--similar_consensus` : 
Correspond au seuil de similarité requis pour fusionner des groupes de reads pour générer une séquence consensus.
Ce paramètre est ajusté en fonction de la qualité médiane des reads (Median read quality). La valeur est extraite du fichier texte généré par NanoPlot, puis convertie en taux d’erreur par base, afin d’obtenir un seuil de similarité cohérent avec la qualité des données du run.

### Étape 7 : Identification

L’identification des séquences est réalisée à l’aide de l'outil `VSEARCH`, qui permet d’effectuer des recherches de similarité contre plusieurs bases de données.

```bash
vsearch --usearch_global "$INPUT_FASTA" --db "${DIRS[ITS_DB]}" --id 0.98 \
        --strand both --maxaccepts 1 --maxrejects 0 \
        --userfields query+target+id+alnlen+qcov --blast6out "$VSEARCH_OUT"
```

L'identification des espèces demeure l'un des défis majeurs en métagénomique fongique. Pour surmonter les limites des bases de données (références manquantes, taxonomie incomplète), nous adoptons une **approche d'identification multiple** en comparant nos consensus finaux à plusieurs bases de données de référence :

  - `NCBI nt`
  - `MycoBank`
  - `SILVA LSU 99%`
  - `Unite 99% (v.7.2)`
  - `Unite dynamic (v.9.0)`

Cette stratégie permet d'affiner le diagnostic et de se rapprocher au mieux de la réalité biologique de l'échantillon. Les *contigs* sont finalement renommés avec l'identification taxonomique la plus probable avec le script python `rename_fasta_with_taxid.py`.


### Étape 8 : Tableau d'Abondance

L'abondance de chaque espèce est estimée par le ré-alignement des reads initiaux sur les contigs finaux avec `Minimap2`.

```bash
# Alignement et conversion au format BAM
minimap2 -ax map-ont -t 32 "$CONTIGS.fasta" "$READS.fastq" | \
samtools sort -@ 32 -o "sorted.bam" -

# Extraction des alignements primaires et statistiques
samtools view -b -F 2304 -@ 32 "sorted.bam" > "primary.sorted.bam"
samtools idxstats "primary.sorted.bam" > "abundance_table.txt"
```
Exemple :

<img width="800" height="500" alt="barcode94 plot" src="https://github.com/user-attachments/assets/50482ca6-a62a-4d34-8c7e-1a76b246ecbd" />

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
