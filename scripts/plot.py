import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def process_taxonomy_data(input_txt, output_tsv, output_plot):
    """
    Traite un fichier texte de données taxonomiques pour extraire les informations,
    regroupe les taxonomies doublons, calcule les pourcentages de reads mappés 
    et génère un fichier TSV ainsi qu'un plot à barres.

    Args:
        input_txt (str): Chemin vers le fichier .txt d'entrée.
        output_tsv (str): Chemin vers le fichier .tsv de sortie.
        output_plot (str): Chemin vers le fichier .png de sortie.
    """
    processed_data = []
    with open(input_txt, 'r') as f:
        for line in f:
            columns = line.strip().split('\t')
            if not columns or columns[0] == '':
                continue
            
            if columns[0] == '*':
                taxonomy = 'Unclassified'
                try:
                    mapped_reads = int(next(filter(lambda x: x.isdigit() and int(x) > 0, reversed(columns))))
                except (StopIteration, ValueError):
                    mapped_reads = 0
                processed_data.append([taxonomy, mapped_reads])
                continue

            col1 = columns[0]
            if '/' in col1:
                taxonomy = col1.split('/')[-1]
            else:
                taxonomy = 'unidentified'
            
            try:
                mapped_reads = int(columns[2])
            except (IndexError, ValueError):
                continue

            processed_data.append([taxonomy, mapped_reads])

    df = pd.DataFrame(processed_data, columns=['Taxonomy', 'Mapped_Reads'])
    
    if df.empty:
        print("Aucune donnée valide n'a été extraite du fichier .tsv d'entrée.")
        return

    df = df.groupby('Taxonomy')['Mapped_Reads'].sum().reset_index()

    total_reads = df['Mapped_Reads'].sum()
    if total_reads > 0:
        df['Percentage'] = (df['Mapped_Reads'] / total_reads) * 100
    else:
        df['Percentage'] = 0

    df.to_csv(output_tsv, sep='\t', index=False, header=True)
    print(f"Fichier TSV de sortie créé : {output_tsv}")

    df_sorted = df.sort_values(by='Percentage', ascending=False)
    
    plt.figure(figsize=(14, 10))
    
    ax = sns.barplot(x='Percentage', y='Taxonomy', data=df_sorted, palette='viridis')
    
    for bar in ax.patches:
        ax.text(
            bar.get_width() + 0.5,
            bar.get_y() + bar.get_height() / 2,
            f"{bar.get_width():.2f}%",
            ha='left',
            va='center',
            fontsize=10,
            color='black'
        )
    # ----------------------------------------------------------------------
    
    plt.title('Pourcentage de Reads Mappés par Taxonomie', fontsize=18)
    plt.xlabel('Pourcentage (%)', fontsize=14)
    plt.ylabel('Taxonomie', fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    
    # Ajuster les limites de l'axe X pour laisser de la place aux pourcentages
    if not df_sorted.empty:
        ax.set_xlim(right=df_sorted['Percentage'].max() * 1.12) # Marge de 12%
    
    plt.tight_layout()
    
    plt.savefig(output_plot)
    print(f"Graphique de sortie créé : {output_plot}")


def main():
    parser = argparse.ArgumentParser(
        description="Script pour analyser les données taxonomiques, générer un fichier TSV et un plot."
    )
    parser.add_argument(
        '--input_txt', 
        type=str, 
        required=True, 
        help="Chemin vers le fichier .txt d'entrée."
    )
    parser.add_argument(
        '--output_tsv', 
        type=str, 
        required=True, 
        help="Chemin vers le fichier .tsv de sortie."
    )
    parser.add_argument(
        '--output_plot', 
        type=str, 
        required=True, 
        help="Chemin vers le fichier .png de sortie."
    )
    
    args = parser.parse_args()
    process_taxonomy_data(args.input_txt, args.output_tsv, args.output_plot)

if __name__ == '__main__':
    main()
