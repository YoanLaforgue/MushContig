import argparse
import sys

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Convertit les scores Phred depuis un fichier .txt NanoStats en valeur base accuracy."
    )
    parser.add_argument(
        "--nanostats", 
        required=True, 
        help="Chemin vers le fichier .txt NanoStats issu de l'outil NanoPlot."
    )
    parser.add_argument(
        "--base_accuracy_report", 
        required=True, 
        help="Chemin vers le fichier .txt en sortie."
    )
    return parser.parse_args()

def calculate_phred_stats(q_score):
    """
    Calcule la probabilité d'erreur et la précision des bases à partir d'un score Phred Q NanoPlot.
    Q = -10 * log10(P)  =>  P = 10^(-Q/10)
    """
    try:
        q = float(q_score)
        error_probability = 10 ** (-q / 10)
        base_accuracy = (1 - error_probability) * 100
        return error_probability, base_accuracy
    except ValueError:
        return None, None

def extract_nanostats_data(filepath):

    extracted_data = {}
    
    try:
        with open(filepath, 'r', encoding='utf-8') as file:
            for line in file:
                line = line.strip()
                if not line:
                    continue
                
                
                if ':' in line:
                    key, value = line.split(':', 1)
                    key = key.strip()
                    
                    value = value.strip().replace(',', '')
                    
                    if key in ["Median read quality", "Mean read quality", "Total bases", "Number of reads"]:
                        extracted_data[key] = float(value)
                        
        return extracted_data
    except FileNotFoundError:
        print(f"Erreur : Le fichier {filepath} n'a pas été trouvé.")
        sys.exit(1)
    except Exception as e:
        print(f"Erreur lors de la lecture du fichier : {e}")
        sys.exit(1)

def generate_report(output_filepath, stats):
    """Génère un rapport texte lisible avec les calculs effectués."""
    try:
        with open(output_filepath, 'w', encoding='utf-8') as out:
            out.write("========================================================\n")
            out.write("                Base Accuracy report\n")
            out.write("========================================================\n\n")
            
            out.write("--- Statistiques générales du Run ---\n")
            if "Number of reads" in stats:
                out.write(f"Nombre total de lectures : {int(stats['Number of reads']):,}\n")
            if "Total bases" in stats:
                out.write(f"Nombre total de bases    : {int(stats['Total bases']):,}\n")
            out.write("\n")
            
            out.write("--- Résultats des calculs de qualité ---\n")
            out.write("Formules utilisées :\n")
            out.write("  - Probabilité d'erreur (P) = 10^(-Q/10)\n")
            out.write("  - Précision de base (%)    = (1 - P) * 100\n\n")

            if "Median read quality" in stats:
                q_median = stats["Median read quality"]
                p_err, accuracy = calculate_phred_stats(q_median)
                
                out.write("[MÉDIANE]\n")
                out.write(f"  Score Phred médian (Q)   : {q_median}\n")
                out.write(f"  Probabilité d'erreur (P) : {p_err:.6f} (soit 1 erreur toutes les {int(1/p_err)} bases)\n")
                out.write(f"  Précision de base        : {accuracy:.4f} %\n\n")
            else:
                out.write("[!] 'Median read quality' non trouvé dans le fichier d'entrée.\n\n")

            if "Mean read quality" in stats:
                q_mean = stats["Mean read quality"]
                p_err, accuracy = calculate_phred_stats(q_mean)
                
                out.write("[MOYENNE]\n")
                out.write(f"  Score Phred moyen (Q)    : {q_mean}\n")
                out.write(f"  Probabilité d'erreur (P) : {p_err:.6f} (soit 1 erreur toutes les {int(1/p_err)} bases)\n")
                out.write(f"  Précision de base        : {accuracy:.4f} %\n\n")

            out.write("========================================================\n")
        
        print(f"Succès : Le rapport a été généré avec succès dans -> '{output_filepath}'")
        
    except Exception as e:
        print(f"Erreur lors de l'écriture du fichier de sortie : {e}")
        sys.exit(1)

def main():
    args = parse_arguments()
    
    nanostats_data = extract_nanostats_data(args.nanostats)
    
    if not nanostats_data:
        print("Avertissement : Aucune donnée trouvée dans le fichier .txt NanoStats.")
        sys.exit(1)
        
    generate_report(args.base_accuracy_report, nanostats_data)

if __name__ == "__main__":
    main()
