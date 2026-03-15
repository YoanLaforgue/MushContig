import argparse

def rename_fasta_data(input_fasta, input_txt, output_fasta):
    """
    Modifie un fichier .fasta en renommant les en-têtes des séquences
    à partir des informations contenues dans un fichier .txt.

    Args:
        input_fasta (str): Chemin vers le fichier .fasta d'entrée.
        input_txt (str): Chemin vers le fichier .txt d'entrée.
        output_fasta (str): Chemin vers le fichier .fasta de sortie.
    """
    taxonomie_dict = {}
    try:
        with open(input_txt, 'r', encoding='utf-8') as f_txt:
            for ligne in f_txt:
                if not ligne.strip():
                    continue
                
                parts = ligne.strip().split()
                if len(parts) < 2:
                    continue
                
                identifiant = parts[0]
                taxonomie_brute = parts[1]
                
                taxo_elements = taxonomie_brute.split(';')
                if len(taxo_elements) >= 2:
                    info_a_recuperer = taxo_elements[-2:]
                    taxonomie_formatee = "/".join(info_a_recuperer)
                    taxonomie_dict[identifiant] = taxonomie_formatee

    except FileNotFoundError:
        print(f"Erreur : Le fichier texte '{input_txt}' n'a pas été trouvé.")
        return
    except Exception as e:
        print(f"Une erreur est survenue lors de la lecture du fichier .txt : {e}")
        return

    try:
        with open(input_fasta, 'r', encoding='utf-8') as f_in, \
             open(output_fasta, 'w', encoding='utf-8') as f_out:
            
            for ligne in f_in:
                if ligne.startswith('>'):
                    titre_initial = ligne.strip()[1:]
                    if titre_initial in taxonomie_dict:
                        taxonomie_a_ajouter = taxonomie_dict[titre_initial]
                        nouveau_titre = f">{titre_initial}_{taxonomie_a_ajouter}\n"
                        f_out.write(nouveau_titre)
                    else:
                        f_out.write(ligne)
                else:
                    f_out.write(ligne)

        print(f"Le fichier .fasta modifié a été sauvegardé sous : {output_fasta}")

    except FileNotFoundError:
        print(f"Erreur : Le fichier .fasta d'entrée '{input_fasta}' n'a pas été trouvé.")
    except Exception as e:
        print(f"Une erreur est survenue lors du traitement du fichier .fasta : {e}")


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
        description="Script pour renommer les séquences d'un fichier .fasta à partir d'informations taxonomiques contenues dans un fichier .txt issu de vsearch."
    )
    
    parser.add_argument(
        '--input_fasta', 
        type=str, 
        required=True, 
        help="Chemin vers le fichier .fasta d'entrée."
    )
    parser.add_argument(
        '--input_txt', 
        type=str, 
        required=True, 
        help="Chemin vers le fichier .txt d'entrée."
    )
    parser.add_argument(
        '--output_fasta', 
        type=str, 
        required=True, 
        help="Chemin vers le fichier .fasta renommé en sortie."
    )
    
   args = parser.parse_args()
   rename_fasta_data(args.input_fasta, args.input_txt, args.output_fasta)
