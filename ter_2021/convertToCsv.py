import sys
from pathlib import Path
from gprofiler import GProfiler


def extract_list_from_line(line: str):
    # on enleve tous les [] et on reconstruit le string avec les crochets initial/final
    sansCrochets: str = line.strip().replace("[", "'").replace("]", "'")
    sansCrochets = f"[{sansCrochets[1:-1]}]"
    return eval(sansCrochets)


def main(indir, outdir):
    data_dir = Path(indir)  # dossier où se trouvent les listes de gènes provenant des maladies recherchées
    output_dir = Path(outdir)  # dossier où sont placés les résultats de g:profiler pour chaque liste

    if not data_dir.is_dir():
        raise ValueError(f"\nInput directory : {data_dir.as_posix()} is not a directory")

    if not output_dir.is_dir():
        raise ValueError(f"\nOutput directory : {output_dir.as_posix()} is not a directory")

    for fichier in list(data_dir.glob("*")):
        if not fichier.parts[-1].startswith(".") and fichier.is_file():
            with open(fichier.as_posix(), 'r') as f:
                print(fichier.absolute().as_posix())
                lines = f.readlines()

            # on appelle l'API g:profiler et on lui donne en entrée notre liste de Str en ciblant la base de données HPO.
            gp = GProfiler(return_dataframe=True)
            df = gp.profile(organism='hsapiens', query=extract_list_from_line(lines[0]), sources=['HP'])
            df.head()

            # conversion de nos données récupérées en csv avec le nom initial des données.
            # ..../output_dir/nom_de_la_maladie.csv
            df.to_csv(f"{output_dir.joinpath(fichier.parts[-1])}.csv")


if __name__ == '__main__':
    if len(sys.argv) == 3:
        main(sys.argv[1], sys.argv[2])
    else:
        print("usage :")
        print("python convertToCsv.py input_dir output_dir")
        print("input_dir  : directory containing files with lists of genes")
        print("output_dir : directory to store the csv files of related genes")
