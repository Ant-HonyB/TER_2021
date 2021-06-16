import sys
from pathlib import Path
import pandas as pd
from pyhpo import stats
from pyhpo.ontology import Ontology


def how_many_symptoms(hposet, current_disease):
    _useful_symptoms = hposet
    counter_symptoms = 0

    for x in hposet:
        counter_symptoms += 1
        disease_model = stats.EnrichmentModel('omim')
        disease_results = disease_model.enrichment('hypergeom', _useful_symptoms)

        if not current_disease.lower() in disease_results[0].get("item").name.lower():
            total_useful_symp = _useful_symptoms
            return total_useful_symp, counter_symptoms

        _useful_symptoms = _useful_symptoms[-1]



def main(indir, outdir):
    """

    :param indir:
    :param outdir:
    :return:
    """
    data_dir = Path(indir)  # dossier où se trouvent les listes de gènes provenant des maladies recherchées
    output_dir = Path(outdir)  # dossier où sont placés les résultats de g:profiler pour chaque liste
    if not data_dir.is_dir():
        raise ValueError(f"\nInput directory : {data_dir.as_posix()} is not a directory")

    if not output_dir.is_dir():
        raise ValueError(f"\nOutput directory : {output_dir.as_posix()} is not a directory")

    for fichier in list(data_dir.glob("*.csv")):
        if not fichier.parts[-1].startswith(".") and fichier.is_file():
            data = pd.read_csv(fichier.absolute().as_posix())
            current_disease = fichier.name.replace(".csv", "")
            _ = Ontology()
            hpo_terms = data.name.to_list()
            hposet = [Ontology.match(q) for q in hpo_terms]
            total_useful_symp, useless_symptoms, counter_symptoms = how_many_symptoms(hposet, current_disease)

if __name__ == '__main__':
    if len(sys.argv) == 3:
        main(sys.argv[1], sys.argv[2])
    else:
        print("usage :")
        print("python howManySymptoms.py input_dir output_dir")
        print("input_dir  : directory containing files with lists of terms")
        print("output_dir : directory to store the csv files of related terms")
