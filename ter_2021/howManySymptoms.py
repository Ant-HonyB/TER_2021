import sys
import time
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from pyhpo import stats
from pyhpo.ontology import Ontology


def how_many_symptoms(hposet, current_disease):
    _useful_symptoms = []
    total_useful_symp = []
    useless_symptoms = []
    counter_symptoms = 0
    didnt_appear = 0

    for symptom in hposet:
        time.sleep(5)
        counter_symptoms += 1
        _useful_symptoms.append(symptom)
        disease_model = stats.EnrichmentModel('omim')
        disease_results = disease_model.enrichment('hypergeom', _useful_symptoms)
        presence_disease = 0
        print(f"la boucle avance {counter_symptoms}")

        for x in disease_results:
            if current_disease.lower() in x.get("item").name.lower():
                presence_disease = 1
                print("la maladie est présente")

        if presence_disease == 1:
            total_useful_symp.append(symptom)

        else:
            useless_symptoms.append(symptom)
            didnt_appear += 1

        if current_disease.lower() in disease_results[0].get("item").name.lower():
            print("la maladie est en tête")
            return total_useful_symp, useless_symptoms, counter_symptoms

        elif didnt_appear == len(hposet):
            print("Cette liste de termes ne reconnait jamais la maladie de départ.\n")
            return total_useful_symp, useless_symptoms, counter_symptoms
        elif counter_symptoms == len(hposet):
            print("La maladie de départ n'est jamais le résultat le plus pertinent.\n")
            return total_useful_symp, useless_symptoms, counter_symptoms



def to_barplot(counter_symptoms, useless_symptoms, total_useful_symp, current_disease,
               hposet, data_dir, output_dir):
    fig = plt.figure()
    ax = fig.add_axes([0,0,1,1])
    langs = ['signes cliniques\n nécéssaires', 'signes cliniques éronnés']
    res = [counter_symptoms, len(useless_symptoms)]
    #fig.legend(f"les symptomes utiles sont : {total_useful_symp}")
    title = "Pas trouvé !"
    if "NCBI" in data_dir.absolute().as_posix():
        title = plt.title(f"{current_disease} (NCBI)\n N = {len(hposet)}")
    elif "CQR" in data_dir.absolute().as_posix():
        title = plt.title(f"{current_disease} (ConQur-Bio)\n N = {len(hposet)}")
    #elif "bruit_test" in data_dir.absolute().as_posix():
       # title = plt.title(f"{current_disease} (test bruit)\n N = {compte_total}")
    title.set_ha("center")
    ax.bar(langs, res)
    if "NCBI" in data_dir.absolute().as_posix():
        plt.savefig(f"{output_dir.joinpath(current_disease).absolute().as_posix()}_NCBI.png")
    elif "CQR" in data_dir.absolute().as_posix():
        plt.savefig(f"{output_dir.joinpath(current_disease).absolute().as_posix()}_CQR-BIO.png")
    elif "bruit_test" in data_dir.absolute().as_posix():
        plt.savefig(f"{output_dir.joinpath(current_disease).absolute().as_posix()}_bruit-test.png")
    plt.figure().clear()


def main(indir, outdir):
    """

    :param indir:
    :param outdir:
    :return:
    """
    data_dir = Path(indir)
    output_dir = Path(outdir)
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
            to_barplot(counter_symptoms, useless_symptoms, total_useful_symp,
                       current_disease, hposet, data_dir, output_dir)



if __name__ == '__main__':
    if len(sys.argv) == 3:
        main(sys.argv[1], sys.argv[2])
    else:
        print("usage :")
        print("python howManySymptoms.py input_dir output_dir")
        print("input_dir  : directory containing files with lists of terms")
        print("output_dir : directory to store the csv files of related terms")
