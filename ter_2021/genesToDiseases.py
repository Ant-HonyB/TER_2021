import sys
import os
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from pyhpo.ontology import Ontology

#Les dictionnaires récupèrent les différents identifiants des symptomes référencés comme associés (Orpha/OMIM) ou non (NA).
diseases_OMIM = {}
diseases_Orpha = {}
diseases_NA = {}

def count_with_twins(hposet, current_disease):
    """

    :param hposet:
    :param current_disease:
    :return:
    """
    counter_OMIM = 0
    counter_Orpha = 0
    counter_NA = 0

    for q in range(len(hposet)):
        presence_OMIM = 0
        presence_Orpha = 0
        _hpo_omim = list(hposet[q].omim_diseases)
        _hpo_orpha = list(hposet[q].orpha_diseases)
        for diseases in range(len(_hpo_omim)):
            if current_disease.lower() in _hpo_omim[diseases].name.lower() and presence_OMIM == 0:
                diseases_OMIM.update({hposet[q].id: hposet[0].name})
                counter_OMIM += 1
                presence_OMIM = 1
        for diseases in range(len(_hpo_orpha)):
            if current_disease.lower() in _hpo_orpha[diseases].name.lower() and presence_Orpha == 0:
                diseases_Orpha.update({hposet[q].id: hposet[0].name})
                counter_Orpha += 1
                presence_Orpha = 1
        if presence_OMIM == 0 and presence_Orpha == 0:
            counter_NA += 1
            diseases_NA.update({hposet[q].id: hposet[q].name})

    compte_total = counter_OMIM + counter_Orpha + counter_NA
    duplicates = compte_total - len(hposet)
    return duplicates #on pourra rajouter les compteurs Orpha et Omim a mettre en legende dans les camemberts

    #print("compte :", compte_total, "et nbr_symptome :", len(hpo_terms))
    #print("doublons :", duplicates)
    #print("omim : ", counter_OMIM)
    #print("orpha : ", counter_Orpha)
    #print("NC : ", counter_NA)



def count_without_twins(hposet, current_disease):
    counter_OMIM = 0
    counter_Orpha = 0
    counter_NA = 0
    counter_found = 0

    for q in range(len(hposet)):
        presence_OMIM = 0
        presence_Orpha = 0
        _hpo_omim = list(hposet[q].omim_diseases)
        _hpo_orpha = list(hposet[q].orpha_diseases)
        for diseases in range(len(_hpo_omim)):
            if current_disease.lower() in _hpo_omim[diseases].name.lower() and presence_OMIM == 0:
                counter_found += 1
                counter_OMIM += 1
                presence_OMIM = 1
        for diseases in range(len(_hpo_orpha)):
            if current_disease.lower() in _hpo_orpha[diseases].name.lower() and presence_Orpha + presence_OMIM == 0 :
                counter_found += 1
                counter_Orpha += 1
                presence_Orpha = 1
        if presence_OMIM == 0 and presence_Orpha == 0:
            counter_NA += 1
    return counter_NA, counter_found

    #print("compte :", counter_found + counter_NA, "et nbr_symptome :", len(hpo_terms))
    #print("symptomes liés à la maladie : ", counter_found)
    #print("lié à orphanet :", counter_Orpha)
    #print("lié à Omim :", counter_OMIM)
    #print("présents dans les deux bases :", duplicates)
    #print("non-liés : ", counter_NA)




def convert_to_piechart(current_disease, output_dir, data_dir, counter_found, duplicates, counter_NA):
    labels = "associé à une BD", "associé aux deux BD", "Non-associé"
    sizes = [counter_found, duplicates, counter_NA]
    colors = ['yellowgreen', 'gold', 'lightskyblue']
    title = "Pas trouvé !"
    if "NCBI" in data_dir.absolute().as_posix():
        title = plt.title(f"{current_disease} (NCBI)")
    elif "CQR" in data_dir.absolute().as_posix():
        title = plt.title(f"{current_disease} (ConRur-Bio)")
    title.set_ha("center")
    plt.pie(sizes, labels=labels, colors=colors,
            autopct='%1.2f%%', startangle=90)
    plt.axis('equal')
    if "NCBI" in data_dir.absolute().as_posix():
        plt.savefig(f"{output_dir.joinpath(current_disease).absolute().as_posix()}_NCBI.png")
    if "NCBI" in data_dir.absolute().as_posix():
        plt.savefig(f"{output_dir.joinpath(current_disease).absolute().as_posix()}_CQR-BIO.png")



def main(indir, outdir):
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
            duplicates = count_with_twins(hposet, current_disease)
            counter_NA, counter_found = count_without_twins(hposet, current_disease)
            convert_to_piechart(current_disease, output_dir, data_dir, counter_found, duplicates, counter_NA)


if __name__ == '__main__':
    if len(sys.argv) == 3:
        main(sys.argv[1], sys.argv[2])
    else:
        print("usage :")
        print("python convertToCsv.py input_dir output_dir")
        print("input_dir  : directory containing files with lists of genes")
        print("output_dir : directory to store the csv files of related genes")
