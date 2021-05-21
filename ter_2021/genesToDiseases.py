import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from pyhpo import stats
from pyhpo.ontology import Ontology
from pyhpo.set import HPOSet, BasicHPOSet

data = pd.read_csv("/Users/Mopchi_44/Casier/Projets/TER_2021/outputs_conqurbio/FirstBatch/NCBI_Pheno/Autoimmune Lymphoproliferative Syndrome.csv")
current_disease = f"{output_dir.joinpath(fichier.parts[-1])}
_ = Ontology()
hpo_terms = data.name.to_list()
hposet = list([Ontology.match(q) for q in hpo_terms])

#Les dictionnaires récupèrent les différents identifiants des symptomes référencés comme associés (Orpha/OMIM) ou non (NA).
diseases_OMIM = {}
diseases_Orpha = {}
diseases_NA = {}
duplicates = 0     #compte les symptomes liés à la maladie de départ dans les deux BDs.
compte_total = 0   #total de la somme des counters
counter_OMIM = 0   #compte les symptomes associés au moins chez OMIM
counter_Orpha = 0  #compte les symptomes associés au moins chez Orphanet
counter_NA = 0     #compte les symptomes non-associés
counter_found = 0  #compte l'ensemble des symptomes liés à la maladie de départ


def count_with_twins():
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
            diseases_NA.update({hposet[q].id:hposet[q].name})

    compte_total = counter_OMIM + counter_Orpha + counter_NA
    duplicates = compte_total - len(hpo_terms)

    #print("compte :", compte_total, "et nbr_symptome :", len(hpo_terms))
    #print("doublons :", duplicates)
    #print("omim : ", counter_OMIM)
    #print("orpha : ", counter_Orpha)
    #print("NC : ", counter_NA)



def count_without_twins():
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

    #print("compte :", counter_found + counter_NA, "et nbr_symptome :", len(hpo_terms))
    #print("symptomes liés à la maladie : ", counter_found)
    #print("lié à orphanet :", counter_Orpha)
    #print("lié à Omim :", counter_OMIM)
    #print("présents dans les deux bases :", duplicates)
    #print("non-liés : ", counter_NA)




def convert_to_piechart():
    labels = "associé à une BD", "associé aux deux BD", "Non-associé"
    sizes = [counter_found, duplicates, counter_NA]
    colors = ['yellowgreen', 'gold', 'lightskyblue']
    title = plt.title("Autoimmune Lymphoproliferative Syndrome (NCBI)")

    title.set_ha("center")
    plt.pie(sizes, labels=labels, colors=colors,
            autopct='%1.2f%%', startangle=90)

    plt.axis('equal')
    plt.savefig("/Users/Mopchi_44/Casier/Projets/TER_2021/outputs_conqurbio/FirstBatch/les_camemberts/Autoimmune_Lymphoproliferative_Syndrome_NCBI.png")
    plt.show()



def main(indir, outdir):
    input_dir = Path(indir)  # dossier où se trouvent les listes de gènes provenant des maladies recherchées
    output_dir = Path(outdir)  # dossier où sont placés les résultats de g:profiler pour chaque liste



if __name__ == '__main__':
    if len(sys.argv) == 3:
        main(sys.argv[1], sys.argv[2])
    else:
        print("usage :")
        print("python convertToCsv.py input_dir output_dir")
        print("input_dir  : directory containing files with lists of genes")
        print("output_dir : directory to store the csv files of related genes")
