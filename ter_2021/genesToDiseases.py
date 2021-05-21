import pandas as pd
import matplotlib.pyplot as plt
from pyhpo import stats
from pyhpo.ontology import Ontology
from pyhpo.set import HPOSet, BasicHPOSet

data = pd.read_csv("/Users/Mopchi_44/Casier/Projets/TER_2021/outputs_conqurbio/FirstBatch/NCBI_Pheno/Autoimmune Lymphoproliferative Syndrome.csv")
current_disease = "Autoimmune Lymphoproliferative Syndrome"
_ = Ontology()
hpo_terms = data.name.to_list()
hposet = list([Ontology.match(q) for q in hpo_terms])

#Les dictionnaires récupèrent les différents identifiants des symptomes référencés comme associés (Orpha/OMIM) ou non (NA).
diseases_OMIM = {}
diseases_Orpha = {}
diseases_NA = {}
doublons = 0        #compte les symptomes liés à la maladie de départ dans les deux BDs.
compte_total = 0    #total de la somme des compteurs
compteur_OMIM = 0   #compte les symptomes associés au moins chez OMIM
compteur_Orpha = 0  #compte les symptomes associés au moins chez Orphanet
compteur_NA = 0     #compte les symptomes non-associés
compteur_found = 0  #compte l'ensemble des symptomes liés à la maladie de départ


def count_with_twins():
    compteur_OMIM = 0
    compteur_Orpha = 0
    compteur_NA = 0

    for q in range(len(hposet)):
        presence_OMIM = 0
        presence_Orpha = 0
        _hpo_omim = list(hposet[q].omim_diseases)
        _hpo_orpha = list(hposet[q].orpha_diseases)
        for diseases in range(len(_hpo_omim)):
            if current_disease.lower() in _hpo_omim[diseases].name.lower() and presence_OMIM == 0:
                diseases_OMIM.update({hposet[q].id: hposet[0].name})
                compteur_OMIM += 1
                presence_OMIM = 1
        for diseases in range(len(_hpo_orpha)):
            if current_disease.lower() in _hpo_orpha[diseases].name.lower() and presence_Orpha == 0:
                diseases_Orpha.update({hposet[q].id: hposet[0].name})
                compteur_Orpha += 1
                presence_Orpha = 1
        if presence_OMIM == 0 and presence_Orpha == 0:
            compteur_NA += 1
            diseases_NA.update({hposet[q].id:hposet[q].name})

    compte_total = compteur_OMIM + compteur_Orpha + compteur_NA
    doublons = compte_total - len(hpo_terms)

    #print("compte :", compte_total, "et nbr_symptome :", len(hpo_terms))
    #print("doublons :", doublons)
    #print("omim : ", compteur_OMIM)
    #print("orpha : ", compteur_Orpha)
    #print("NC : ", compteur_NA)



def count_without_twins():
    compteur_OMIM = 0
    compteur_Orpha = 0
    compteur_NA = 0
    compteur_found = 0

    for q in range(len(hposet)):
        presence_OMIM = 0
        presence_Orpha = 0
        _hpo_omim = list(hposet[q].omim_diseases)
        _hpo_orpha = list(hposet[q].orpha_diseases)
        for diseases in range(len(_hpo_omim)):
            if current_disease.lower() in _hpo_omim[diseases].name.lower() and presence_OMIM == 0:
                compteur_found += 1
                compteur_OMIM += 1
                presence_OMIM = 1
        for diseases in range(len(_hpo_orpha)):
            if current_disease.lower() in _hpo_orpha[diseases].name.lower() and presence_Orpha + presence_OMIM == 0 :
                compteur_found += 1
                compteur_Orpha += 1
                presence_Orpha = 1
        if presence_OMIM == 0 and presence_Orpha == 0:
            compteur_NA += 1

    #print("compte :", compteur_found + compteur_NA, "et nbr_symptome :", len(hpo_terms))
    #print("symptomes liés à la maladie : ", compteur_found)
    #print("lié à orphanet :", compteur_Orpha)
    #print("lié à Omim :", compteur_OMIM)
    #print("présents dans les deux bases :", doublons)
    #print("non-liés : ", compteur_NA)





labels = "associé à une BD", "associé aux deux BD", "Non-associé"
sizes = [compteur_found, doublons, compteur_NA]
colors = ['yellowgreen', 'gold', 'lightskyblue']
title = plt.title("Autoimmune Lymphoproliferative Syndrome (NCBI)")

title.set_ha("center")
plt.pie(sizes, labels=labels, colors=colors,
        autopct='%1.2f%%', startangle=90)

plt.axis('equal')
plt.savefig("/Users/Mopchi_44/Casier/Projets/TER_2021/outputs_conqurbio/FirstBatch/les_camemberts/Autoimmune_Lymphoproliferative_Syndrome_NCBI.png")
plt.show()









def main(indir, outdir):
    # something


if __name__ == '__main__':
    if len(sys.argv) == 3:
        main(sys.argv[1], sys.argv[2])
    else:
        print("usage :")
        print("python convertToCsv.py input_dir output_dir")
        print("input_dir  : directory containing files with lists of genes")
        print("output_dir : directory to store the csv files of related genes")
