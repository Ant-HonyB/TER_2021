import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from pyhpo.ontology import Ontology

#Les dictionnaires récupèrent les différents identifiants des symptomes référencés comme associés (Orpha/OMIM) ou non (NA).
#Ils ne sont pas exploité par le programme pour le moment.
diseases_OMIM = {}
diseases_Orpha = {}
diseases_NA = {}


def count_with_twins(hposet, current_disease):
    """
    Interroge l'API HPO et demande, pour chaque signe clinique, s'il est affilié à la maladie de départ.
    Fait le compte total des signes cliniques affilié à la base de données OMIM et/ou Orphanet,
    et celui des signes n'apparaissant dans aucune base de données.
    :param hposet: ensemble des signes cliniques du fichier initial converti en des termes référencés par HPO
    :param current_disease: maladie de départ recherchée sur Entrez Gene de NCBI ou sur ConQuR-Bio.
    :return: le compte des signes cliniques respectivement référencés soit par OMIM, par Orphanet,
    par les deux ou par aucun.
    """

    counter_OMIM = 0
    counter_Orpha = 0
    counter_NA = 0
    duplicates = 0

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
            if current_disease.lower() in _hpo_orpha[diseases].name.lower() and presence_Orpha == 0 and presence_OMIM == 1:
                duplicates += 1
                if not counter_OMIM == 0:
                    counter_OMIM -= 1
                diseases_Orpha.update({hposet[q].id: hposet[0].name})
                presence_Orpha = 1
            elif current_disease.lower() in _hpo_orpha[diseases].name.lower() and presence_Orpha + presence_OMIM == 0:
                counter_Orpha += 1
                diseases_Orpha.update({hposet[q].id: hposet[0].name})
                presence_Orpha = 1
        if presence_OMIM + presence_Orpha == 0:
            counter_NA += 1
            diseases_NA.update({hposet[q].id: hposet[q].name})

    compte_total = len(hposet)
    return duplicates, counter_OMIM, counter_Orpha, counter_NA, compte_total


def convert_to_graphs(current_disease, output_dir, data_dir, counter_OMIM, duplicates,
                        counter_Orpha, counter_NA, compte_total):
    """
    Exploite les différents comptes obtenus par count_with_twins() sous forme de camemberts
    et de diagramme de Venn.
    :param current_disease: maladie de départ recherchée sur Entrez Gene de NCBI ou sur ConQuR-Bio.
    :param output_dir: dossier où sont déposés les différents graphiques obtenus sur la base de l'ensemble
    des fichiers provenant du dossier initial data_dir.
    :param data_dir: dossier où se trouve les fichiers contenant les signes cliniques récupérés à partir de g:profiler
    sur la base des recherches de maladie de départ entrées sur Entrez Gene de NCBI ou sur ConQuR-Bio.
    :param counter_OMIM: compte des signes cliniques affiliées à la maladie de départ chez OMIM.
    :param duplicates: compte des signes cliniques affiliées à la maladie de départ chez OMIM ET Orphanet.
    :param counter_Orpha: compte des signes cliniques affiliées à la maladie de départ chez Orphanet.
    :param counter_NA: compte des signes cliniques affiliées à la maladie de départ chez aucune base de données.
    :return: créer des fichiers images contenant respectivement les camemberts et les diagrammes de Venn.
    """

    labels = "", "", "", ""
    sizes = [counter_OMIM, duplicates, counter_Orpha, counter_NA]
    colors = ['khaki', 'springgreen', "turquoise", 'gainsboro']
    title = "Pas trouvé !"
    if "NCBI" in data_dir.absolute().as_posix():
        title = plt.title(f"{current_disease} (NCBI)\n N = {compte_total}")
    elif "CQR" in data_dir.absolute().as_posix():
        title = plt.title(f"{current_disease} (ConQur-Bio)\n N = {compte_total}")
    elif "bruit_test" in data_dir.absolute().as_posix():
        title = plt.title(f"{current_disease} (test bruit)\n N = {compte_total}")
    title.set_ha("center")
    plt.pie(sizes, labels=labels, colors=colors, pctdistance=0.5, labeldistance=1.1,
            autopct='%1.1f%%', startangle=90)
    plt.legend([f"OMIM : {counter_OMIM}", f"Orphanet + OMIM : {duplicates}",
                f"Orphanet : {counter_Orpha}", f"Absent : {counter_NA}"], loc='lower right')
    plt.axis('equal')
    if "NCBI" in data_dir.absolute().as_posix():
        plt.savefig(f"{output_dir.joinpath(current_disease).absolute().as_posix()}_NCBI.png")
    elif "CQR" in data_dir.absolute().as_posix():
        plt.savefig(f"{output_dir.joinpath(current_disease).absolute().as_posix()}_CQR-BIO.png")
    elif "bruit_test" in data_dir.absolute().as_posix():
        plt.savefig(f"{output_dir.joinpath(current_disease).absolute().as_posix()}_bruit-test.png")
    plt.figure().clear()

    v = venn3(subsets=(counter_OMIM+1, counter_Orpha+1, duplicates+1, counter_NA+1, 0, 0, 0), set_labels=("", "", ""),
              set_colors=("khaki", "turquoise", "gainsboro"),
              alpha=0.5)
    v.get_label_by_id('100').set_text(f"{counter_OMIM}")
    v.get_label_by_id('010').set_text(f"{counter_Orpha}")
    v.get_label_by_id('001').set_text(f"{counter_NA}")
    v.get_label_by_id('110').set_text(f"{duplicates}")
    v.get_patch_by_id('110').set_color('springgreen')
    plt.title("on a un problème huston")
    plt.legend([f"OMIM : {counter_OMIM}", f"Orphanet : {counter_Orpha}", f"Orphanet + OMIM : {duplicates}",
                f"Absent : {counter_NA}"], loc='lower right')
    if "NCBI" in data_dir.absolute().as_posix():
        plt.title(f"{current_disease} (NCBI)\n N = {compte_total}")
    elif "CQR" in data_dir.absolute().as_posix():
        plt.title(f"{current_disease} (ConQur-Bio)\n N = {compte_total}")
    elif "bruit_test" in data_dir.absolute().as_posix():
        plt.title(f"{current_disease} (test bruit)\n N = {compte_total}")
    if "NCBI" in data_dir.absolute().as_posix():
        plt.savefig(f"{output_dir.joinpath(current_disease).absolute().as_posix()}_Venn_NCBI.png")
    elif "CQR" in data_dir.absolute().as_posix():
        plt.savefig(f"{output_dir.joinpath(current_disease).absolute().as_posix()}_Venn_CQR-BIO.png")
    elif "bruit_test" in data_dir.absolute().as_posix():
        plt.savefig(f"{output_dir.joinpath(current_disease).absolute().as_posix()}_Venn_bruit-test.png")
    plt.figure().clear()







def main(indir, outdir):
    """
    Le Main s'occupe d'automatiser le traitement de chaque fichier du dossier de départ. A partir des listes
    de signes cliniques de chaque maladie de départ, il exploite les deux fonctions du programme pour déposer
    dans le dossier final les résultats des comptes sous forme de camemberts et de diagramme de Venn.
    :param indir: chemin d'accès au dossier initial où se trouve les fichiers exploités
    :param outdir: chemin d'accès au dossier final où sont déposés les résultats du programme.
    :return: les résultats composés de camemberts et de diagramme de Venn mettant en évidence
    les proportions d'affiliation entre signes cliniques et maladie de départ.
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
            duplicates, counter_OMIM, counter_Orpha, counter_NA, compte_total = count_with_twins(hposet, current_disease)
            convert_to_graphs(current_disease, output_dir, data_dir, counter_OMIM,
                                duplicates, counter_Orpha, counter_NA, compte_total)


if __name__ == '__main__':
    if len(sys.argv) == 3:
        main(sys.argv[1], sys.argv[2])
    else:
        print("usage :")
        print("python convertToCsv.py input_dir output_dir")
        print("input_dir  : directory containing files with lists of genes")
        print("output_dir : directory to store the csv files of related genes")
