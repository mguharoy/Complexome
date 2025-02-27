from typing import Optional
import io
import csv
import math
import sys
from dataclasses import dataclass
from urllib.error import URLError
from urllib.request import urlopen
from pathlib import Path
from datetime import date
import time
import pandas as pd
import numpy as np

from matplotlib.axes import Axes
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import textwrap

ONLY_REGULATED_SUBUNITS = True

ComplexT = dict[str, list[str]]
FileName = str
FileContents = bytes


@dataclass(frozen=True)
class SubunitInfo:
    complex_id: str
    name: str
    subunit: str
    log2fc: float
    apvalue: float


@dataclass(frozen=True)
class Complexome:
    file: Path
    complexes: ComplexT
    complex_names: dict[str, str]
    complex_GO_terms: ComplexT
    proteomics_data: dict[str, tuple[float, float]]
    all_perturbed_complexes: list[SubunitInfo]
    LOG2FC_THRESHOLD: float = 2.0
    ADJP_THRESHOLD: float = 0.05
    topN_GOterms_to_plot: int = 10


def _urlretrieve_with_retries(
    url: str, filename: str, retries: int = 3, delay: float = 1.5
) -> None:
    for attempt in range(retries):
        try:
            with urlopen(url) as response, open(filename, mode="wb") as out_file:
                out_file.write(response.read())
            return
        except URLError:
            if attempt < retries - 1:
                time.sleep(delay)
            else:
                raise


def setup(
    proteomics_data: dict[
        FileName, FileContents
    ],  # This is a mapping that we get from Colab
    LOG2FC_THRESHOLD: float,
    ADJP_THRESHOLD: float,
    topN_GOterms_to_plot: int,
    organism_taxon_id: str = "9606",
) -> Complexome:
    today = date.today()
    ComplexomeSavedFile = Path(
        organism_taxon_id + ".tsv"
    )  # Downloaded complexome file in ComplexTAB format.
    ComplexomeFile = Path(
        organism_taxon_id + "Complexome_" + str(today) + ".tsv"
    )  # Rename complexome file with a date stamp for future reference.

    if not ComplexomeFile.exists():
        _urlretrieve_with_retries(
            f"https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/{ComplexomeSavedFile}",
            filename=str(ComplexomeFile),
        )

    complexes, complex_names, complex_GO_terms = parse_complexome_data(
        organism_taxon_id, str(ComplexomeFile)
    )
    if len(proteomics_data) > 1:
        raise RuntimeError("Cannot accept multiple input files")
    parsed_proteomics_data = parse_user_proteomics_data(
        list(proteomics_data.values())[0]
    )

    return Complexome(
        file=ComplexomeFile,
        complexes=complexes,
        complex_names=complex_names,
        complex_GO_terms=complex_GO_terms,
        LOG2FC_THRESHOLD=LOG2FC_THRESHOLD,
        ADJP_THRESHOLD=ADJP_THRESHOLD,
        topN_GOterms_to_plot=topN_GOterms_to_plot,
        proteomics_data=parsed_proteomics_data,
        all_perturbed_complexes=identify_perturbed_complexes(
            complexes,
            complex_names,
            parsed_proteomics_data,
            LOG2FC_THRESHOLD,
            ADJP_THRESHOLD,
        ),
    )


def addComplexParticipants(participantsList: str, members=None) -> list[str]:
    if members is None:
        members = []

    for participant in participantsList.split("|"):
        participantId = participant.split("(")[0]
        if "[" in participantId:
            # This is in the case of molecule sets (paralogs that cannot be distinguished in this context).
            participant_ids = (
                str(participantId).replace("[", "").replace("]", "").split(",")
            )
            for paralog in participant_ids:
                if paralog not in members:
                    members.append(paralog)
        else:
            if participantId not in members:
                members.append(participantId)

    return members


def ComplexSubunitNumbersDistribution(complexesDict):
    numSubunitsPerComplex = {}
    for key, value in complexesDict.items():
        if len(value) not in numSubunitsPerComplex:
            numSubunitsPerComplex[len(value)] = [key]
        else:
            numSubunitsPerComplex[len(value)].append(key)

    numSubunitsPerComplexDist = {
        key: len(value) for key, value in numSubunitsPerComplex.items()
    }

    return numSubunitsPerComplexDist


def write_to_csv(moleculeList, outputFileName):
    with open(outputFileName, "w") as outputFile:
        for molecule in moleculeList:
            outputFile.write(molecule + "\n")


def parse_complexome_data(
    organism_taxon_id: str, complexome_file: str
) -> tuple[ComplexT, dict[str, str], ComplexT]:
    # Open and parse the Complex Portal data file and collect the per complex list of participants.
    Complexes = {}
    ComplexNames = {}
    ComplexesGOTerms = {}
    with open(complexome_file) as cmplx:
        isHeader = cmplx.read(1) == "#"

    with open(complexome_file) as fd:
        rd = csv.reader(fd, delimiter="\t")
        for row in rd:
            if isHeader:
                isHeader = False
                continue

            # Each row corresponds to a specific, annotated complex.
            complexID = row[0]
            complexName = row[1]
            complexTaxonID = row[3]
            complexParticipants = row[4]
            complexExtendedParticipants = row[-1]
            complexAnnotatedGOTerms = row[7]

            # Make sure that the organism taxon id is matching.
            if complexTaxonID != organism_taxon_id:
                print(row)
                print("Organism taxon id does not match with the expected organism!")
                sys.exit()

            # Extract information about the participants of the complex. Participants can include proteins (UniProtKB), chemical entities (ChEBI), RNA (RNAcentral) and complexes (Complex Portal).
            complexMembers = addComplexParticipants(complexParticipants)

            # Check if one or more participants are themselves complexes. In that case, the expanded list of protein members are contained in the Expanded participant list (last) column.
            if "CPX-" in complexParticipants:
                complexMembers = addComplexParticipants(
                    complexExtendedParticipants, complexMembers
                )

            # Add information about this complex and its participant molecules into the global dictionary of complexes.
            if complexID not in Complexes:
                Complexes[complexID] = complexMembers
                ComplexNames[complexID] = complexName
            else:
                print("Complex ID already exists:", complexID)
                print("Stopping for checks!!")
                sys.exit()

            # Extract information about the annotated GO terms of the complex.
            complexGOTerms = complexAnnotatedGOTerms.split("|")

            # Add information about this complex and its associated GO terms into the global dictionary of GO terms.
            if complexID not in ComplexesGOTerms:
                ComplexesGOTerms[complexID] = complexGOTerms
            else:
                print("Complex ID already exists:", complexID)
                print("Stopping for checks!!")
                sys.exit()
    return (Complexes, ComplexNames, ComplexesGOTerms)


def summary_statistics(complexome: Complexome) -> None:
    print("Total number of annotated complexes:", len(complexome.complexes))
    print("Total number of complex names:", len(complexome.complex_names))
    print(
        "Total number of complexes with associated GO terms:",
        len(complexome.complex_GO_terms),
    )


def unique_identities(complexome: Complexome) -> tuple[list[str], list[str]]:
    # Iterate over the list of complexes and store the unique proteins, metabolites and RNA molecules contained in that complexome.
    uniqueProteins = []
    uniqueMetabolites = []
    for complex in complexome.complexes.values():
        for subunit in complex:
            if "CPX-" in subunit:
                continue
            elif "URS" in subunit:
                continue
            elif "CHEBI:" in subunit:
                if subunit not in uniqueMetabolites:
                    uniqueMetabolites.append(subunit)
            else:
                if subunit not in uniqueProteins:
                    uniqueProteins.append(subunit)

    # print("Total number of unique proteins:", len(uniqueProteins))
    # print("Total number of unique metabolites:", len(uniqueMetabolites))

    # write_to_csv(
    #    uniqueMetabolites, organism_taxon_id + "_metabolites_" + str(today) + ".csv"
    # )
    # write_to_csv(uniqueProteins, organism_taxon_id + "_proteins_" + str(today) + ".csv")
    return uniqueProteins, uniqueMetabolites


def plot(complexome: Complexome, axis: Optional[Axes] = None) -> Axes:
    if axis is None:
        axis = plt.subplot()

    numAllSubunitsPerComplexDist = ComplexSubunitNumbersDistribution(
        complexome.complexes
    )

    axis.bar(
        list(numAllSubunitsPerComplexDist.keys()),
        numAllSubunitsPerComplexDist.values(),
        color="g",
    )
    axis.set_xlabel("Number of subunits", fontsize=16)
    axis.set_ylabel("Number of complexes", fontsize=16)
    axis.set_xticks(axis.get_xticks(), size=14)
    axis.set_yticks(axis.get_yticks(), size=14)
    axis.set_title("Subunit distribution (proteins, metabolites, RNA)", fontsize=18)

    return axis


def proteins_only(complexome: Complexome, axis: Optional[Axes] = None) -> Axes:
    if axis is None:
        axis = plt.subplot()
    # Iterate over the list of complexes and store the protein subunits only per complex.
    proteinSubunitsPerComplex = {}
    for complex_id, complex in complexome.complexes.items():
        proteinSubunits = []
        for subunit in complex:
            if "CPX-" in subunit:
                continue
            elif "URS" in subunit:
                continue
            elif "CHEBI:" in subunit:
                continue
            else:
                if subunit not in proteinSubunits:
                    proteinSubunits.append(subunit)

        proteinSubunitsPerComplex[complex_id] = proteinSubunits

    numProteinSubunitsPerComplexDist = ComplexSubunitNumbersDistribution(
        proteinSubunitsPerComplex
    )

    axis.bar(
        list(numProteinSubunitsPerComplexDist.keys()),
        numProteinSubunitsPerComplexDist.values(),
        color="g",
    )
    axis.set_xlabel("Number of protein subunits", fontsize=16)
    axis.set_ylabel("Number of complexes", fontsize=16)
    axis.set_title("Subunit distribution (proteins only)", fontsize=18)
    return axis


def shared_protein_subunits(
    complexome: Complexome, axis: Optional[Axes] = None
) -> Axes:
    if axis is None:
        axis = plt.subplot()
    # Compute the distibution of shared protein subunits among the different complexes.
    uniqueProteins, uniqueMetabolites = unique_identities(complexome)
    proteinSubunitsPerComplex = {}
    for complex_id, complex in complexome.complexes.items():
        proteinSubunits = []
        for subunit in complex:
            if "CPX-" in subunit:
                continue
            elif "URS" in subunit:
                continue
            elif "CHEBI:" in subunit:
                continue
            else:
                if subunit not in proteinSubunits:
                    proteinSubunits.append(subunit)

        proteinSubunitsPerComplex[complex_id] = proteinSubunits
    proteinsInNumComplexes = {}
    for protein in uniqueProteins:
        for proteinSubunits in proteinSubunitsPerComplex.values():
            if protein in proteinSubunits:
                if protein not in proteinsInNumComplexes:
                    proteinsInNumComplexes[protein] = 1
                else:
                    proteinsInNumComplexes[protein] += 1

    sharedSubunitsPerComplex = {}
    for key, value in proteinsInNumComplexes.items():
        # Useful to collate all values >10 into the last bar on the plot (see below).
        if value > 10:
            value = 10
        if value not in sharedSubunitsPerComplex:
            sharedSubunitsPerComplex[value] = 1
        else:
            sharedSubunitsPerComplex[value] += 1

    axis.bar(
        list(sharedSubunitsPerComplex.keys()),
        list(sharedSubunitsPerComplex.values()),
        color="g",
    )
    axis.set_xlabel("Number of complexes", fontsize=16)
    axis.set_ylabel("Number of proteins", fontsize=16)
    axis.set_xticks(
        range(1, 11), ["1", "2", "3", "4", "5", "6", "7", "8", "9", ">10"], size=14
    )
    axis.set_title("Distribution of shared protein subunits", fontsize=18)

    return axis


def proteomics_coverage_of_complexome(
    complexome: Complexome, axis: Optional[Axes] = None
) -> Axes:
    if axis is None:
        axis = plt.subplot()
    proteomicsCoveragePerComplex: dict[str, float] = {}
    for complex_id, complex in complexome.complexes.items():
        numSubunits = 0
        measuredSubunits = 0
        for subunit in complex:
            if "CPX-" in subunit:
                continue
            elif "URS" in subunit:
                continue
            elif "CHEBI:" in subunit:
                continue
            else:
                numSubunits += 1
                canonicalUniProtID = subunit[:6]
                if canonicalUniProtID in complexome.proteomics_data:
                    measuredSubunits += 1
        proteomicsCoveragePerComplex[complex_id] = measuredSubunits / numSubunits

    _min = min(proteomicsCoveragePerComplex.values())
    _max = max(proteomicsCoveragePerComplex.values())
    if _min == _max:
        _max = _max + sys.float_info.epsilon
        bins = 1
        delta = 1.0
    else:
        bins = min(
            len(proteomicsCoveragePerComplex),
            max(10, len(proteomicsCoveragePerComplex) // 160),
        )
        delta = (_max - _min) / bins

    histogram: list[int] = []
    xs: list[float] = []
    _base = _min
    i = 1
    while _base < _max:
        if i == bins:
            top = _max + sys.float_info.epsilon
        else:
            top = _base + delta
        histogram.append(
            sum(
                1
                for _ in filter(
                    lambda x: _base <= x < top, proteomicsCoveragePerComplex.values()
                )
            )
        )
        xs.append(_base)
        if i == bins:
            _base = _max
        else:
            _base += delta
        i += 1

    assert sum(histogram) == len(proteomicsCoveragePerComplex), (
        f"{sum(histogram)=} /= {len(proteomicsCoveragePerComplex)}"
    )
    axis.bar(x=xs, height=histogram, width=delta, align="edge", color="g")
    axis.set_xbound(lower=0.0, upper=1.0)
    axis.set_xlabel("Proteomics Coverage")
    axis.set_ylabel("Count")
    return axis


def plot_volcano2(
    complexome: Complexome, axis: Optional[Axes] = None, marker_size: Optional[int] = 20
) -> Axes:
    if axis is None:
        axis = plt.subplot()

    def colour(log2fc: float, pval: float) -> str:
        if pval > -math.log10(complexome.ADJP_THRESHOLD):
            if log2fc > complexome.LOG2FC_THRESHOLD:
                return "tab:red"
            elif log2fc < -complexome.LOG2FC_THRESHOLD:
                return "tab:blue"
            else:
                return "tab:gray"
        else:
            return "tab:gray"

    transformed_data = [
        (log2FC, -math.log10(adj_pval))
        for (log2FC, adj_pval) in complexome.proteomics_data.values()
    ]
    colours = [colour(*datum) for datum in transformed_data]
    (xvals, yvals) = zip(*transformed_data)
    axis.scatter(xvals, yvals, s=marker_size, c=colours)
    axis.set_title("Volcano plot")
    axis.set_xlabel("log2 (FC)")
    axis.set_ylabel("-log10 (adjPval)")
    return axis


def plot_volcano(complexome: Complexome, axis: Optional[Axes] = None) -> Axes:
    dfProteomics = pd.DataFrame.from_dict(
        complexome.proteomics_data, orient="index", columns=["log2FC", "adjPval"]
    )
    dfProteomics.index.names = ["UniProt ID"]
    dfProteomics = dfProteomics.reset_index()
    dfProteomics[["log2FC", "adjPval"]] = dfProteomics[["log2FC", "adjPval"]].apply(
        pd.to_numeric
    )
    p_values = dfProteomics["adjPval"].tolist()
    transformed_pvals = []  # Volcano plot y-axis needs -log10(adjPval)
    for value in p_values:
        transformed_pvals.append(-1 * np.log10(value))
    dfProteomics["qvals"] = transformed_pvals

    fig = go.Figure()
    fig.update_layout(
        title="Volcano plot",
        xaxis_title="log2 FC",
        yaxis_title="-log10 (adjPval)",
        paper_bgcolor="white",
        plot_bgcolor="white",
    )

    # We can play around with this. For example, when Complex Viewer is used to highlight a specific complex,
    # this part can be used to show the complete Volcano plot with the specific subunits of the selected complex highlighted.
    colors = []
    foldchanges = dfProteomics["log2FC"].tolist()
    for i in range(0, len(foldchanges)):
        if transformed_pvals[i] > -1 * np.log10(complexome.ADJP_THRESHOLD):
            if foldchanges[i] > complexome.LOG2FC_THRESHOLD:
                colors.append("red")
            elif foldchanges[i] < -1.0 * complexome.LOG2FC_THRESHOLD:
                colors.append("blue")
            else:
                colors.append("grey")
        else:
            colors.append("grey")

    fig.add_trace(
        go.Scattergl(
            x=foldchanges,
            y=transformed_pvals,
            mode="markers",
            text=dfProteomics["UniProt ID"].tolist(),
            hovertemplate="%{text}: %{x}<br>",  # put the UniProtID: log2FC, adjPval, later use the gene name.
            marker={
                "color": colors,
                "size": 6,
            },
        )
    )

    fig.show()


def wrap_labels(ax, width, topN_GOterms_to_plot, break_long_words=False):
    """
    Word wrap long GO term names (used as axis labels).
    Based on https://medium.com/dunder-data/automatically-wrap-graph-labels-in-matplotlib-and-seaborn-a48740bc9ce
    """
    labels = []
    for label in ax.get_xticklabels():
        text = label.get_text()
        labels.append(
            textwrap.fill(text, width=width, break_long_words=break_long_words)
        )
    ax.set_xticks(range(topN_GOterms_to_plot))  # new addition -- checking.
    ax.set_xticklabels(labels, rotation=40, ha="right", va="top")


def parse_user_proteomics_data(data: bytes) -> dict[str, tuple[float, float]]:
    isHeader = True  # The proteomics csv data file contains a header line. Set to 'False' otherwise.
    uniprotID_column = 0
    log2FC_column = 1
    adjPval_column = 2

    proteomicsData = {}
    csvReader = csv.reader(io.StringIO(data.decode("utf8")), delimiter=",")
    for row in csvReader:
        if isHeader:
            isHeader = False
            continue
        else:
            uniprotID = row[uniprotID_column]
            log2FC = row[log2FC_column]
            adjPVal = row[adjPval_column]
            proteomicsData[uniprotID] = (float(log2FC), float(adjPVal))

    return proteomicsData


def identify_perturbed_complexes(
    Complexes: ComplexT,
    ComplexNames: dict[str, str],
    proteomicsData: dict[str, tuple[float, float]],
    LOG2FC_THRESHOLD: float,
    ADJP_THRESHOLD: float,
) -> list[SubunitInfo]:
    allPerturbedComplexes = []
    for complexId in Complexes:
        complexMembers = Complexes[complexId]
        complexName = ComplexNames[complexId]

        # Check if at least one member of this complex has an omics measurement.
        OmicsMeasuredComplex = False
        measuredSubunits = []  # Has an omics measurement.
        regulatedSubunits = []  # Fulfills the criteria for a differentially expressed subunit.

        for member in complexMembers:
            if member in proteomicsData:
                measuredSubunits.append(member)
                if ONLY_REGULATED_SUBUNITS:
                    log2FCValue = abs(float(proteomicsData[member][0]))
                    adjPValue = round(float(proteomicsData[member][1]), 2)
                    if log2FCValue >= LOG2FC_THRESHOLD and adjPValue <= ADJP_THRESHOLD:
                        OmicsMeasuredComplex = True
                        regulatedSubunits.append(member)
                    else:
                        OmicsMeasuredComplex = True

        if OmicsMeasuredComplex:
            for subunit in regulatedSubunits:
                (log2fc, apval) = proteomicsData[subunit]
                subunit_info = SubunitInfo(
                    complex_id=complexId,
                    name=complexName,
                    subunit=subunit,
                    # measure=list(proteomicsData.values()),
                    log2fc=log2fc,
                    apvalue=apval,
                )
                allPerturbedComplexes.append(subunit_info)
    return allPerturbedComplexes
