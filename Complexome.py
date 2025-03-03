from typing import Optional
import io
import csv
import math
import sys
import json
import time
from dataclasses import dataclass
from functools import cache
from urllib.error import URLError
from urllib.request import urlopen, Request
from urllib.parse import urlencode
from pathlib import Path
from datetime import date
import textwrap

from matplotlib.axes import Axes
import matplotlib.pyplot as plt

MAX_RETRIES = 20
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
    taxon: str
    file: Path
    complexes: ComplexT
    complex_names: dict[str, str]
    complex_GO_terms: ComplexT
    proteomics_data: dict[str, tuple[float, float]]


def _urlretrieve_with_retries(
    request: Request, retries: int = 3, delay: float = 1.5
) -> tuple[dict[str, str], bytes]:
    for attempt in range(retries):
        try:
            with urlopen(request) as response:
                headers = dict(response.headers)
                return headers, response.read()

        except URLError:
            if attempt < retries - 1:
                time.sleep(delay)
            else:
                raise

    raise RuntimeError("Should never be raised")


def _fetch_genename_mapping(ids: set[str]) -> dict[str, str]:
    params = urlencode(
        {"from": "UniProtKB_AC-ID", "to": "Gene_Name", "ids": ",".join(ids)}
    ).encode("utf-8")
    request = Request(
        "https://rest.uniprot.org/idmapping/run", data=params, method="POST"
    )
    request.add_header("Content-Type", "application/x-www-form-urlencoded")
    _, data = _urlretrieve_with_retries(request)
    response = json.loads(data.decode("utf-8"))
    have_results = False
    retries = 0
    if (job := response.get("jobId")) is not None:
        request = Request(
            f"https://rest.uniprot.org/idmapping/status/{job}", method="GET"
        )
        while not have_results and retries < MAX_RETRIES:
            headers, data = _urlretrieve_with_retries(request)
            response = json.loads(data.decode("utf-8"))
            if "results" in response:
                get_results_url = headers.get("Link", "").split(";")[0].strip("<>")
                total_results = int(headers.get("X-Total-Results", "0"))
                have_results = True
            else:
                time.sleep(0.5)
            retries += 1

    expected = math.ceil(len(ids) / 25) + 10
    counter = 0
    results: dict[str, str] = {}
    while have_results and len(results) < total_results and counter < expected:
        headers, data = _urlretrieve_with_retries(
            Request(get_results_url, method="GET")
        )
        response = json.loads(data.decode("utf-8"))
        get_results_url = headers.get("Link", "").split(";")[0].strip("<>")
        if get_results_url == "":
            have_results = False
        results.update({el.get("from"): el.get("to") for el in response.get("results")})
        counter += 1

    return results


def setup(
    proteomics_data: dict[
        FileName, FileContents
    ],  # This is a mapping that we get from Colab
    organism_taxon_id: str = "9606",
) -> Complexome:
    if len(proteomics_data) > 1:
        raise RuntimeError("Cannot accept multiple input files")

    complexome_file = Path(
        organism_taxon_id + "Complexome_" + str(date.today()) + ".tsv"
    )  # Rename complexome file with a date stamp for future reference.

    if not complexome_file.exists():
        _, data = _urlretrieve_with_retries(
            Request(
                f"https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/{organism_taxon_id}.tsv"
            )
        )
        complexome_file.write_bytes(data)

    complexes, complex_names, complex_GO_terms = _parse_complexome_data(
        organism_taxon_id, str(complexome_file)
    )

    # protein_ids, _ = _unique_identities(list(complexes.values()))
    # gene_names = _fetch_genename_mapping(
    #     {protein_id.split("-")[0] for protein_id in protein_ids}
    # )

    parsed_proteomics_data = _parse_user_proteomics_data(
        list(proteomics_data.values())[0]
    )

    return Complexome(
        taxon=organism_taxon_id,
        file=complexome_file,
        complexes=complexes,
        complex_names=complex_names,
        complex_GO_terms=complex_GO_terms,
        proteomics_data=parsed_proteomics_data,
    )


def _add_complex_participants(participants: str, members=None) -> list[str]:
    if members is None:
        members = []

    for participant in participants.split("|"):
        participant_id = participant.split("(")[0]
        if "[" in participant_id:
            # This is in the case of molecule sets (paralogs that cannot be distinguished in this context).
            participant_ids = (
                str(participant_id).replace("[", "").replace("]", "").split(",")
            )
            for paralog in participant_ids:
                if paralog not in members:
                    members.append(paralog)
        else:
            if participant_id not in members:
                members.append(participant_id)

    return members


def _complex_subunit_numbers_distribution(complexes: ComplexT) -> dict[int, int]:
    num_subunits_per_complex = {}
    for key, value in complexes.items():
        if len(value) not in num_subunits_per_complex:
            num_subunits_per_complex[len(value)] = [key]
        else:
            num_subunits_per_complex[len(value)].append(key)

    return {key: len(value) for key, value in num_subunits_per_complex.items()}


def _parse_complexome_data(
    organism_taxon_id: str, complexome_file: str
) -> tuple[ComplexT, dict[str, str], ComplexT]:
    "Parse the Complex Portal data file and collect the per complex list of participants."
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
                raise RuntimeError(
                    f"Organism taxon id does not match with the expected organism! {row}"
                )

            # Extract information about the participants of the complex. Participants can include proteins (UniProtKB), chemical entities (ChEBI), RNA (RNAcentral) and complexes (Complex Portal).
            complexMembers = _add_complex_participants(complexParticipants)

            # Check if one or more participants are themselves complexes. In that case, the expanded list of protein members are contained in the Expanded participant list (last) column.
            if "CPX-" in complexParticipants:
                complexMembers = _add_complex_participants(
                    complexExtendedParticipants, complexMembers
                )

            # Add information about this complex and its participant molecules into the global dictionary of complexes.
            if complexID not in Complexes:
                Complexes[complexID] = complexMembers
                ComplexNames[complexID] = complexName
            else:
                raise RuntimeError(f"Complex ID already exists: {complexID}")

            # Extract information about the annotated GO terms of the complex.
            complexGOTerms = complexAnnotatedGOTerms.split("|")

            # Add information about this complex and its associated GO terms into the global dictionary of GO terms.
            if complexID not in ComplexesGOTerms:
                ComplexesGOTerms[complexID] = complexGOTerms
            else:
                raise RuntimeError(f"Complex ID already exists: {complexID}")

    return (Complexes, ComplexNames, ComplexesGOTerms)


def summary_statistics(complexome: Complexome) -> None:
    print("Total number of annotated complexes:", len(complexome.complexes))
    print("Total number of complex names:", len(complexome.complex_names))
    print(
        "Total number of complexes with associated GO terms:",
        len(complexome.complex_GO_terms),
    )


def _unique_identities(complexes: list[list[str]]) -> tuple[list[str], list[str]]:
    "Find the unique proteins, metabolites and RNA molecules contained in that complexome."
    unique_proteins = set()
    unique_metabolites = set()
    for cpx in complexes:
        for subunit in cpx:
            if "CPX-" in subunit:
                continue
            elif "URS" in subunit:
                continue
            elif "CHEBI:" in subunit:
                unique_metabolites.add(subunit)
            else:
                unique_proteins.add(subunit)

    return list(unique_proteins), list(unique_metabolites)


def plot(complexome: Complexome, axis: Optional[Axes] = None) -> Axes:
    if axis is None:
        axis = plt.subplot()

    x, height = zip(
        *_complex_subunit_numbers_distribution(complexome.complexes).items()
    )

    axis.bar(list(x), list(height), color="g")
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
    protein_subunits_per_complex = {}
    for complex_id, cplx in complexome.complexes.items():
        protein_subunits = {
            subunit
            for subunit in cplx
            if "CPX-" not in subunit
            and "URS" not in subunit
            and "CHEBI:" not in subunit
        }

        protein_subunits_per_complex[complex_id] = list(protein_subunits)

    x, height = zip(
        *_complex_subunit_numbers_distribution(protein_subunits_per_complex).items()
    )

    axis.bar(list(x), list(height), color="g")
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
    uniqueProteins, uniqueMetabolites = _unique_identities(
        list(complexome.complexes.values())
    )
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


def plot_volcano(
    complexome: Complexome,
    axis: Optional[Axes] = None,
    marker_size: Optional[int] = 20,
    log2fc_threshold: float = 2.0,
    adjp_threshold: float = 0.05,
) -> Axes:
    if axis is None:
        axis = plt.subplot()

    def colour(log2fc: float, pval: float) -> str:
        if pval > -math.log10(adjp_threshold):
            if log2fc > log2fc_threshold:
                return "tab:red"
            elif log2fc < -log2fc_threshold:
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


@cache
def protein_to_gene_name_mapping(proteins: tuple[str, ...]) -> dict[str, str]:
    return _fetch_genename_mapping(set(proteins))


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


def _parse_user_proteomics_data(data: bytes) -> dict[str, tuple[float, float]]:
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
