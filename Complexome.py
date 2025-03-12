from typing import Optional
import io
import csv
import math
import sys
import json
import time
import operator
from dataclasses import dataclass
from functools import cache
from collections import Counter
from urllib.error import URLError
from urllib.request import urlopen, Request
from urllib.parse import urlencode
from pathlib import Path
from datetime import date
import textwrap

from matplotlib.axes import Axes
import matplotlib.pyplot as plt
from matplotlib_venn import venn2  # type: ignore[import-untyped]

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
class OutputTableRow:
    complex_id: str
    complex_name: str
    coverage: float
    subunit: str
    genename: str
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
    results: dict[str, str] = {}
    params = urlencode(
        {"from": "UniProtKB_AC-ID", "to": "Gene_Name", "ids": ",".join(ids)}
    ).encode("utf-8")
    request = Request(
        "https://rest.uniprot.org/idmapping/run", data=params, method="POST"
    )
    request.add_header("Content-Type", "application/x-www-form-urlencoded")
    headers, data = _urlretrieve_with_retries(request)
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
                results = {
                    el.get("from"): el.get("to") for el in response.get("results")
                }
                total_results = int(headers.get("X-Total-Results", "0"))
                have_results = True
            else:
                time.sleep(0.5)
            retries += 1

    expected = math.ceil(len(ids) / 25) + 10
    counter = 0

    while (
        have_results
        and get_results_url != ""
        and len(results) < total_results
        and counter < expected
    ):
        headers, data = _urlretrieve_with_retries(
            Request(get_results_url, method="GET")
        )
        response = json.loads(data.decode("utf-8"))
        get_results_url = headers.get("Link", "").split(";")[0].strip("<>")
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
    else:
        data = complexome_file.read_bytes()

    complexes, complex_names, complex_GO_terms = _parse_complexome_data(
        organism_taxon_id, data.decode("utf-8")
    )

    parsed_proteomics_data = _parse_user_proteomics_data(
        list(proteomics_data.values())[0]
    )

    if len(complexes) != len(complex_names) or len(complexes) != len(complex_GO_terms):
        raise ValueError(
            "Error: Inconsistency in parsing complex IDs, complex names and complex GO terms."
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
    organism_taxon_id: str, complexome_data: str
) -> tuple[ComplexT, dict[str, str], ComplexT]:
    "Parse the Complex Portal data file and collect the per complex list of participants."
    complexes = {}
    complex_names = {}
    complexes_GO_terms = {}

    for row in csv.reader(complexome_data.splitlines(), delimiter="\t"):
        if len(row) == 0 or row[0].startswith("#"):  # Skip the header line
            continue

        # Each row corresponds to a specific, annotated complex.
        complex_id = row[0]
        complex_name = row[1]
        complex_taxon_id = row[3]
        complex_participants = row[4]
        complex_extended_participants = row[-1]
        complex_annotated_GO_terms = row[7]

        # Make sure that the organism taxon id is matching.
        if complex_taxon_id != organism_taxon_id:
            raise RuntimeError(
                f"Organism taxon id does not match with the expected organism! {row}"
            )

        # Extract information about the participants of the complex.
        # Participants can include proteins (UniProtKB), chemical entities (ChEBI), RNA (RNAcentral) and complexes (Complex Portal).
        complex_members = _add_complex_participants(complex_participants)

        # Check if one or more participants are themselves complexes.
        # In that case, the expanded list of protein members are contained in the Expanded participant list (last) column.
        if "CPX-" in complex_participants:
            complex_members = _add_complex_participants(
                complex_extended_participants, complex_members
            )

        # Add information about this complex and its participant molecules into the global dictionary of complexes.
        if complex_id not in complexes:
            complexes[complex_id] = complex_members
            complex_names[complex_id] = complex_name
        else:
            raise RuntimeError(f"Complex ID already exists: {complex_id}")

        # Extract information about the annotated GO terms of the complex.
        complex_GO_terms = complex_annotated_GO_terms.split("|")

        # Add information about this complex and its associated GO terms into the global dictionary of GO terms.
        if complex_id not in complexes_GO_terms:
            complexes_GO_terms[complex_id] = complex_GO_terms
        else:
            raise RuntimeError(f"Complex ID already exists: {complex_id}")

    return (complexes, complex_names, complexes_GO_terms)


def summary_statistics(complexome: Complexome) -> None:
    print(
        "Total number of annotated complexes in the downloaded dataset:",
        len(complexome.complexes),
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
    axis.tick_params(axis="both", which="major", labelsize=14)
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
    axis.tick_params(axis="both", which="major", labelsize=14)
    axis.set_xticks(
        range(1, 11), ["1", "2", "3", "4", "5", "6", "7", "8", "9", ">10"], size=14
    )
    axis.set_title("Distribution of shared protein subunits", fontsize=18)

    return axis


def proteomics_coverage_of_complexome(complexome: Complexome) -> dict[str, float]:
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
                if "-" in subunit or "_" in subunit:
                    canonicalUniProtID = subunit[:6]
                else:
                    canonicalUniProtID = subunit

                if canonicalUniProtID in complexome.proteomics_data:
                    measuredSubunits += 1
        proteomicsCoveragePerComplex[complex_id] = measuredSubunits / numSubunits
    return proteomicsCoveragePerComplex

def plot_proteomics_coverage_of_complexome(
    complexome: Complexome, axis: Optional[Axes] = None
) -> Axes:
    if axis is None:
        axis = plt.subplot()
    proteomicsCoveragePerComplex: dict[str, float] = proteomics_coverage_of_complexome(complexome)

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
    axis.set_xlabel("Proteomics Coverage", fontsize=16)
    axis.set_ylabel("Count", fontsize=16)
    axis.tick_params(axis="both", which="major", labelsize=14)
    return axis


def plot_venn_diagram(complexome: Complexome, axis: Optional[Axes] = None) -> Axes:
    "Visualize the overlap between the complexome proteins and the dataset of proteins measured in the proteomics experiment."
    if axis is None:
        axis = plt.subplot()
    all_canonical_subunits_complexome: set[str] = set()
    for complex_id, complex in complexome.complexes.items():
        for subunit in complex:
            if "CPX-" in subunit:
                continue
            elif "URS" in subunit:
                continue
            elif "CHEBI:" in subunit:
                continue
            else:
                if "-" in subunit or "_" in subunit:
                    all_canonical_subunits_complexome.add(subunit[:6])
                else:
                    all_canonical_subunits_complexome.add(subunit)

    all_measured_protein_subunits = set(complexome.proteomics_data.keys())
    venn2(
        [all_canonical_subunits_complexome, all_measured_protein_subunits],
        set_labels=("Complexome proteins", "Proteomics dataset"),
        ax=axis
    )

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


def wrap_labels(ax, width, topN_GOterms_to_plot, break_long_words=False, fontsize=10):
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
    ax.set_xticklabels(labels, rotation=40, ha="right", va="top", fontsize=fontsize)


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
                    log2fc=log2fc,
                    apvalue=apval,
                )
                allPerturbedComplexes.append(subunit_info)
    return allPerturbedComplexes


def GO_analysis_perturbed_complexes(
    complexome: Complexome,
    all_perturbed_complexes: list[SubunitInfo],
    top_n: int,
    axis: Optional[Axes] = None,
):
    if axis is None:
        axis = plt.subplot()
    perturbed_complex_ids = {subunit.complex_id for subunit in all_perturbed_complexes}

    # Iterate over the perturbed complexes and identify their associated GO terms.
    affected_complexes_all_GO_terms = Counter(
        (
            go_term
            for complex_id in perturbed_complex_ids
            for go_term in complexome.complex_GO_terms[complex_id]
        )
    )
    top = dict(affected_complexes_all_GO_terms.most_common(top_n))

    axis.bar(
        list(top.keys()),
        list(top.values()),
        color="g",
    )
    wrap_labels(axis, 30, top_n, fontsize=6)
    axis.set_ylabel("# occurrences", fontsize=10)
    axis.set_title(
        f"Top {top_n} most frequent GO terms associated with the perturbed complexes",
        fontsize=10,
    )

    return axis


def perturbation_score_calculator(
    complexome: Complexome,
    all_perturbed_complexes: list[SubunitInfo],
) -> dict[str, list[str, float, float]]:

    perturbation_scores = {}
    perturbed_complex_ids = {subunit.complex_id for subunit in all_perturbed_complexes}
    for complex_id, complex in complexome.complexes.items():
        if complex_id not in perturbed_complex_ids:
            continue

        log2_fc_values = []
        adjp_values = []
        for subunit in complex:
            if "CPX-" in subunit:
                continue
            elif "URS" in subunit:
                continue
            elif "CHEBI:" in subunit:
                continue
            else:
                if "-" in subunit or "_" in subunit:
                    subunit_protein_id = subunit[:6] # Maybe a bit of an assumption, but ok for now.
                else:
                    subunit_protein_id = subunit

            if subunit_protein_id in complexome.proteomics_data:
                log2_fc_values.append(complexome.proteomics_data[subunit_protein_id][0])
                adjp_values.append(complexome.proteomics_data[subunit_protein_id][1])

        allDownRegulated = all(n < 0 for n in log2_fc_values)
        allUpRegulated = all(n > 0 for n in log2_fc_values)
        alteredComplex = False
        if not allDownRegulated and not allUpRegulated:
            alteredComplex = True

        if allDownRegulated:
            perturbationType='Down-regulated'
        elif allUpRegulated:
            perturbationType='Up-regulated'
        elif alteredComplex:
            perturbationType='Altered'

        regulation_score=sum([abs(a * -math.log10(b)) for a, b in zip(log2_fc_values, adjp_values)])
        regulation_score_normalized = regulation_score/float(len(log2_fc_values))

        perturbation_scores[complex_id] = [perturbationType, float(regulation_score), float(regulation_score_normalized)]
    return perturbation_scores

def format_output_table_data(complexome: Complexome, perturbed_complex_subunits: list[SubunitInfo], key: str) -> list[list[str]]:
    genes = protein_to_gene_name_mapping(tuple({info.subunit for info in perturbed_complex_subunits}))
    coverage = proteomics_coverage_of_complexome(complexome)
    data = sorted([OutputTableRow(complex_id = info.complex_id, complex_name = info.name, coverage = coverage.get(info.complex_id), subunit = info.subunit, genename = genes.get(info.subunit), log2fc = info.log2fc, apvalue = info.apvalue) for info in perturbed_complex_subunits], key=operator.attrgetter(key), reverse=True)
    return [["Complex ID", "Complex Name", "Coverage", "Subunit Id", "Gene Name", "Log2FC","Adj P-value"]] + [[info.complex_id, info.complex_name, f"{info.coverage:.2f}", info.subunit, info.genename, f"{info.log2fc:.2f}", f"{info.apvalue:.3f}"] for info in data]

