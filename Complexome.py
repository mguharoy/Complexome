"""Identify perturbed protein complexes."""

from __future__ import annotations

import csv
import enum
import io
import json
import math
import operator
import sys
import textwrap
import time
import typing
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import UTC, datetime
from functools import cache, reduce
from pathlib import Path
from urllib.error import URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import matplotlib.pyplot as plt
from matplotlib_venn import venn2  # type: ignore[import-untyped]

if typing.TYPE_CHECKING:
    from matplotlib.axes import Axes

MAX_RETRIES = 20
ONLY_REGULATED_SUBUNITS = True
UNIPROT_ID_COLUMN = 0
LOG2FC_COLUMN = 1
ADJPVAL_COLUMN = 2

ComplexT = dict[str, list[str]]
FileName = str
FileContents = bytes


class Perturbation(enum.StrEnum):
    """State of regulation."""

    UP_REGULATED = "Up-regulated"
    DOWN_REGULATED = "Down-regulated"
    ALTERED_COMPLEX = "Altered"
    UNINITIALIZED = "Unknown"

    def __add__(self, other: Perturbation) -> Perturbation:  # type: ignore[override]
        """Combine regulation states."""
        if self == self.UNINITIALIZED:
            return other
        if other in (self.UNINITIALIZED, self):
            return self
        return self.ALTERED_COMPLEX  # type: ignore[return-value]

    @classmethod
    def from_value(cls, value: float) -> Perturbation:
        """Create a regulation state from a value."""
        if value < 0:
            return cls.DOWN_REGULATED
        if value > 0:
            return cls.UP_REGULATED
        return cls.ALTERED_COMPLEX


@dataclass(frozen=True)
class PerturbationScore:
    """Complexome ranked scoring."""

    perturbation: Perturbation
    score: float
    score_normalized: float


@dataclass(frozen=True)
class SubunitInfo:
    """Information on a single sub-unit within a complex."""

    complex_id: str
    name: str
    subunit: str
    log2fc: float
    apvalue: float


@dataclass(frozen=True)
class OutputTableRow:
    """Data for a row in the complexes table."""

    complex_id: str
    complex_name: str
    coverage: float
    perturbationType: str
    perturbationScore: float
    perturbationScoreNormalized: float
    subunit: str
    genename: str
    log2fc: float
    apvalue: float


@dataclass(frozen=True)
class Complexome:
    """Stores information on the complexome for a taxon."""

    taxon: str
    file: Path
    complexes: ComplexT
    complex_names: dict[str, str]
    gene_ontology_terms: ComplexT
    proteomics_data: dict[str, tuple[float, float]]


class ComplexExistsError(Exception):
    """Error when adding a complex id that already exists."""

    def __init__(self, complex_id: str) -> None:
        """Init."""
        super().__init__(f"Complex ID already exists: {complex_id}")


class TaxonIdMismatchError(Exception):
    """Error when taxon ids mismatch while parsing."""

    def __init__(self, row: str, expected: str) -> None:
        super().__init__(
            "While parsing these data:\n"
            f"{row}\n"
            "I found an organism taxon ID does not "
            f"match with the expected organism ({expected})!"
        )


class ParsingConsistencyError(Exception):
    """There should be equal number of complexes, names, and GO terms."""

    def __init__(
        self, num_complexes: int, num_names: int, num_gene_ontology_terms: int
    ) -> None:
        super().__init__(f"""After parsing, I found an inconsistency:
        \t{num_complexes=}, {num_names=}, {num_gene_ontology_terms=}""")


class SetupError(Exception):
    """There was an error in the use of the setup() function."""

    def __init__(self) -> None:
        super().__init__("Cannot accept multiple input files.")


def _urlretrieve_with_retries(
    request: Request,
    retries: int = 3,
    delay: float = 1.5,
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

    raise RuntimeError


def _fetch_genename_mapping(ids: set[str]) -> dict[str, str]:
    results: dict[str, str] = {}
    params = urlencode(
        {"from": "UniProtKB_AC-ID", "to": "Gene_Name", "ids": ",".join(ids)},
    ).encode("utf-8")
    request = Request(
        "https://rest.uniprot.org/idmapping/run",
        data=params,
        method="POST",
    )
    request.add_header("Content-Type", "application/x-www-form-urlencoded")
    headers, data = _urlretrieve_with_retries(request)
    response = json.loads(data.decode("utf-8"))
    have_results = False
    retries = 0
    if (job := response.get("jobId")) is not None:
        request = Request(
            f"https://rest.uniprot.org/idmapping/status/{job}",
            method="GET",
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
            Request(get_results_url, method="GET"),
        )
        response = json.loads(data.decode("utf-8"))
        get_results_url = headers.get("Link", "").split(";")[0].strip("<>")
        results.update({el.get("from"): el.get("to") for el in response.get("results")})
        counter += 1

    return results


def setup(
    proteomics_data: dict[
        FileName,
        FileContents,
    ],  # This is a mapping that we get from Colab
    organism_taxon_id: str = "9606",
) -> Complexome:
    if len(proteomics_data) > 1:
        raise SetupError

    complexome_file = Path(
        organism_taxon_id + "Complexome_" + str(datetime.now(tz=UTC).today()) + ".tsv",
    )  # Rename complexome file with a date stamp for future reference.

    if not complexome_file.exists():
        _, data = _urlretrieve_with_retries(
            Request(
                f"https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/{organism_taxon_id}.tsv",
            ),
        )
        complexome_file.write_bytes(data)
    else:
        data = complexome_file.read_bytes()

    complexes, complex_names, gene_ontology_terms = _parse_complexome_data(
        organism_taxon_id,
        data.decode("utf-8"),
    )

    parsed_proteomics_data = _parse_user_proteomics_data(
        next(iter(proteomics_data.values())),
    )

    if len(complexes) != len(complex_names) or len(complexes) != len(
        gene_ontology_terms
    ):
        raise ParsingConsistencyError(
            num_complexes=len(complexes),
            num_names=len(complex_names),
            num_gene_ontology_terms=len(gene_ontology_terms),
        )

    return Complexome(
        taxon=organism_taxon_id,
        file=complexome_file,
        complexes=complexes,
        complex_names=complex_names,
        gene_ontology_terms=gene_ontology_terms,
        proteomics_data=parsed_proteomics_data,
    )


def _add_complex_participants(
    participants: str, members: list[str] | None = None
) -> list[str]:
    if members is None:
        members = []

    for participant in participants.split("|"):
        participant_id = participant.split("(")[0]
        if "[" in participant_id:
            # This is in the case of molecule sets
            # (paralogs that cannot be distinguished in this context).
            participant_ids = (
                str(participant_id).replace("[", "").replace("]", "").split(",")
            )
            for paralog in participant_ids:
                if paralog not in members:
                    members.append(paralog)
        elif participant_id not in members:
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
    organism_taxon_id: str,
    complexome_data: str,
) -> tuple[ComplexT, dict[str, str], ComplexT]:
    """Parse Complex Portal data and collect the per complex list of participants."""
    complexes = {}
    complex_names = {}
    complexes_gene_ontology_terms = {}

    for row in csv.reader(complexome_data.splitlines(), delimiter="\t"):
        if len(row) == 0 or row[0].startswith("#"):  # Skip the header line
            continue

        # Each row corresponds to a specific, annotated complex.
        complex_id = row[0]
        complex_name = row[1]
        complex_taxon_id = row[3]
        complex_participants = row[4]
        complex_extended_participants = row[-1]
        complex_annotated_gene_ontology_terms = row[7]

        # Make sure that the organism taxon id is matching.
        if complex_taxon_id != organism_taxon_id:
            raise TaxonIdMismatchError(str(row), organism_taxon_id)

        # Extract information about the participants of the complex.
        # Participants can include proteins (UniProtKB),
        # chemical entities (ChEBI), RNA (RNAcentral)
        # and complexes (Complex Portal).
        complex_members = _add_complex_participants(complex_participants)

        # Check if one or more participants are themselves complexes.
        # In that case, the expanded list of protein members are
        # contained in the Expanded participant list (last) column.
        if "CPX-" in complex_participants:
            complex_members = _add_complex_participants(
                complex_extended_participants,
                complex_members,
            )

        # Add information about this complex and its participant molecules
        # into the dictionary of complexes for this complexome.
        if complex_id not in complexes:
            complexes[complex_id] = complex_members
            complex_names[complex_id] = complex_name
        else:
            raise ComplexExistsError(complex_id)

        # Extract information about the annotated GO terms of the complex.
        gene_ontology_terms = complex_annotated_gene_ontology_terms.split("|")

        # Add information about this complex and its associated
        # gene ontology terms into the dictionary of gene
        # ontology terms for this complexome.
        if complex_id not in complexes_gene_ontology_terms:
            complexes_gene_ontology_terms[complex_id] = gene_ontology_terms
        else:
            raise ComplexExistsError(complex_id)

    return (complexes, complex_names, complexes_gene_ontology_terms)


def summary_statistics(complexome: Complexome) -> None:
    print(
        "Total number of annotated complexes in the downloaded dataset:",
        len(complexome.complexes),
    )


def _unique_identities(complexes: list[list[str]]) -> tuple[list[str], list[str]]:
    """Find the unique proteins, metabolites and RNA molecules in complexome."""
    unique_proteins = set()
    unique_metabolites = set()
    for cpx in complexes:
        for subunit in cpx:
            if "CPX-" in subunit or "URS" in subunit:
                continue
            if "CHEBI:" in subunit:
                unique_metabolites.add(subunit)
            else:
                unique_proteins.add(subunit)

    return list(unique_proteins), list(unique_metabolites)


def plot(complexome: Complexome, axis: Axes | None = None) -> Axes:
    if axis is None:
        axis = plt.subplot()

    x, height = zip(
        *_complex_subunit_numbers_distribution(complexome.complexes).items(),
        strict=False,
    )

    axis.bar(list(x), list(height), color="g")
    axis.set_xlabel("Number of subunits", fontsize=16)
    axis.set_ylabel("Number of complexes", fontsize=16)
    axis.set_xticks(axis.get_xticks(), size=14)
    axis.set_yticks(axis.get_yticks(), size=14)
    axis.set_title("Subunit distribution (proteins, metabolites, RNA)", fontsize=18)

    return axis


def proteins_only(complexome: Complexome, axis: Axes | None = None) -> Axes:
    if axis is None:
        axis = plt.subplot()
    # Store the protein subunits only per complex.
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
        *_complex_subunit_numbers_distribution(protein_subunits_per_complex).items(),
        strict=False,
    )

    axis.bar(list(x), list(height), color="g")
    axis.set_xlabel("Number of protein subunits", fontsize=16)
    axis.set_ylabel("Number of complexes", fontsize=16)
    axis.tick_params(axis="both", which="major", labelsize=14)
    axis.set_title("Subunit distribution (proteins only)", fontsize=18)
    return axis


def shared_protein_subunits(
    complexome: Complexome,
    axis: Axes | None = None,
) -> Axes:
    if axis is None:
        axis = plt.subplot()
    # Compute the distibution of shared protein subunits among the different complexes.
    unique_proteins, unique_metabolites = _unique_identities(
        list(complexome.complexes.values()),
    )
    protein_subunits_per_complex: dict[str, list[str]] = {}
    for complex_id, cmplx in complexome.complexes.items():
        all_protein_subunits = set()
        for subunit in cmplx:
            if "CPX-" in subunit or "URS" in subunit or "CHEBI:" in subunit:
                continue
            all_protein_subunits.add(subunit)
        protein_subunits_per_complex[complex_id] = list(all_protein_subunits)

    proteins_in_num_complexes: dict[str, int] = defaultdict(int)
    for protein in unique_proteins:
        for unique_protein_subunits in protein_subunits_per_complex.values():
            if protein in unique_protein_subunits:
                proteins_in_num_complexes[protein] += 1

    shared_subunits_per_complex: dict[int, int] = defaultdict(int)
    for value in proteins_in_num_complexes.values():
        # Useful to collate all values >10 into the last bar on the plot (see below).
        shared_subunits_per_complex[min(value, 10)] += 1

    axis.bar(
        list(shared_subunits_per_complex.keys()),
        list(shared_subunits_per_complex.values()),
        color="g",
    )
    axis.set_xlabel("Number of complexes", fontsize=16)
    axis.set_ylabel("Number of proteins", fontsize=16)
    axis.tick_params(axis="both", which="major", labelsize=14)
    axis.set_xticks(
        range(1, 11),
        ["1", "2", "3", "4", "5", "6", "7", "8", "9", ">10"],
        size=14,
    )
    axis.set_title("Distribution of shared protein subunits", fontsize=18)

    return axis


def proteomics_coverage_of_complexome(complexome: Complexome) -> dict[str, float]:
    proteomics_coverage_per_complex = {}
    for complex_id, cmplx in complexome.complexes.items():
        num_subunits = 0
        measured_subunits = 0
        for subunit in cmplx:
            if "CPX-" in subunit or "URS" in subunit or "CHEBI:" in subunit:
                continue
            num_subunits += 1
            if "-" in subunit or "_" in subunit:
                canonical_uniprot_id = subunit[:6]
            else:
                canonical_uniprot_id = subunit

            if canonical_uniprot_id in complexome.proteomics_data:
                measured_subunits += 1
        proteomics_coverage_per_complex[complex_id] = measured_subunits / num_subunits
    return proteomics_coverage_per_complex


def plot_proteomics_coverage_of_complexome(
    complexome: Complexome,
    axis: Axes | None = None,
) -> Axes:
    if axis is None:
        axis = plt.subplot()
    proteomics_coverage_of_complex = proteomics_coverage_of_complexome(complexome)

    _min = min(proteomics_coverage_of_complex.values())
    _max = max(proteomics_coverage_of_complex.values())
    if _min == _max:
        _max = _max + sys.float_info.epsilon
        bins = 1
        delta = 1.0
    else:
        bins = min(
            len(proteomics_coverage_of_complex),
            max(10, len(proteomics_coverage_of_complex) // 160),
        )
        delta = (_max - _min) / bins

    histogram: list[int] = []
    xs: list[float] = []
    _base = _min
    i = 1
    while _base < _max:
        top = (_max + sys.float_info.epsilon) if i == bins else (_base + delta)
        histogram.append(
            sum(
                1
                for _ in filter(
                    lambda x: _base <= x < top,
                    proteomics_coverage_of_complex.values(),
                )
            ),
        )
        xs.append(_base)
        if i == bins:
            _base = _max
        else:
            _base += delta
        i += 1

    assert sum(histogram) == len(proteomics_coverage_of_complex), (
        f"{sum(histogram)=} /= {len(proteomics_coverage_of_complex)}"
    )
    axis.bar(x=xs, height=histogram, width=delta, align="edge", color="g")
    axis.set_xbound(lower=0.0, upper=1.0)
    axis.set_xlabel("Proteomics Coverage", fontsize=16)
    axis.set_ylabel("Count", fontsize=16)
    axis.tick_params(axis="both", which="major", labelsize=14)
    return axis


def plot_venn_diagram(complexome: Complexome, axis: Axes | None = None) -> Axes:
    """Venn diagram of complexome proteins and proteins in the proteomics experiment."""
    if axis is None:
        axis = plt.subplot()
    all_canonical_subunits_complexome: set[str] = set()
    for cmplx in complexome.complexes.values():
        for subunit in cmplx:
            if "CPX-" in subunit or "URS" in subunit or "CHEBI:" in subunit:
                continue
            if "-" in subunit or "_" in subunit:
                all_canonical_subunits_complexome.add(subunit[:6])
            else:
                all_canonical_subunits_complexome.add(subunit)

    all_measured_protein_subunits = set(complexome.proteomics_data.keys())
    venn2(
        [all_canonical_subunits_complexome, all_measured_protein_subunits],
        set_labels=("Complexome proteins", "Proteomics dataset"),
        ax=axis,
    )

    return axis


def plot_volcano(
    complexome: Complexome,
    axis: Axes | None = None,
    marker_size: int | None = 20,
    log2fc_threshold: float = 2.0,
    adjp_threshold: float = 0.05,
) -> Axes:
    if axis is None:
        axis = plt.subplot()

    def colour(log2fc: float, pval: float) -> str:
        if pval > -math.log10(adjp_threshold):
            if log2fc > log2fc_threshold:
                return "tab:red"
            if log2fc < -log2fc_threshold:
                return "tab:blue"
            return "tab:gray"
        return "tab:gray"

    transformed_data = [
        (log2FC, -math.log10(adj_pval))
        for (log2FC, adj_pval) in complexome.proteomics_data.values()
    ]
    colours = [colour(*datum) for datum in transformed_data]
    (xvals, yvals) = zip(*transformed_data, strict=False)
    axis.scatter(xvals, yvals, s=marker_size, c=colours)
    axis.set_title("Volcano plot")
    axis.set_xlabel("log2 (FC)")
    axis.set_ylabel("-log10 (adjPval)")
    return axis


@cache
def protein_to_gene_name_mapping(proteins: tuple[str, ...]) -> dict[str, str]:
    return _fetch_genename_mapping(set(proteins))


def wrap_labels(
    ax: Axes, width: int, top_n: int, break_long_words: bool = False, fontsize: int = 10
) -> None:
    """Word wrap long GO term names (used as axis labels).

    Based on https://medium.com/dunder-data/automatically-wrap-graph-labels-in-matplotlib-and-seaborn-a48740bc9ce
    """
    labels = []
    for label in ax.get_xticklabels():
        text = label.get_text()
        labels.append(
            textwrap.fill(text, width=width, break_long_words=break_long_words),
        )
    if len(labels) < top_n:
        ax.set_xticks(range(len(labels)))
    else:
        ax.set_xticks(range(top_n))  # new addition -- checking. was needed for cases where the #GO terms to be plotted were less than the specified top_n.
    ax.set_xticklabels(labels, rotation=40, ha="right", va="top", fontsize=fontsize)


def _parse_user_proteomics_data(data: bytes) -> dict[str, tuple[float, float]]:
    is_header = True  # The proteomics csv data file contains a header line.
    toIgnore = ['NA', 'NaN', float('NaN')]

    proteomics_data = {}
    reader = csv.reader(io.StringIO(data.decode("utf8")), delimiter=",")
    for row in reader:
        if is_header:
            is_header = False
            continue
        uniprot_id = row[UNIPROT_ID_COLUMN]
        log2fc = row[LOG2FC_COLUMN]
        adj_pval = row[ADJPVAL_COLUMN]
        if log2fc in toIgnore or adj_pval in toIgnore:
            continue
        proteomics_data[uniprot_id] = (float(log2fc), float(adj_pval))

    return proteomics_data


def identify_perturbed_complexes(
    complexes: ComplexT,
    names: dict[str, str],
    proteomics_data: dict[str, tuple[float, float]],
    log2fc_threshold: float,
    adjp_threshold: float,
) -> list[SubunitInfo]:
    all_perturbed_complexes = []
    for complex_id, members in complexes.items():
        complex_name = names[complex_id]

        # Check if at least one member of this complex has an omics measurement.
        omics_measured_complex = False
        measured_subunits = []  # Has an omics measurement.

        # Fulfills the criteria for a differentially expressed subunit.
        regulated_subunits = []

        for member in members:
            if member in proteomics_data:
                measured_subunits.append(member)
                if ONLY_REGULATED_SUBUNITS:
                    log2_fc_value = abs(float(proteomics_data[member][0]))
                    adj_pvalue = round(float(proteomics_data[member][1]), 2)
                    if (
                        log2_fc_value >= log2fc_threshold
                        and adj_pvalue <= adjp_threshold
                    ):
                        omics_measured_complex = True
                        regulated_subunits.append(member)
                    else:
                        omics_measured_complex = True

        if omics_measured_complex:
            for subunit in regulated_subunits:
                (log2fc, apval) = proteomics_data[subunit]
                subunit_info = SubunitInfo(
                    complex_id=complex_id,
                    name=complex_name,
                    subunit=subunit,
                    log2fc=log2fc,
                    apvalue=apval,
                )
                all_perturbed_complexes.append(subunit_info)
    return all_perturbed_complexes


def plot_gene_ontology_analysis_perturbed_complexes(
    complexome: Complexome,
    all_perturbed_complexes: list[SubunitInfo],
    top_n: int,
    axis: Axes | None = None,
) -> Axes:
    """Plot the top_n most perturbed complexes."""
    if axis is None:
        axis = plt.subplot()
    perturbed_complex_ids = {subunit.complex_id for subunit in all_perturbed_complexes}

    # Iterate over the perturbed complexes and identify their associated GO terms.
    affected_complexes_all_gene_ontology_terms = Counter(
        go_term
        for complex_id in perturbed_complex_ids
        for go_term in complexome.gene_ontology_terms[complex_id]
    )
    top = dict(affected_complexes_all_gene_ontology_terms.most_common(top_n))

    axis.bar(
        list(top.keys()),
        list(top.values()),
        color="g",
    )
    wrap_labels(axis, 30, top_n, fontsize=8)
    axis.set_ylabel("# occurrences", fontsize=12)
    axis.set_title(
        "Most frequent GO terms associated with the perturbed complexes",
        fontsize=14,
    )

    return axis


def regulation(log2fc: list[float]) -> Perturbation:
    """Compute up or down regulation.

    If all log2fc values are negative then down regulated.
    If all log2fc values are positive then up regulated.
    Otherwise altered complex.
    """
    return reduce(
        lambda acc, val: acc + Perturbation.from_value(val),
        log2fc,
        Perturbation.UNINITIALIZED,
    )


def perturbation_scores(
    complexome: Complexome,
    all_perturbed_complexes: list[SubunitInfo],
) -> dict[str, PerturbationScore]:
    """Scoring function."""
    perturbation_scores = {}
    perturbed_complex_ids = {subunit.complex_id for subunit in all_perturbed_complexes}
    for complex_id, cmplx in complexome.complexes.items():
        if complex_id not in perturbed_complex_ids:
            continue

        log2_fc_values = []
        adjp_values = []
        for subunit in cmplx:
            if "CPX-" in subunit or "URS" in subunit or "CHEBI:" in subunit:
                continue
            if "-" in subunit or "_" in subunit:
                subunit_protein_id = subunit[
                    :6
                ]  # Maybe a bit of an assumption, but ok for now.
            else:
                subunit_protein_id = subunit

            if subunit_protein_id in complexome.proteomics_data:
                log2_fc_values.append(complexome.proteomics_data[subunit_protein_id][0])
                adjp_values.append(complexome.proteomics_data[subunit_protein_id][1])

        perturbation = regulation(log2_fc_values)

        regulation_score = sum(
            [
                abs(fc * -math.log10(pval))
                for fc, pval in zip(log2_fc_values, adjp_values, strict=False)
            ],
        )
        regulation_score_normalized = regulation_score / float(len(log2_fc_values))

        perturbation_scores[complex_id] = PerturbationScore(
            perturbation=perturbation,
            score=regulation_score,
            score_normalized=regulation_score_normalized,
        )
    return perturbation_scores


def format_output_table_data(
    complexome: Complexome, perturbed_complex_subunits: list[SubunitInfo], key: str
) -> list[list[str]]:
    genes = protein_to_gene_name_mapping(
        tuple({info.subunit for info in perturbed_complex_subunits})
    )
    coverage = proteomics_coverage_of_complexome(complexome)
    perturbation = perturbation_scores(complexome, perturbed_complex_subunits)
    data = sorted(
        [
            OutputTableRow(
                complex_id=info.complex_id,
                complex_name=info.name,
                coverage=coverage.get(info.complex_id) or 0.0,
                perturbationType=perturbation.get(info.complex_id).perturbation,
                perturbationScore=perturbation.get(info.complex_id).score,
                perturbationScoreNormalized=perturbation.get(info.complex_id).score_normalized,
                subunit=info.subunit,
                genename=genes.get(info.subunit) or "",
                log2fc=info.log2fc,
                apvalue=info.apvalue,
            )
            for info in perturbed_complex_subunits
        ],
        key=operator.attrgetter(key),
        reverse=True,
    )
    return [
        [
            "Complex ID",
            "Complex Name",
            "Coverage",
            "Perturbation Type",
            "Perturbation Score",
            "Normalized Score",
            "Subunit Id",
            "Gene Name",
            "Log2FC",
            "Adj P-value",
        ]
    ] + [
        [
            info.complex_id,
            info.complex_name,
            f"{info.coverage:.2f}",
            info.perturbationType,
            f"{info.perturbationScore:.2f}",
            f"{info.perturbationScoreNormalized:.2f}",
            info.subunit,
            info.genename,
            f"{info.log2fc:.2f}",
            f"{info.apvalue:.3f}",
        ]
        for info in data
    ]
