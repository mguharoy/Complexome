// @ts-check
/**
 * Set up the workflow structure with event listeners.
 * Also run the worker thread to decode CSV documents.
 * And the service worker to cache everything.
 */

import {
  barPlot,
  histogramPlot,
  table,
  vennPlot,
  volcanoPlot,
} from "./plot.mjs";

/**
 * @import { Sorting, TableRow } from "./plot.mjs";
 */

/**
 * Handle changes to the user supplied proteoimcs file.
 * @param {Worker} csv - A reference to the worker thread for parsing user data.
 */
async function onProteomicsFile(csv) {
  const userfile = /** @type {HTMLInputElement | null} */ (
    /** @type {unknown} */ document.querySelector("#proteomics-file")
  );
  const files = userfile?.files;
  if (files && files.length === 1) {
    const text = await files[0]?.text();
    if (text) {
      csv.postMessage({ op: "csv", data: text });
    }
  }
}

async function registerServiceWorker() {
  if ("serviceWorker" in navigator) {
    try {
      const registration = await navigator.serviceWorker.register("cache.mjs", {
        scope: "./",
      });
      if (registration.installing) {
        console.log("Service worker installing...");
      } else if (registration.waiting) {
        console.log("Service worker installed.");
      } else if (registration.active) {
        console.log("Service worker active.");
      }
    } catch (error) {
      console.error(`Service Worker registration failed with: ${error}`);
    }
  }
}

/**
 * @typedef Mapping
 * @property from {string}
 * @property to {string}
 */

/**
 * @param {Subunit[]} perturbedSubunits
 * @returns {Promise<Map<string, string>>}
 */
async function fetchGeneNameMapping(perturbedSubunits) {
  const results = new Map();
  let haveResults = false;
  let retries = 0;
  const ids = new Set(perturbedSubunits.map((info) => info.subunit));
  const request = new Request("https://rest.uniprot.org/idmapping/run", {
    method: "POST",
    headers: {
      "Content-Type": "application/x-www-form-urlencoded; charset=UTF-8",
    },
    body: new URLSearchParams({
      from: "UniProtKB_AC-ID",
      to: "Gene_Name",
      ids: ids.values().toArray().join(","),
    }),
  });
  const response = await fetch(request);
  const job = await response.json();
  /** @type {string | undefined} */
  let resultsURL = undefined;

  let totalResults = 0;
  if ("jobId" in job) {
    while (!haveResults && retries < 20) {
      const status = await fetch(
        `https://rest.uniprot.org/idmapping/status/${job.jobId}`,
      );
      const headers = status.headers;
      const statusResponse = await status.json();
      if ("results" in statusResponse) {
        resultsURL = (headers.get("Link") ?? "")
          .split(";")[0]
          ?.replace(/^<|>$/g, "");
        statusResponse.results.forEach((/** @type {Mapping} */ mapping) =>
          results.set(mapping.from, mapping.to),
        );
        totalResults = parseInt(headers.get("X-Total-Results") ?? "0");
        haveResults = true;
      } else {
        await new Promise((resolve) => setTimeout(resolve, 500));
      }
      retries++;
    }
  }

  const maxExpected = Math.ceil(ids.size / 25) + 10;
  let counter = 0;

  while (
    haveResults &&
    resultsURL &&
    results.size < totalResults &&
    counter < maxExpected
  ) {
    const resultsResponse = await fetch(resultsURL);
    const mappings = await resultsResponse.json();
    resultsURL = (resultsResponse.headers.get("Link") ?? "")
      .split(";")[0]
      ?.replace(/^<|>$/g, "");
    mappings.results.forEach((/** @type {Mapping} */ mapping) =>
      results.set(mapping.from, mapping.to),
    );
    counter++;
  }

  return results;
}

/**
 * Enable or disable controls.
 * @param {boolean} enable
 */
function controls(enable) {
  const log2fc = document.getElementById("log2fc-threshold");
  const adjp = document.getElementById("adjp-threshold");
  const goterms = document.getElementById("top-n-go-terms");

  if (log2fc && "disabled" in log2fc) {
    log2fc.disabled = !enable;
  }
  if (adjp && "disabled" in adjp) {
    adjp.disabled = !enable;
  }
  if (goterms && "disabled" in goterms) {
    goterms.disabled = !enable;
  }
}

/**
 * Get a value form a Map with a default if it doesn't exist.
 * @template K
 * @template V
 * @param {Map<K, V>} map
 * @param {K} key
 * @param {V} dflt
 * @returns {V}
 */
function mapGetWithDefault(map, key, dflt) {
  return map.get(key) ?? dflt;
}

/**
 * Set a value in a Map with a default if it doesn't exist.
 * @template K
 * @template V
 * @param {Map<K, V>} map
 * @param {K} key
 * @param {(value: V) => V} updater
 * @param {V} dflt
 * @returns {Map<K, V>}
 */
function mapSetWithDefault(map, key, updater, dflt) {
  const current = map.get(key);
  if (current === undefined) {
    map.set(key, dflt);
    return map;
  } else {
    const value = updater(map.get(key) ?? dflt);
    map.set(key, value);
    return map;
  }
}

function subunitDistributionPlot() {
  const count = new Map();
  const complexes = window.complexome ? window.complexome[0] : [];

  for (const cplx of complexes.values()) {
    let subunits = Array.from(cplx.values()).filter(
      (value) =>
        !value.includes("CPX-") &&
        !value.includes("URS") &&
        !value.includes("CHEBI:"),
    );
    count.set(
      subunits.length,
      mapGetWithDefault(count, subunits.length, 0) + 1,
    );
  }
  return barPlot(
    Array.from(count.entries())
      .sort((a, b) => +a[0] - +b[0])
      .map(([x, y]) => [x.toString(), y]),
    {
      margin: [30, 0, 30, 30],
      width: 600,
      height: 400,
      xlabel: "Number of complexes",
      ylabel: "↑ Frequency",
      title: "Subunit distribution (proteins only)",
    },
  );
}

function sharedSubunitsPlot() {
  /** @type {Map<string, number>} */
  const count = new Map();
  const proteinObsCount = new Map();
  const complexes = window.complexome ? window.complexome[0] : [];

  for (const cplx of complexes.values()) {
    let subunits = Array.from(cplx.values()).filter(
      (value) =>
        !value.includes("CPX-") &&
        !value.includes("URS") &&
        !value.includes("CHEBI:"),
    );
    for (const subunit of subunits) {
      proteinObsCount.set(
        subunit,
        mapGetWithDefault(proteinObsCount, subunit, 0) + 1,
      );
    }
  }

  for (const value of proteinObsCount.values()) {
    const key = Math.min(value, 10) <= 9 ? `${value}` : ">10";
    count.set(key, mapGetWithDefault(count, key, 0) + 1);
  }

  return barPlot(Array.from(count.entries()).sort(), {
    margin: [30, 0, 30, 40],
    width: 600,
    height: 400,
    xlabel: "Number of protein subunits",
    ylabel: "↑ Frequency",
    title: "Shared protein subunits",
  });
}

function coveragePlot() {
  const coverage = new Map();
  const complexes = window.complexome ? window.complexome[0] : [];

  for (const [complexID, cplx] of complexes) {
    const subunits = Array.from(cplx.values())
      .filter(
        (subunit) =>
          !subunit.includes("CPX-") &&
          !subunit.includes("URS") &&
          !subunit.includes("CHEBI:"),
      )
      .map((subunit) =>
        subunit.includes("-") || subunit.includes("_")
          ? subunit.slice(6)
          : subunit,
      );
    const numSubunits = subunits.length;
    const measuredSubunits = subunits.filter((subunit) =>
      window?.userdata?.has(subunit),
    ).length;
    coverage.set(complexID, measuredSubunits / numSubunits);
  }

  return histogramPlot(Array.from(coverage.values()), {
    margin: [30, 10, 30, 50],
    width: 600,
    height: 400,
    xlabel: "Proteomics coverage",
    ylabel: "↑ Frequency",
    title: "",
  });
}

function vennDiagram() {
  const allCanonicalSubunits = new Set();
  const complexes = window.complexome ? window.complexome[0] : [];

  for (const cplx of complexes.values()) {
    cplx
      .values()
      .filter(
        (subunit) =>
          !subunit.includes("CPX-") &&
          !subunit.includes("URS") &&
          !subunit.includes("CHEBI:"),
      )
      .map((subunit) =>
        subunit.includes("-") || subunit.includes("_")
          ? subunit.slice(6)
          : subunit,
      )
      .forEach((subunit) => {
        allCanonicalSubunits.add(subunit);
      });
  }
  return vennPlot(
    { A: allCanonicalSubunits, B: new Set(window?.userdata?.keys()) },
    {
      margin: [0, 0, 20, 20],
      width: 600,
      height: 400,
      alabel: "Complexome proteins",
      blabel: "Proteomics dataset",
      title: "",
    },
  );
}

function volcano() {
  const data = Array.from(
    window?.userdata?.entries().map(([name, [log2fc, adjpval]]) => {
      return { adjPval: adjpval, log2fc: log2fc, name: name };
    }) ?? [],
  );
  const log2fc = /** @type {HTMLInputElement | null} */ (
    /** @type {unknown} */ document.getElementById("log2fc-threshold")
  );
  const adjp = /** @type {HTMLInputElement | null} */ (
    /** @type {unknown} */ document.getElementById("adjp-threshold")
  );

  return volcanoPlot(data, {
    margin: [30, 10, 30, 30],
    width: 600,
    height: 400,
    xlabel: "log2 (FC)",
    ylabel: "-log10 (adjPval)",
    title: "",
    log2fcThreshold: parseFloat(log2fc?.value ?? "0"),
    adjpThreshold: parseFloat(adjp?.value ?? "0"),
  });
}

/**
 * @typedef Subunit
 * @property name {string}
 * @property subunit {string}
 * @property log2fc {number}
 * @property apvalue {number}
 */

function goTerms() {
  const log2fc = /** @type {HTMLInputElement | null} */ (
    /** @type {unknown} */ document.getElementById("log2fc-threshold")
  );
  const adjp = /** @type {HTMLInputElement | null} */ (
    /** @type {unknown} */ document.getElementById("adjp-threshold")
  );
  const numTermsEl = /** @type {HTMLInputElement | null} */ (
    /** @type {unknown} */ document.getElementById("top-n-go-terms")
  );

  const numTerms = parseInt(numTermsEl?.value ?? "10");

  const log2fcThreshold = parseFloat(log2fc?.value ?? "0");
  const adjpThreshold = parseFloat(adjp?.value ?? "0");

  /** @type {Map<string, Set<string>>} */
  const complexes = window.complexome ? window.complexome[0] : new Map();

  /** @type {Map<string, string[]>} */
  const goTerms = window.complexome ? window.complexome[2] : new Map();

  /** @type {Map<string, [number, number]>} */
  const proteomics = window.userdata ? window.userdata : new Map();

  let observedComplexes = new Set();

  /** @type {Array<[string, number]>} */
  const counts = complexes
    .entries()
    .flatMap(
      /** @type {([cpxid, members]: [string, Set<string>]) => Iterator<[string, string]>} */ ([
        cpxid,
        members,
      ]) => members.values().map((member) => [cpxid, member]),
    )
    .map(([cpxid, member]) => {
      const [measured_log2fc, measured_adjp] = proteomics.get(member) ?? [
        0,
        Infinity,
      ];
      return {
        name: cpxid,
        subunit: member,
        log2fc: measured_log2fc,
        apvalue: measured_adjp,
      };
    })
    .filter(
      (/** @type {Subunit} */ subunit) =>
        Math.abs(subunit.log2fc) >= log2fcThreshold &&
        subunit.apvalue <= adjpThreshold,
    )
    .flatMap((subunit) => {
      if (observedComplexes.has(subunit.name)) {
        return [];
      } else {
        observedComplexes.add(subunit.name);
        return goTerms.get(subunit.name) ?? [];
      }
    })
    .reduce((/** @type {Map<string, number>} */ acc, goTerm) => {
      return mapSetWithDefault(acc, goTerm, (count) => count + 1, 1);
    }, new Map())
    .entries()
    .toArray()
    .sort((a, b) => +b[1] - +a[1])
    .slice(0, numTerms);

  console.log(counts);

  const bottomMargin = counts.reduce(
    (acc, label) => Math.max(acc, label[0].length * 4.5),
    0,
  );
  const leftMargin =
    (counts[0] ?? ["", 0])[0].length * 4.5 * Math.cos(Math.PI / 3) - 30;
  console.log(bottomMargin);

  return barPlot(counts, {
    margin: [25, 0, bottomMargin, leftMargin],
    width: leftMargin + 600,
    height: bottomMargin + 400,
    xlabel: "",
    ylabel: "↑ Frequency",
    title: "Most frequent GO terms",
  });
}

/**
 * @typedef Perturbation
 * @property perturbation {string}
 * @property score {number}
 * @property normalizedScore {number}
 */

/**
 * @param {string} other
 * @param {number} log2fc
 * @returns {string}
 */
function whatPerturbation(other, log2fc) {
  let self = "";
  switch (Math.sign(log2fc)) {
    case -1:
      self = "Down-regulated";
      break;
    case 1:
      self = "Up-regulated";
      break;
    case 0:
      self = "Altered";
      break;
  }
  if (self === "Unknown") {
    return other;
  } else if (other === "Unknown" || other === self) {
    return self;
  } else {
    return "Altered";
  }
}

/**
 * @param {keyof TableRow} column
 * @param {"asc" | "desc"} order
 * @returns {(a: TableRow, b: TableRow) => number}
 */
function compare(column, order) {
  console.log(column, order);
  if (order == "asc") {
    return (a, b) => {
      if (typeof a[column] === "number") {
        return +a[column] - +b[column];
      } else {
        return a[column] < b[column] ? -1 : 1;
      }
    };
  } else {
    return (a, b) => {
      if (typeof a[column] === "number") {
        return +b[column] - +a[column];
      } else {
        return b[column] < a[column] ? -1 : 1;
      }
    };
  }
}

/**
 * Display the data table.
 * @param {(cid: string) => void} viewComplex - Render a complex viewer.
 * @param {Sorting} sorting - Which column to sort in which order.
 * @returns {Promise<Node[]>}
 */
async function dataTable(viewComplex, sorting) {
  const log2fc = /** @type {HTMLInputElement | null} */ (
    /** @type {unknown} */ document.getElementById("log2fc-threshold")
  );
  const adjp = /** @type {HTMLInputElement | null} */ (
    /** @type {unknown} */ document.getElementById("adjp-threshold")
  );

  const log2fcThreshold = parseFloat(log2fc?.value ?? "0");
  const adjpThreshold = parseFloat(adjp?.value ?? "0");

  /** @type {Map<string, Set<string>>} */
  const complexes = window.complexome ? window.complexome[0] : new Map();

  /** @type {Map<string, string>} */
  const complexNames = window.complexome ? window.complexome[1] : new Map();

  /** @type {Map<string, [number, number]>} */
  const proteomics = window.userdata ? window.userdata : new Map();

  /** @type {Subunit[]} */
  const perturbedSubunits = complexes
    .entries()
    .flatMap(
      /** @type {([cpxid, members]: [string, Set<string>]) => Iterator<[string, string]>} */ ([
        cpxid,
        members,
      ]) => members.values().map((member) => [cpxid, member]),
    )
    .map(([cpxid, member]) => {
      const [measured_log2fc, measured_adjp] = proteomics.get(member) ?? [
        0,
        Infinity,
      ];
      return {
        name: cpxid,
        subunit: member,
        log2fc: measured_log2fc,
        apvalue: measured_adjp,
      };
    })
    .filter(
      (/** @type {Subunit} */ subunit) =>
        Math.abs(subunit.log2fc) >= log2fcThreshold &&
        subunit.apvalue <= adjpThreshold,
    )
    .toArray();

  const perturbedIDs = new Set(
    perturbedSubunits.map((subunit) => subunit.name),
  );

  const geneNames = await fetchGeneNameMapping(perturbedSubunits);
  const coverage = new Map(
    complexes.entries().map(([cid, members]) => {
      const [num_subunits, measured_subunits] = members
        .values()
        .filter(
          (subunit) =>
            !subunit.includes("CPX-") &&
            !subunit.includes("URS") &&
            !subunit.includes("CHEBI:"),
        )
        .map((subunit) =>
          subunit.includes("-") || subunit.includes("_")
            ? subunit.slice(0, 6)
            : subunit,
        )
        .reduce(
          (acc, subunit) => [
            acc[0] + 1,
            proteomics.has(subunit) ? acc[1] + 1 : acc[1],
          ],
          [0, 0],
        );
      return [cid, measured_subunits / num_subunits];
    }),
  );

  /** @type {Map<string, Perturbation>} */
  const perturbations = new Map(
    complexes
      .entries()
      .filter(([cid, _]) => perturbedIDs.has(cid))
      .map(([cid, members]) => {
        const [regulationScore, count, perturbation] = members
          .values()
          .filter(
            (subunit) =>
              !subunit.includes("CPX-") &&
              !subunit.includes("URS") &&
              !subunit.includes("CHEBI:"),
          )
          .map((subunit) =>
            subunit.includes("-") || subunit.includes("_")
              ? subunit.slice(0, 6)
              : subunit,
          )
          .reduce(
            (acc, subunit) => {
              const proteomicsData = proteomics.get(subunit);
              if (proteomicsData) {
                const [log2fc, adjp] = proteomicsData;
                return [
                  acc[0] + Math.abs(log2fc * -Math.log10(adjp)),
                  acc[1] + 1,
                  whatPerturbation(acc[2], log2fc),
                ];
              } else {
                return acc;
              }
            },
            [0, 0, "Unknown"],
          );
        return [
          cid,
          {
            perturbation,
            score: regulationScore,
            normalizedScore: regulationScore / count,
          },
        ];
      }),
  );

  /** @type {TableRow[]} */
  const tableRows = perturbedSubunits
    .map((/** @type {Subunit} */ subunit) => {
      const perturbation = perturbations.get(subunit.name);
      return {
        cid: subunit.name,
        name: complexNames.get(subunit.name) ?? "",
        coverage: coverage.get(subunit.name) ?? 0.0,
        type: perturbation?.perturbation ?? "",
        score: perturbation?.score ?? 0.0,
        normalizedScore: perturbation?.normalizedScore ?? 0.0,
        subunitID: subunit.subunit,
        geneName: geneNames.get(subunit.subunit) ?? "",
        log2fc: subunit.log2fc,
        adjpval: subunit.apvalue,
      };
    })
    .sort(compare(sorting.column, sorting.order));

  return table(tableRows, viewComplex, sorting);
}

/**
 * Display the Complexome Viewer.
 * @param {string} cid - The Complexome complex identifier.
 */
function viewComplexome(cid) {
  const nameHeader = /** @type {HTMLElement | null} */ (
    /** @type {unknown} */ document.getElementById("complex-name")
  );
  const complexomeLink = /** @type {HTMLAnchorElement | null} */ (
    /** @type {unknown} */ document.getElementById("complexome-link")
  );
  const saveImageButton = /** @type {HTMLButtonElement} */ (
    /** @type {unknown} */ document.getElementById("save-image")
  );
  const iframe = /** @type {HTMLIFrameElement | null} */ (
    /** @type {unknown} */ document.getElementById("my-complexome-viewer")
  );
  const proteomics = window.userdata ? window.userdata : new Map();
  const log2fc = proteomics
    .entries()
    .map(([protein, [log2fc, _]]) => ({ protein, log2fc }))
    .toArray();
  if (iframe) {
    iframe.style.display = "block";
    iframe.contentWindow?.postMessage({ cid, log2fc }, "*");
  }

  if (saveImageButton) {
    saveImageButton.style.display = "block";
  }

  if (nameHeader) {
    nameHeader.innerText = cid;
  }

  if (complexomeLink) {
    complexomeLink.style.display = "block";
    complexomeLink.href = `https://www.ebi.ac.uk/complexportal/complex/${cid}`;
  }
}

function drawComplexomePlots() {
  document
    .getElementById("subunit-dist")
    ?.replaceChildren(...subunitDistributionPlot());
  document
    .getElementById("shared-subunits")
    ?.replaceChildren(...sharedSubunitsPlot());
}

async function drawPlots() {
  document
    .getElementById("proteomics-coverage")
    ?.replaceChildren(...coveragePlot());
  document.getElementById("venn")?.replaceChildren(...vennDiagram());
  document.getElementById("volcano")?.replaceChildren(...volcano());
  document.getElementById("goterms")?.replaceChildren(...goTerms());
  document.getElementById("perturbed-complexes-table")?.replaceChildren(
    ...(await dataTable(viewComplexome, {
      column: "coverage",
      order: "desc",
    })),
  );
}

function clearPlots() {
  document.getElementById("subunit-dist")?.replaceChildren();
  document.getElementById("shared-subunits")?.replaceChildren();
  document.getElementById("proteomics-coverage")?.replaceChildren();
  document.getElementById("venn")?.replaceChildren();
  document.getElementById("volcano")?.replaceChildren();
  document.getElementById("goterms")?.replaceChildren();
  document.getElementById("data-table")?.replaceChildren();
}

/**
 * Handle messages from the worker thread.
 * @param {MessageEvent} event
 */
async function handleMessage(event) {
  if ("userdata" in event.data) {
    controls(true);
    window.userdata = event.data["userdata"];
    if (window.complexome) {
      await drawPlots();
    }
  } else if ("complex" in event.data) {
    window.complexome = event.data["complex"];
    drawComplexomePlots();
    if (window.userdata) {
      await drawPlots();
    }
  } else {
    console.error("Unknown message:", event.data);
  }
}

async function setup() {
  controls(false);
  clearPlots();
  const species = /** @type {HTMLSelectElement | null} */ (
    /** @type {unknown} */ document.querySelector("#species")
  );
  const userfile = /** @type {HTMLInputElement | null} */ (
    /** @type {unknown} */ document.querySelector("#proteomics-file")
  );
  const log2fc = /** @type {HTMLInputElement | null} */ (
    /** @type {unknown} */ document.querySelector("#log2fc-threshold")
  );
  const adjpv = /** @type {HTMLInputElement | null} */ (
    /** @type {unknown} */ document.querySelector("#adjp-threshold")
  );
  const goterms = /** @type {HTMLInputElement | null} */ (
    /** @type {unknown} */ document.querySelector("#top-n-go-terms")
  );

  // setup the proxy worker
  //await registerServiceWorker();

  // setup the CSV parsing worker
  if (window.Worker) {
    const csvWorker = new Worker("/worker/parser.mjs", { type: "module" });
    csvWorker.onmessage = handleMessage;

    // Handle proteomics file changes
    userfile?.addEventListener("change", () => onProteomicsFile(csvWorker));

    // Handle species changes
    species?.addEventListener("change", () =>
      csvWorker.postMessage({ op: "getcomplex", data: species?.value }),
    );

    // Preload cache with the default species
    csvWorker.postMessage({ op: "getcomplex", data: species?.value });

    // If the user has already selected a file, load the data.
    await onProteomicsFile(csvWorker);
  } else {
    console.error(
      "Cannot complete setup: browser does not support web workers.",
    );
    throw Error("Cannot complete setup: browser does not support web workers.");
  }

  log2fc?.addEventListener("change", drawPlots);
  adjpv?.addEventListener("change", drawPlots);

  goterms?.addEventListener("change", () =>
    document.getElementById("goterms")?.replaceChildren(...goTerms()),
  );
}

document.addEventListener("DOMContentLoaded", setup);
document.querySelector("#save-image")?.addEventListener("click", () => {
  const iframe = /** @type {HTMLIFrameElement | null} */ (
    /** @type {unknown} */ document.getElementById("my-complexome-viewer")
  );
  if (iframe) {
    iframe.style.display = "block";
    iframe.contentWindow?.postMessage({ screenshot: true }, "*");
  }
});

window.addEventListener("message", (event) => {
  console.log(event);
  const { action, data } = event.data;
  console.log(action, data);
  switch (action) {
    case "download":
      const { name, url } = data;
      let a = document.createElement("a");
      a.setAttribute("download", `${name}.png`);
      a.setAttribute("href", url);
      console.log("About to download...");
      a.click();
      break;
  }
});

window.addEventListener(
  "custom:table-sort",
  async (/** @type {CustomEvent<Sorting>} */ event) => {
    document.getElementById("perturbed-complexes-table")?.replaceChildren(
      ...(await dataTable(viewComplexome, {
        column: event.detail.column,
        order: event.detail.order,
      })),
    );
  },
);
