// @ts-check
/**
 * Set up the workflow structure with event listeners.
 * Also run the worker thread to decode CSV documents.
 * And the service worker to cache everything.
 */

import { barPlot, histogramPlot, vennPlot, volcanoPlot } from "./plot.mjs";

/**
 * Handle changes to the user supplied proteoimcs file.
 * @param {Event & {currentTarget: HTMLInputElement & EventTarget}} event
 * @param {Worker} csv - A reference to the worker thread for parsing user data.
 */
async function onProteomicsFile(event, csv) {
	if (event.currentTarget.files && event.currentTarget.files.length === 1) {
		const files = event.currentTarget.files;
		const text = await files[0]?.text();
		if (text) {
			csv.postMessage({ op: "csv", data: text });
		}
	}
}

async function registerServiceWorker() {
	if ("serviceWorker" in navigator) {
		try {
			const registration = await navigator.serviceWorker.register("cache.mjs", { scope: "./" });
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

function enableControls() {
	const species = document.getElementById("species");

	const log2fc = document.getElementById("log2fc-threshold");
	const adjp = document.getElementById("adjp-threshold");
	const goterms = document.getElementById("top-n-go-terms");

	if (species && "disabled" in species) {
		species.disabled = false;
	}
	if (log2fc && "disabled" in log2fc) {
		log2fc.disabled = false;
	}
	if (adjp && "disabled" in adjp) {
		adjp.disabled = false;
	}
	if (goterms && "disabled" in goterms) {
		goterms.disabled = false;
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
 * @returns {V}
 */
function mapSetWithDefault(map, key, updater, dflt) {
	const current = map.get(key);
	if (current === undefined) {
		map.set(key, dflt);
		return dflt;
	} else {
		const value = updater(map.get(key) ?? dflt);
		map.set(key, value);
		return value;
	}
}

function subunitDistributionPlot() {
	const count = new Map();
	const complexes = (window.complexome) ? window.complexome[0] : [];

	for (const cplx of complexes.values()) {
		let subunits = Array.from(cplx.values())
			.filter((value) => !value.includes("CPX-") && !value.includes("URS") && !value.includes("CHEBI:"));
		count.set(subunits.length, mapGetWithDefault(count, subunits.length, 0) + 1);
	}
	return barPlot(
		Array.from(count.entries())
			.sort((a, b) => +a[0] - +b[0])
			.map(([x, y]) => [x.toString(), y]),
		{
			hmargin: 30,
			vmargin: 30,
			width: 600,
			height: 400,
			xlabel: "Number of complexes",
			ylabel: "↑ Frequency",
			title: "Subunit distribution (proteins only)",
		}
	)
}

function sharedSubunitsPlot() {
	/** @type {Map<string, number>} */
	const count = new Map();
	const proteinObsCount = new Map();
	const complexes = (window.complexome) ? window.complexome[0] : [];

	for (const cplx of complexes.values()) {
		let subunits = Array.from(cplx.values())
			.filter((value) => !value.includes("CPX-") && !value.includes("URS") && !value.includes("CHEBI:"));
		for (const subunit of subunits) {
			proteinObsCount.set(subunit, mapGetWithDefault(proteinObsCount, subunit, 0) + 1);
		}
	}

	for (const value of proteinObsCount.values()) {
		const key = Math.min(value, 10) <= 9 ? `${value}` : ">10";
		count.set(key, mapGetWithDefault(count, key, 0) + 1);
	}

	return barPlot(Array.from(count.entries()).sort(),
		{
			hmargin: 40,
			vmargin: 30,
			width: 600,
			height: 400,
			xlabel: "Number of protein subunits",
			ylabel: "↑ Frequency",
			title: "Shared protein subunits",
		}
	)
}

function coveragePlot() {
	const coverage = new Map();
	const complexes = (window.complexome) ? window.complexome[0] : [];

	for (const [complexID, cplx] of complexes) {
		const subunits = Array.from(cplx.values())
			.filter((subunit) => !subunit.includes("CPX-") && !subunit.includes("URS") && !subunit.includes("CHEBI:"))
			.map((subunit) => (subunit.includes("-") || subunit.includes("_")) ? subunit.slice(6) : subunit);
		const numSubunits = subunits.length;
		const measuredSubunits = subunits
			.filter((subunit) => window?.userdata?.has(subunit))
			.length;
		coverage.set(complexID, measuredSubunits / numSubunits);
	}

	return histogramPlot(Array.from(coverage.values()),
		{
			hmargin: 20,
			vmargin: 20,
			width: 600,
			height: 400,
			xlabel: "Proteomics coverage",
			ylabel: "↑ Frequency",
			title: "",
		}
	)
}

function vennDiagram() {
	const allCanonicalSubunits = new Set();
	const complexes = (window.complexome) ? window.complexome[0] : [];

	for (const cplx of complexes.values()) {
		cplx.values()
			.filter((subunit) => !subunit.includes("CPX-") && !subunit.includes("URS") && !subunit.includes("CHEBI:"))
			.map((subunit) => (subunit.includes("-") || subunit.includes("_")) ? subunit.slice(6) : subunit)
			.forEach((subunit) => {
				allCanonicalSubunits.add(subunit);
			});
	}
	return vennPlot({ A: allCanonicalSubunits, B: new Set(window?.userdata?.keys()) },
		{
			hmargin: 20,
			vmargin: 20,
			width: 600,
			height: 400,
			alabel: "Complexome proteins",
			blabel: "Proteomics dataset",
			title: "",
		})
}

function volcano() {
	const data = Array.from(window?.userdata?.entries().map(([name, [log2fc, adjpval]]) => {
		return { "adjPval": adjpval, "log2fc": log2fc, "name": name };
	}) ?? []);
	const log2fc = /** @type {HTMLInputElement | null} */ (/** @type {unknown} */ document.getElementById("log2fc-threshold"));
	const adjp = /** @type {HTMLInputElement | null} */ (/** @type {unknown} */ document.getElementById("adjp-threshold"));

	return volcanoPlot(data, {
		hmargin: 20,
		vmargin: 30,
		width: 600,
		height: 400,
		xlabel: "log2 (FC)",
		ylabel: "-log10 (adjPval)",
		title: "Volcano plot",
		log2fcThreshold: parseFloat(log2fc?.value ?? "0"),
		adjpThreshold: parseFloat(adjp?.value ?? "0"),
	});
}

/**
 * @typedef SubunitInfo Data applicable to complex subunits (proteins).
 * @property cid {string} The complex identifier
 * @property name {string} The name of the complex
 * @property subunit {string} The name of the subunit
 * @property log2fc {number} The log2(Fold change)
 * @property apvalue {number} The adjusted p-value
 */

function goTerms() {
	const log2fc = /** @type {HTMLInputElement | null} */ (/** @type {unknown} */ document.getElementById("log2fc-threshold"));
	const adjp = /** @type {HTMLInputElement | null} */ (/** @type {unknown} */ document.getElementById("adjp-threshold"));
	const numTermsEl = /** @type {HTMLInputElement | null} */ (/** @type {unknown} */ document.getElementById("top-n-go-terms"));

	const numTerms = parseInt(numTermsEl?.value ?? "10");

	const log2fcThreshold = parseFloat(log2fc?.value ?? "0");
	const adjpThreshold = parseFloat(adjp?.value ?? "0");

	/** @type {Map<string, Set<string>>} */
	const complexes = (window.complexome) ? window.complexome[0] : new Map();

	/** @type {Map<string, string>} */
	const names = (window.complexome) ? window.complexome[1] : new Map();

	/** @type {Map<string, [number, number]>} */
	const proteomics = (window.userdata) ? window.userdata : new Map();

	/** @type {Set<string>} */
	const perturbed = new Set();

	/** @type {Array<[string, number]>} */
	const counts = complexes.values().flatMap((/** @type {string} */ member) => {
		const [measured_log2fc, measured_adjp] = proteomics.get(member) ?? [0, Infinity];
		return [member, measured_log2fc, measured_adjp];
	}).filter(([_, measured_log2fc, measured_adjp]) => (measured_log2fc >= log2fcThreshold) && (measured_adjp <= adjpThreshold))
				.reduce((acc, [member, _x, _y]) => mapSetWithDefault(acc, member, (count) => count + 1, 1), new Map()).entries().toArray().sort((a, b) => +a[1] - +b[1]).slice(numTerms - 1);

	return barPlot(counts, {
		hmargin: 20,
		vmargin: 20,
		width: 600,
		height: 400,
		xlabel: "",
		ylabel: "↑ Frequency",
		title: "Most frequent GO terms associated with the perturbed complexes",
	});
}

function drawComplexomePlots() {
	document.getElementById("subunit-dist")?.replaceChildren(...subunitDistributionPlot());
	document.getElementById("shared-subunits")?.replaceChildren(...sharedSubunitsPlot());
}

function drawPlots() {
	document.getElementById("proteomics-coverage")?.replaceChildren(...coveragePlot());
	document.getElementById("venn")?.replaceChildren(...vennDiagram());
	document.getElementById("volcano")?.replaceChildren(...volcano());
	document.getElementById("goterms")?.replaceChildren(...goTerms());
}

/**
 * Handle messages from the worker thread.
 * @param {MessageEvent} event
 */
function handleMessage(event) {
	if ("userdata" in event.data) {
		enableControls();
		window.userdata = event.data["userdata"];
		if (window.complexome) {
			drawPlots();
		}
	} else if ("complex" in event.data) {
		console.log(`Got data from complex: ${event.data['complex']}`)
		window.complexome = event.data["complex"];
		drawComplexomePlots();
		if (window.userdata) {
			drawPlots();
		}
	} else {
		console.error("Unknown message:", event.data);
	}
}

async function setup() {
	// Preload default species
	const species = /** @type {HTMLInputElement | null} */ (/** @type {unknown} */ document.querySelector("#species"));
	const userfile = /** @type {HTMLInputElement | null} */ (/** @type {unknown} */ document.querySelector("#proteomics-file"));

	// setup the proxy worker
	//await registerServiceWorker();

	// setup the CSV parsing worker
	if (window.Worker) {
		const csvWorker = new Worker("/worker/parser.mjs", { type: "module" });
		csvWorker.onmessage = handleMessage;
		// Handle proteomics file changes
		userfile?.addEventListener("change", (event) => onProteomicsFile((/** @type {Event & {currentTarget: HTMLInputElement & EventTarget}} */ (/** @type {unknown} */ event)), csvWorker));
		// cache the currently selected species
		csvWorker.postMessage({ op: "getcomplex", data: species?.value });
	} else {
		console.error("Cannot complete setup: browser does not support web workers.");
		throw Error("Cannot complete setup: browser does not support web workers.");
	}
}

document.addEventListener("DOMContentLoaded", setup);
