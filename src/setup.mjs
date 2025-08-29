/**
 * Set up the workflow structure with event listeners.
 * Also run the worker thread to decode CSV documents.
 * And the service worker to cache everything.
 */

import * as d3 from "https://esm.run/d3";

import { minmax, bisect } from "./numeric.mjs";
import { barPlot, histogramPlot, vennPlot, volcanoPlot } from "./plot.mjs";

async function onProteomicsFile(files, csv) {
	if (files && files.length === 1) {
		const text = await files[0].text();
		csv.postMessage({ op: "csv", data: text });
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
	document.getElementById("species").disabled = false;
	document.getElementById("log2fc-threshold").disabled = false;
	document.getElementById("adjp-threshold").disabled = false;
	document.getElementById("top-n-go-terms").disabled = false;
}

function mapGetWithDefault(map, key, dflt) {
	if (map.has(key)) {
		return map.get(key);
	} else {
		return dflt;
	}
}

function subunitDistributionPlot() {
	const count = new Map();

	for (const [complexID, cplx] of window.complexome[0]) {
		let subunits = Array.from(cplx.values())
				.filter((value) => !value.includes("CPX-") && !value.includes("URS") && !value.includes("CHEBI:"));
		count.set(subunits.length, mapGetWithDefault(count, subunits.length, 0) + 1);
	}
	return barPlot(Array.from(count.entries()).sort((a,b) => +a[0] - +b[0]),
								 {
									 hmargin: 30,
									 vmargin: 30,
									 width: 600,
									 height: 400,
									 xlabel: "Number of complexes",
									 ylabel: "↑ Frequency",
									 title: "Subunit distribution (proteins only)",
									 scale: d3.scaleLog,
								 }
								)
}

function sharedSubunitsPlot() {
	const count = new Map();
	const proteinObsCount = new Map();

	for (const [complexID, cplx] of window.complexome[0]) {
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
									 scale: d3.scaleLinear,
									 xticks: Array.from(count.keys()).sort(),
								 }
								)
}

function drawComplexomePlots() {
  document.getElementById("subunit-dist").replaceChildren(...subunitDistributionPlot());
	document.getElementById("shared-subunits").replaceChildren(...sharedSubunitsPlot());
}

function coveragePlot() {
	const coverage = new Map();
	for (const [complexID, cplx] of window.complexome[0]) {
		const subunits = Array.from(cplx.values())
					.filter((subunit) => !subunit.includes("CPX-") && !subunit.includes("URS") && !subunit.includes("CHEBI:"))
					.map((subunit) => (subunit.includes("-") || subunit.includes("_")) ? subunit.slice(6) : subunit);
		const numSubunits = subunits.length;
		const measuredSubunits = subunits
				.filter((subunit) => window.userdata.has(subunit))
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
	for (const [complexID, cplx] of window.complexome[0]) {
		cplx.values()
			.filter((subunit) => !subunit.includes("CPX-") && !subunit.includes("URS") && !subunit.includes("CHEBI:"))
			.map((subunit) => (subunit.includes("-") || subunit.includes("_")) ? subunit.slice(6) : subunit)
			.forEach((subunit) => {
				allCanonicalSubunits.add(subunit);
			});
	}
	return vennPlot({A: allCanonicalSubunits, B: new Set(window.userdata.keys())},
									{
										hmargin: 20,
										vmargin: 20,
										width: 600,
										height: 400,
										alabel: "Complexome proteins",
										blabel: "Proteomics dataset",
									})
}

function volcano() {
	const data = Array.from(window.userdata?.entries().map(([name, [log2fc, adjpval]]) => {
		return {"adjPval": adjpval, "log2fc": log2fc, "name": name};
	}));
	return volcanoPlot(data, {
		hmargin: 20,
		vmargin: 30,
		width: 600,
		height: 400,
		xlabel: "log2 (FC)",
		ylabel: "-log10 (adjPval)",
		title: "Volcano plot",
		log2fcThreshold: document.getElementById("log2fc-threshold").value,
		adjpThreshold: document.getElementById("adjp-threshold").value,
	});
}

function drawPlots() {
	document.getElementById("proteomics-coverage")?.replaceChildren(...coveragePlot());
	document.getElementById("venn")?.replaceChildren(...vennDiagram());
	document.getElementById("volcano")?.replaceChildren(...volcano());
}

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
	console.log('setup() running...')
	// Preload default species
	const species = document.querySelector("#species").value;

	// setup the proxy worker
	//await registerServiceWorker();

	// setup the CSV parsing worker
	if (window.Worker) {
		const csvWorker = new Worker("parser.js", {type: "module"});
		csvWorker.onmessage = handleMessage;
		// Handle proteomics file changes
		document.querySelector("#proteomics-file")?.addEventListener("change", (event) => onProteomicsFile(event.target.files, csvWorker));
		// cache the currently selected species
		csvWorker.postMessage({op: "getcomplex", data: species});
	} else {
		console.error("Cannot complete setup: browser does not support web workers.");
		throw Error("Cannot complete setup: browser does not support web workers.");
	}

	
}

document.addEventListener("DOMContentLoaded", setup);
