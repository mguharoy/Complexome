/**
 * Set up the workflow structure with event listeners.
 * Also run the worker thread to decode CSV documents.
 */

//import Plotly from "https://cdn.skypack.dev/plotly.js-dist";
import * as d3 from "https://esm.run/d3";

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

/**
 * @param {Array<number>} data - the data to summarize
 * @return [number, number]    - [minimum value, maximum value]
 */
function minmax(data) {
	return data.reduce(([min, max], value) => [Math.min(min, value), Math.max(max, value)], [Infinity, -Infinity]);
}

function histogramPlot(data, options) {
	let [min, max] = minmax(data);
	const numBins = (min === max) ? 1 : Math.min(data.length, Math.max(10, Math.floor(data.length / 160)));

	const x = d3.scaleLinear()
				.domain([0, 1.0])
				.range([2*options.hmargin, options.width - options.hmargin]);
	const histogram = d3.histogram()
				.value((d) => d)
				.domain(x.domain())
				.thresholds(x.ticks(numBins));

	const bins = histogram(data).reduce((acc, val) => {
		const [bin_min, bin_max] = minmax(val);
		if (acc.length === 0) {
			return [
				{
					from: 0,
					to: Math.max(bin_max, 0.01),
					height: val.length,
				}
			];
		} else {
			return [...acc, {
				from: acc[acc.length - 1].to,
				to: bin_max,
				height: val.length,
			}];
		}
	}, []);
	const maxHeight = d3.max(bins, (d) => d.height);

	const y = d3.scaleLinear()
				.domain([0, d3.max(bins, (d) => d.height)])
				.range([options.height - 2*options.vmargin, options.vmargin])

	const svg = d3.create("svg")
				.attr("width", options.width)
				.attr("height", options.height)
				.attr("viewBox", [0, 0, options.width, options.height])
				.attr("style", "max-width: 100%; height: auto;");
	const tooltip = d3.create("div").classed("tooltip", true);

	svg.append("g")
		.attr("transform", `translate(0, ${options.height - options.vmargin - 10})`)
		.call(d3.axisBottom(x))
		.call((g) => g.append("text")
					.attr("x", options.width / 2)
					.attr("y", 25)
					.attr("fill", "currentColor")
					.attr("text-anchor", "middle")
					.text(options.xlabel));
	svg.append("g")
		.attr("transform", `translate(${2*options.hmargin}, ${options.vmargin - 10})`)
		.call(d3.axisLeft(y))
		.call(g => g.append("text")
          .attr("x", -options.hmargin)
          .attr("y", options.vmargin - 10)
          .attr("fill", "currentColor")
          .attr("text-anchor", "start")
          .text(options.ylabel));
	svg.append("g")
		.selectAll()
		.data(bins)
		.join("rect")
		.attr("class", "bar")
		.attr("x", 1)
		.attr("transform", (d) => `translate(${x(d.from)}, ${y(d.height) + options.vmargin - 10})`)
		.attr("width", (d) => x(d.to) - x(d.from))
		.attr("height", (d) => options.height - y(d.height) - 2*options.vmargin)
		.on("mouseover", (event, d) => {
			console.log(d.height, max / 5);
			tooltip.transition().duration(200).style("opacity", 1);
			tooltip.html(`${d.height}`)
				.style("left", `${x(d.from)}px`)
				.style("top", `${y(Math.max(d.height, maxHeight / 10)) + 10}px`);
		})
		.on("mouseout", () => tooltip.transition().duration(200).style("opacity", 0));

	svg.append("text")
		.attr("x", options.width / 2)
		.attr("y", options.vmargin)
		.attr("text-anchor", "middle")
		.style("font-size", "1.5em")
		.text(options.title);

	return [svg.node(), tooltip.node()];
}

/**
 * @param {Array<[x, y]>} data - the data to plot
 * @param {Object} options
 * @param {string} options.parent
 * @param {number} options.hmargin
 * @param {number} options.vmargin
 * @param {number} options.width
 * @param {number} options.height
 * @param {string} options.title
 * @param {string} options.xlabel
 * @param {string} options.ylabel
 * @param {d3.scale} options.scale
 * @param {string[] | undefined} options.xticks
 */
function barPlot(data, options) {
	const x = d3.scaleBand()
				.domain(data.map((d) => d[0]))
				.range([options.hmargin, options.width - options.hmargin])
				.padding(0.1);
	// const y = d3.scaleLinear()
	// 			.domain([0, d3.max(data, (d) => d[1])])
	// 			.range([options.height - options.margin, options.margin]);
	const max = d3.max(data, (d) => d[1]);
  const y = options.scale()
				.domain([1, max])
				.range([options.height - options.vmargin, options.vmargin]);
	const svg = d3.create("svg")
				.attr("width", options.width)
				.attr("height", options.height)
				.attr("viewBox", [0, 0, options.width, options.height])
				.attr("style", "max-width: 100%; height: auto;");
	const tooltip = d3.create("div").classed("tooltip", true);

	const xaxis = options.xticks ?
				d3.axisBottom(x).tickValues(options.xticks) :
				d3.axisBottom(x).tickValues(x.domain().filter((_, d) => (d+1) % 5 === 0));

	svg.append("g")
    .selectAll()
    .data(data)
    .join("rect")
		.attr("class", "bar")
    .attr("x", (d) => x(d[0]))
    .attr("y", (d) => y(d[1]))
    .attr("height", (d) => {return y(1) - y(d[1]);})
    .attr("width", x.bandwidth())
		.on("mouseover", (event, d) => {
			tooltip.transition().duration(200).style("opacity", 1);
			tooltip.html(`${d[1]}`)
				.style("left", `${x(d[0])}px`)
				.style("top", `${y(Math.max(d[1], max / 10))}px`);
		})
		.on("mouseout", () => tooltip.transition().duration(200).style("opacity", 0));

	svg.append("g")
    .attr("transform", `translate(0, ${options.height - options.vmargin})`)
    .call(xaxis)
		.call((g) => g.append("text")
					.attr("x", options.width / 2)
					.attr("y", 25)
					.attr("fill", "currentColor")
					.attr("text-anchor", "middle")
					.text(options.xlabel));

  svg.append("g")
    .attr("transform", `translate(${options.hmargin}, 0)`)
    .call(d3.axisLeft(y))
		.call(g => g.append("text")
          .attr("x", -options.hmargin)
          .attr("y", options.vmargin - 10)
          .attr("fill", "currentColor")
          .attr("text-anchor", "start")
          .text(options.ylabel));

	svg.append("text")
		.attr("x", options.width / 2)
		.attr("y", options.vmargin)
		.attr("text-anchor", "middle")
		.style("font-size", "1.5em")
		.text(options.title);

  return [svg.node(), tooltip.node()];
}

function vennPlot(data, options) {
	const intersection = data.A.intersection(data.B);
	const total = data.A.size + data.B.size;
	const [Ab, aB, AB] = [data.A.size / total, data.B.size / total, intersection.size / total];
	const [R, r] = [Math.sqrt(Aa / Math.PI), Math.sqrt(Ab / Math.PI)];
	if (intersection.size > 0) {
		
	} else {
		const d = R + r + Math.max((r + R / 2), 0.2);
	}
	const svg = d3.create("svg")
				.attr("width", options.width)
				.attr("height", options.height)
				.attr("viewBox", [0, 0, options.width, options.height])
				.attr("style", "max-width: 100%; height: auto;");
	const tooltip = d3.create("div").classed("tooltip", true);
	return [svg.node(), tooltip.node()];
}

function mapGetWithDefault(map, key, dflt) {
	if (map.has(key)) {
		return map.get(key);
	} else {
		return dflt;
	}
}

function subunitDistributionPlot(selector) {
	const count = new Map();

	for (const [complexID, cplx] of window.complexome[0]) {
		let subunits = Array.from(cplx.values())
				.filter((value) => !value.includes("CPX-") && !value.includes("URS") && !value.includes("CHEBI:"));
		count.set(subunits.length, mapGetWithDefault(count, subunits.length, 0) + 1);
	}
	return barPlot(Array.from(count.entries()).sort((a,b) => +a[0] - +b[0]),
								 {
									 parent: selector,
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

function sharedSubunitsPlot(selector) {
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
								 {parent: selector,
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
  document.getElementById("subunit-dist").replaceChildren(...subunitDistributionPlot("#subunit-dist"));
	document.getElementById("shared-subunits").replaceChildren(...sharedSubunitsPlot("#shared-subunits"));
}

function coveragePlot(selector) {
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
								 {parent: selector,
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

function drawPlots() {
	document.getElementById("proteomics-coverage").replaceChildren(...coveragePlot("#proteomics-coverage"));
}

function handleMessage(event) {
	if ("userdata" in event.data) {
		enableControls();
		window.userdata = event.data["userdata"];
		if (window.complexome) {
			drawPlots();
		}
	} else if ("complex" in event.data) {
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
