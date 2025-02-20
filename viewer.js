import { App } from "https://cdn.skypack.dev/pin/complexviewer@v2.2.4-mcKlBbXJL2XODAKhG35J/mode=imports,min/optimized/complexviewer.js";
import {interpolateRgb} from "https://cdn.skypack.dev/d3-interpolate@3";

function mkScale(protein_log2fc) {
	// Get the range of log2fc values.
	const [min_log2fc, max_log2fc] = protein_log2fc.reduce(
		([cmin, cmax], {log2fc}) => [Math.min(cmin, log2fc), Math.max(cmax, log2fc)],
		[Infinity, -Infinity]
	);

	const scale = max_log2fc === min_log2fc ? 1.0 : (1.0 / (max_log2fc - min_log2fc));
	const interp = interpolateRgb("purple", "orange");

	return (fc) => {
		if (fc === undefined) {
			return "#fff";
		}
		//normalize the fc then get the interpolated colour
		return interp(scale * (fc - min_log2fc));
	}
}

/**
 * @param {HtmlElement}                      el             - The DOM element to draw the viewer in
 * @param {string}                           complex_id     - The complexome complex identifier to draw
 * @param {{protein: str, log2fc: number}[]} protein_log2fc - A mapping of protein names to log2fc values
 * @param {number}                           width          - Width to make the SVG (default 800)
 * @param {number}                           height         - height to make the SVG (default 800)
 */
export async function draw(el, complex_id, protein_log2fc, width, height) {
	const complexviewer = new App(el);
	const theSVG = document.querySelector(".complexViewerSVG");
	if (theSVG) {
		theSVG.style.width = `${width ?? 800}px`;
		theSVG.style.height = `${height ?? 800}px`;
	}
	const res = await fetch(`https://www.ebi.ac.uk/intact/complex-ws/export/${complex_id}`);
	const data = await res.json();
	complexviewer.readMIJSON(data);
	// map proteins to DOM nodes
	const allText = Array.from(el.querySelectorAll("svg text"));
	const labels = new Map();
	data.data.forEach((entry) => {
		if (entry.identifier) {
			labels.set(entry.identifier.id,
								 [
									 allText.find((el) => el.textContent === entry.label)?.parentNode,
									 protein_log2fc.find((protein) => protein.protein === entry.identifier.id)?.log2fc
								 ]
								);
		}
	});
	console.log(labels);
	const colour = mkScale(protein_log2fc);
	labels
		.entries()
		.forEach(([prot, [parent, log2fc]]) => {
			const children = Array.from(parent.childNodes);
			// remove the groups
			children.filter((el) => el.tagName.toLowerCase() === 'g')
				.forEach((el) => parent.removeChild(el));
			// colour
			children[1].setAttribute("fill", colour(log2fc));
		})

	return labels;
}

