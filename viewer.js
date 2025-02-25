import { App } from "https://cdn.skypack.dev/pin/complexviewer@v2.2.4-mcKlBbXJL2XODAKhG35J/mode=imports,min/optimized/complexviewer.js";
import { interpolateRgb } from "https://cdn.skypack.dev/d3-interpolate@3";

function mkScale(protein_log2fc) {
  // Get the range of log2fc values.
  const [min_log2fc, max_log2fc] = protein_log2fc.reduce(
    ([cmin, cmax], { log2fc }) => [
      Math.min(cmin, log2fc),
      Math.max(cmax, log2fc),
    ],
    [Infinity, -Infinity],
  );

  const scale =
    max_log2fc === min_log2fc ? 1.0 : 1.0 / (max_log2fc - min_log2fc);
  const interp = interpolateRgb("purple", "orange");

  return [
    min_log2fc,
    max_log2fc,
    (fc) => {
      if (fc === undefined) {
        return "#fff";
      }
      //normalize the fc then get the interpolated colour
      return interp(scale * (fc - min_log2fc));
    },
  ];
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
  const res = await fetch(
    `https://www.ebi.ac.uk/intact/complex-ws/export/${complex_id}`,
  );
  const data = await res.json();
  complexviewer.readMIJSON(data);

  const allText = Array.from(el.querySelectorAll("svg text"));
  const [, , colour] = mkScale(protein_log2fc);
  data.data.forEach((entry) => {
    if (entry.identifier) {
      const parent = allText.find(
        (el) => el.textContent === entry.label,
      )?.parentNode;
      const log2fc = protein_log2fc.find(
        (protein) => protein.protein === entry.identifier.id,
      )?.log2fc;
      const children = Array.from(parent.childNodes);
      // remove the groups
      children
        .filter((el) => el.tagName.toLowerCase() === "g")
        .forEach((el) => parent.removeChild(el));
      // colour
      children[1].setAttribute("fill", colour(log2fc));
    }
  });
}

/**
 * @param {{protein: str, log2fc: number}[]} protein_log2fc - A mapping of protein names to log2fc values
 * @return {HTMLElement}
 */
export function drawScale(protein_log2fc) {
  const [min, max, colour] = mkScale(protein_log2fc);

  const container = document.createElement("div");
  container.style.width = "100%";
  container.style.maxWidth = "42rem";
  const subContainer = document.createElement("div");
  subContainer.style.position = "relative";
  container.appendChild(subContainer);
  const inner = document.createElement("div");
  inner.style.height = "4rem";
  inner.style.display = "flex";
  subContainer.appendChild(inner);

  protein_log2fc
    .map((val) => val.log2fc)
    .toSorted((a, b) => a - b)
    .forEach((val) => {
      const bar = document.createElement("div");
      bar.style.height = "100%";
      bar.style.flex = "1 1 0%";
      bar.style.backgroundColor = colour(val);
      inner.appendChild(bar);
    });

  const labels = document.createElement("div");
  labels.style.display = "flex";
  labels.style.marginTop = "0.5rem";
  labels.style.justifyContent = "space-between";
  subContainer.appendChild(labels);

  const minLabel = document.createElement("span");
  minLabel.style.fontSize = ".875rem";
  minLabel.style.lineHeight = "1.25rem";
  minLabel.textContent = min.toString();
  labels.appendChild(minLabel);

  const maxLabel = document.createElement("span");
  maxLabel.style.fontSize = ".875rem";
  maxLabel.style.lineHeight = "1.25rem";
  maxLabel.textContent = max.toString();
  labels.appendChild(maxLabel);

  return container;
}
