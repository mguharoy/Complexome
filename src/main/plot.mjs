// @ts-check
/**
 * Plotting functions
 */

import * as d3 from "https://esm.run/d3";

import { minmax, bisect, identity } from "../shared/numeric.mjs";

/**
 * @typedef TableRow A row from the data table.
 * @property cid {string}
 * @property name {string}
 * @property coverage {number}
 * @property type {string}
 * @property score {number}
 * @property normalizedScore {number}
 * @property subunitID {string}
 * @property geneName {string}
 * @property log2fc {number}
 * @property adbpval {number}
 */

/**
 * @typedef PlotOptions Configure how the plot should be displayed
 *
 * @property margin {[number, number, number, number]} The margins: top, right, bottom, left
 * @property width {number} The total width of the plot image (in px).
 * @property height {number} The total height of the plot image (in px).
 * @property title {string} The title text.
 * @property xlabel {string} The x-axis label text.
 * @property ylabel {string} The y-axis label text.
 */

/**
 * @typedef Interval A histogram bin.
 *
 * @property from {number} The start (left side) of the bin.
 * @property to {number} The end (right side) of the bin.
 * @property height {number} The height (count) of the bin.
 */

/**
 * Draw a histogram from numeric data.
 * @param {number[]} data - The data to draw the histogram from.
 * @param {PlotOptions} options - Configure the plot visuals.
 * @returns {Node[]} - Return the plot and tooltip DOM nodes.
 */
function histogramPlot(data, options) {
  const [mTop, mRight, mBottom, mLeft] = options.margin;
  const [min, max] = minmax(data);
  const numBins =
    min === max
      ? 1
      : Math.min(data.length, Math.max(10, Math.floor(data.length / 160)));

  const x = d3
    .scaleLinear()
    .domain([0, 1.0])
    .range([mLeft, options.width - mRight]);

  /** @type {Array<d3.Bin<number, number>>} */
  const binnedData = d3
    .bin()
    .value((/** @type {number} */ d) => d)
    .domain([0, 1.0])
    .thresholds(x.ticks(numBins))(data);

  /** @type {Interval[]} */
  const bins = binnedData.reduce(
    (/** @type {Array<Interval>} */ acc, /** @type {number[]} */ val) => {
      if (val.length === 0) return acc;
      const [bin_min, bin_max] = minmax(val);
      return [
        ...acc,
        {
          from: acc[acc.length - 1]?.to ?? bin_min,
          to: Math.max(bin_max, 0.01),
          height: val.length,
        },
      ];
    },
    [],
  );

  /** @type {number} */
  const maxHeight = d3.max(bins, (/** @type {Interval} */ d) => d.height) ?? 0;

  const y = d3
    .scaleLinear()
    .domain([0, maxHeight])
    .range([options.height - options.margin[2], options.margin[0]]);

  const svg = d3
    .create("svg")
    .attr("width", options.width)
    .attr("height", options.height)
    .attr("viewBox", [0, 0, options.width, options.height])
    .attr("style", "max-width: 100%; height: auto;");
  const tooltip = d3.create("div").classed("tooltip", true);

  svg
    .append("g")
    .selectAll()
    .data(bins)
    .join("rect")
    .attr("class", "bar")
    .attr("x", 1)
    .attr(
      "transform",
      /** @param {Interval} d */ (/** @type {Interval} */ d) =>
        `translate(${x(d.from)}, ${y(d.height)})`,
    )
    .attr("width", (/** @type {Interval} */ d) => x(d.to) - x(d.from))
    .attr(
      "height",
      (/** @type {Interval} */ d) =>
        options.height - y(d.height) - options.margin[0],
    )
    .on("mouseover", (/** @type {unknown} */ _, /** @type {Interval} */ d) => {
      tooltip.transition().duration(200).style("opacity", 1);
      tooltip
        .html(`${d.height}`)
        .style("left", `${x(d.from)}px`)
        .style("top", `${y(Math.max(d.height, maxHeight / 10)) + 10}px`);
    })
    .on("mouseout", () =>
      tooltip.transition().duration(200).style("opacity", 0),
    );
  svg
    .append("g")
    .attr("transform", `translate(0, ${options.height - mBottom})`)
    .call(d3.axisBottom(x))
    .call((g) =>
      g
        .append("text")
        .attr("x", options.width / 2)
        .attr("y", mBottom)
        .attr("fill", "currentColor")
        .attr("text-anchor", "middle")
        .text(options.xlabel),
    );
  svg
    .append("g")
    .attr("transform", `translate(${mLeft}, 0)`)
    .call(d3.axisLeft(y))
    .call((g) =>
      g
        .append("text")
        .attr("x", -mLeft)
        .attr("y", mTop - 10)
        .attr("fill", "currentColor")
        .attr("text-anchor", "start")
        .text(options.ylabel),
    );

  svg
    .append("text")
    .attr("x", options.width / 2)
    .attr("y", options.margin[0])
    .attr("text-anchor", "middle")
    .style("font-size", "1.5em")
    .text(options.title);

  const result = [svg.node(), tooltip.node()];
  return result[0] && result[1] ? [result[0], result[1]] : [];
}

/**
 * @param {Array<[string, number]>} data - the data to plot shaped as [x, y].
 * @param {PlotOptions} options - Configure the plot display.
 * @returns {Node[]} The plot SVG node and the tooltip node.
 */
function barPlot(data, options) {
  const [mTop, mRight, mBottom, mLeft] = options.margin;
  const x = d3
    .scaleBand()
    .domain(data.map((d) => d[0]))
    .range([mLeft, options.width - mRight])
    .padding(0.1);
  const max = d3.max(data, (/** @type {[unknown, number]}  */ d) => d[1]) ?? 0;
  const y = d3
    .scaleLinear()
    .domain([0, max])
    .range([options.height - mBottom, mTop]);
  const svg = d3
    .create("svg")
    .attr("width", options.width)
    .attr("height", options.height)
    .attr("viewBox", [0, 0, options.width, options.height])
    .attr("style", "max-width: 100%; height: auto;");
  const tooltip = d3.create("div").classed("tooltip", true);

  svg
    .append("g")
    .selectAll()
    .data(data)
    .join("rect")
    .attr("class", "bar")
    .attr("x", (/** @type {[string, unknown]} */ d) => x(d[0]) ?? "0")
    .attr("y", (/** @type {[unknown, number]} */ d) => y(d[1]))
    .attr("height", (/** @type {[unknown, number]} */ d) => {
      return y(0) - y(d[1]);
    })
    .attr("width", x.bandwidth())
    .on(
      "mouseover",
      (/** @type {unknown} */ _, /** @type {[string, number]} */ d) => {
        tooltip.transition().duration(200).style("opacity", 1);
        tooltip
          .html(`${d[1]}`)
          .style("left", `${x(d[0])}px`)
          .style("top", `${y(Math.max(d[1], max / 10))}px`);
      },
    )
    .on("mouseout", () =>
      tooltip.transition().duration(200).style("opacity", 0),
    );

  svg
    .append("g")
    .attr("transform", `translate(0, ${options.height - mBottom})`)
    .call(d3.axisBottom(x))
    .selectAll("text")
    .style("text-anchor", "end")
    .attr("dx", "-.8em")
    .attr("dy", ".15em")
    .attr("transform", "rotate(-65)")
    .call((g) =>
      g
        .append("text")
        .attr("x", options.width / 2)
        .attr("y", 25)
        .attr("fill", "currentColor")
        .attr("text-anchor", "middle")
        .text(options.xlabel),
    );

  svg
    .append("g")
    .attr("transform", `translate(${mLeft}, 0)`)
    .call(d3.axisLeft(y))
    .call((g) =>
      g
        .append("text")
        .attr("x", -mLeft)
        .attr("y", mTop - 10)
        .attr("fill", "currentColor")
        .attr("text-anchor", "start")
        .text(options.ylabel),
    );

  svg
    .append("text")
    .attr("x", options.width / 2)
    .attr("y", mTop)
    .attr("text-anchor", "middle")
    .style("font-size", "1.5em")
    .text(options.title);

  const result = [svg.node(), tooltip.node()];
  return result[0] && result[1] ? [result[0], result[1]] : [];
}

/**
 * @typedef VennData Data shape to supply the vennPlot function.
 *
 * @property A {Set<unknown>} The data in the left set.
 * @property B {Set<unknown>} The data in the right set.
 */

/**
 * @typedef VennPlotOptions Special properties for Venn diagram plots.
 *
 * @property alabel {string} The label for the left set.
 * @property blabel {string} The label for the right set.
 */

/**
 * Plot a Venn Diagram.
 * wikipedia.org/wiki/Venn_diagram
 *
 * @param {VennData} data - Data to plot.
 * @param {Omit<PlotOptions, "xlabel" | "ylabel"> & VennPlotOptions} options - Configure the plot visuals.
 * @returns {Node[]}
 */
function vennPlot(data, options) {
  const [mTop, mRight, mBottom, mLeft] = options.margin;
  const intersection = data.A.intersection(data.B);
  const total = data.A.size + data.B.size;
  const [Ab, aB, AB] = [
    data.A.size / total,
    data.B.size / total,
    intersection.size / total,
  ];

  const [R, r] = [Math.sqrt(Ab / Math.PI), Math.sqrt(aB / Math.PI)];
  const scale = Math.min(
    (options.height - (mTop + mBottom)) / (2 * Math.max(R, r)),
    (options.width - (mLeft + mRight)) / (2 * Math.max(R, r)),
  );

  /** @type {[number, number]} */
  const bracket = AB < aB ? [R, 1] : [0.001, R];
  const d =
    intersection.size > 0
      ? data.B.isSubsetOf(data.A)
        ? Math.abs(R - r) + 0.0001
        : bisect(R, r, AB, bracket)
      : R + r + Math.max(r + R / 2, 0.2);
  const xoff = options.width / 2 - (d / 2) * scale;

  const xintReal = -(r * r - R * R - d * d) / (2 * d);
  const xint = xintReal * scale + xoff;
  const y1int =
    Math.sqrt(R * R - xintReal * xintReal) * scale + options.height / 2;
  const y2int =
    -Math.sqrt(R * R - xintReal * xintReal) * scale + options.height / 2;

  const svg = d3
    .create("svg")
    .attr("width", options.width)
    .attr("height", options.height)
    .attr("viewBox", [0, 0, options.width, options.height])
    .attr("style", "max-width: 100%; height: auto;");
  const tooltip = d3.create("div").classed("tooltip", true);

  svg
    .append("path")
    .attr("class", "complexome-venn")
    .attr(
      "d",
      `M ${xint} ${y2int} A ${R * scale} ${R * scale} 0 1 0 ${xint} ${y1int} A ${r * scale} ${r * scale} 0 0 1 ${xint} ${y2int}Z`,
    )
    .on("mouseenter", (/** @type {MouseEvent} */ event) => {
      const [mx, my] = d3.pointer(event);
      tooltip.transition().duration(200).style("opacity", 1);
      tooltip
        .html((data.A.size - intersection.size).toString())
        .style("left", `${mx}px`)
        .style("top", `${my - 22}px`);
    })
    .on("mousemove", (/** @type {MouseEvent} */ event) => {
      const [mx, my] = d3.pointer(event);
      tooltip.style("left", `${mx}px`).style("top", `${my - 22}px`);
    })
    .on("mouseout", () =>
      tooltip.transition().duration(200).style("opacity", 0),
    );
  svg
    .append("path")
    .attr("class", "intersection-venn")
    .attr(
      "d",
      `M ${xint} ${y2int} A ${r * scale} ${r * scale} 0 0 0 ${xint} ${y1int} A ${R * scale} ${R * scale} 0 0 0 ${xint} ${y2int}Z`,
    )
    .on("mouseenter", (/** @type {MouseEvent} */ event) => {
      const [mx, my] = d3.pointer(event);
      tooltip.transition().duration(200).style("opacity", 1);
      tooltip
        .html(intersection.size.toString())
        .style("left", `${mx}px`)
        .style("top", `${my - 22}px`);
    })
    .on("mousemove", (/** @type {MouseEvent} */ event) => {
      const [mx, my] = d3.pointer(event);
      tooltip.style("left", `${mx}px`).style("top", `${my - 22}px`);
    })
    .on("mouseout", () =>
      tooltip.transition().duration(200).style("opacity", 0),
    );
  svg
    .append("path")
    .attr("class", "proteomics-venn")
    .attr(
      "d",
      `M ${xint} ${y1int} A ${r * scale} ${r * scale} 0 1 0 ${xint} ${y2int} A ${R * scale} ${R * scale} 0 0 1 ${xint} ${y1int}Z`,
    )
    .on("mouseenter", (/** @type {MouseEvent} */ event) => {
      const [mx, my] = d3.pointer(event);
      tooltip.transition().duration(200).style("opacity", 1);
      tooltip
        .html((data.B.size - intersection.size).toString())
        .style("left", `${mx}px`)
        .style("top", `${my - 22}px`);
    })
    .on("mousemove", (/** @type {MouseEvent} */ event) => {
      const [mx, my] = d3.pointer(event);
      tooltip.style("left", `${mx}px`).style("top", `${my - 22}px`);
    })
    .on("mouseout", () =>
      tooltip.transition().duration(200).style("opacity", 0),
    );
  svg
    .append("text")
    .attr("class", "venn-label")
    .attr("x", options.width / 3)
    .attr("y", (4 * options.height) / 5)
    .attr("fill", "currentColor")
    .attr("text-anchor", "end")
    .text(options.alabel);
  svg
    .append("text")
    .attr("class", "venn-label")
    .attr("x", (2 * options.width) / 3)
    .attr("y", (4 * options.height) / 5)
    .attr("fill", "currentColor")
    .attr("text-anchor", "start")
    .text(options.blabel);

  const result = [svg.node(), tooltip.node()];
  return result[0] && result[1] ? [result[0], result[1]] : [];
}

/**
 * @typedef VolcanoPlotOptions Configure the volcano plot.
 *
 * @property log2fcThreshold {number} Vertical lines to plot.
 * @property adjpThreshold {number} A horizonal line.
 */

/**
 * @typedef VolcanoDatum A single observation to plot.
 *
 * @property adjPval {number}  The adjusted P-value.
 * @property log2fc {number} log2(fold change).
 * @property name {string} The protein name.
 */

/**
 * Draw a Volcano plot.
 * https://en.wikipedia.org/wiki/Volcano_plot_(statistics)
 *
 * @param {VolcanoDatum[]} data - The data to plot.
 * @param {PlotOptions & VolcanoPlotOptions} options - Configure the plot visuals.
 * @returns {Node[]} - Return the plot and the tooltip DOM nodes.
 */
function volcanoPlot(data, options) {
  const [mTop, mRight, mBottom, mLeft] = options.margin;
  const [xmin, xmax] = minmax(
    data,
    (/** @type {VolcanoDatum} */ datum) => datum.adjPval,
  );
  const [ymin, ymax] = minmax(
    data,
    (/** @type {VolcanoDatum} */ datum) => datum.log2fc,
  );
  const xScale = d3
    .scaleLinear()
    .domain([Math.min(-options.log2fcThreshold, ymin), ymax])
    .range([mLeft, options.width - mRight])
    .nice();
  const yScale = d3
    .scaleLinear()
    .domain([-Math.log10(xmin), Math.min(-Math.log10(xmax), 0.0)])
    .range([mTop, options.height - mBottom])
    .nice();

  /** @type {[number, number]} */
  const [countBlue, countRed] = data.reduce(
    (acc, datum) => {
      if (-Math.log10(datum.adjPval) > -Math.log10(options.adjpThreshold)) {
        if (datum.log2fc > options.log2fcThreshold) {
          return [acc[0], acc[1] + 1];
        } else if (datum.log2fc < -options.log2fcThreshold) {
          return [acc[0] + 1, acc[1]];
        }
      }
      return acc;
    },
    [0, 0],
  );

  /** @type {(log2fc: number, pval: number) => string} */
  const colour = (log2fc, pval) => {
    if (-Math.log10(pval) > -Math.log10(options.adjpThreshold)) {
      if (log2fc > options.log2fcThreshold) {
        return "#cc0000";
      } else if (log2fc < -options.log2fcThreshold) {
        return "#003399";
      } else {
        return "#8a8a8a";
      }
    } else {
      return "#8a8a8a";
    }
  };

  const svg = d3
    .create("svg")
    .attr("width", options.width)
    .attr("height", options.height)
    .attr("viewBox", [0, 0, options.width, options.height])
    .attr("style", "max-width: 100%; height: auto;");
  const tooltip = d3.create("div").classed("tooltip", true);

  svg
    .append("g")
    .attr("transform", `translate(0, ${options.height - mBottom})`)
    .call(d3.axisBottom(xScale))
    .call((g) =>
      g
        .append("text")
        .attr("x", options.width / 2)
        .attr("y", 25)
        .attr("fill", "currentColor")
        .attr("text-anchor", "middle")
        .text(options.xlabel),
    );
  svg
    .append("g")
    .attr("transform", `translate(${mLeft}, 0)`)
    .call(d3.axisLeft(yScale))
    .call((g) =>
      g
        .append("text")
        .attr("x", -mLeft)
        .attr("y", 15)
        .attr("fill", "currentColor")
        .attr("text-anchor", "start")
        .text(options.ylabel),
    );

  svg
    .append("g")
    .selectAll()
    .data(data)
    .join("circle")
    .attr("cx", (datum) => xScale(datum.log2fc))
    .attr("cy", (datum) => yScale(-Math.log10(datum.adjPval)))
    .attr("r", 3)
    .attr("opacity", 0.5)
    .attr("fill", (datum) => colour(datum.log2fc, datum.adjPval))
    .on(
      "mouseenter",
      (/** @type {MouseEvent} */ event, /** @type {VolcanoDatum} */ datum) => {
        const [mx, my] = d3.pointer(event);
        tooltip.transition().duration(200).style("opacity", 1);
        tooltip
          .html(
            `${datum.name}<br/>log2(FC): ${datum.log2fc.toFixed(3)}<br/>-log10(pval): ${-Math.log10(datum.adjPval).toFixed(3)}`,
          )
          .style("left", `${mx + 5}px`)
          .style("top", `${my - 46}px`);
      },
    )
    .on("mouseout", () =>
      tooltip.transition().duration(200).style("opacity", 0),
    );

  // The threshold markers
  const yThresh = yScale(-Math.log10(options.adjpThreshold));
  svg
    .append("polyline")
    .attr(
      "points",
      `${mLeft},${yThresh} ${xScale(-options.log2fcThreshold)},${yThresh} ${xScale(-options.log2fcThreshold)},${mTop}`,
    )
    .attr("fill", "none")
    .attr("stroke", "red")
    .attr("stroke-width", "1px")
    .attr("stroke-dasharray", "4");

  svg
    .append("polyline")
    .attr(
      "points",
      `${options.width - mRight},${yThresh} ${xScale(options.log2fcThreshold)},${yThresh} ${xScale(options.log2fcThreshold)},${mTop}`,
    )
    .attr("fill", "none")
    .attr("stroke", "red")
    .attr("stroke-width", "1px")
    .attr("stroke-dasharray", "4");

  // Count labels
  svg
    .append("text")
    .attr("x", options.width / 10)
    .attr("y", options.height / 5)
    .attr("text-anchor", "start")
    .attr("font-size", "0.5em")
    .attr("font-weight", "bold")
    .attr("fill", "#003399")
    .text(`${countBlue}`);
  svg
    .append("text")
    .attr("x", (9 * options.width) / 10)
    .attr("y", options.height / 5)
    .attr("text-anchor", "end")
    .attr("font-size", "0.5em")
    .attr("font-weight", "bold")
    .attr("fill", "#cc0000")
    .text(`${countRed}`);

  svg
    .append("text")
    .attr("x", options.width / 2)
    .attr("y", mTop)
    .attr("text-anchor", "middle")
    .style("font-size", "1.5em")
    .text(options.title);

  const result = [svg.node(), tooltip.node()];
  return result[0] && result[1] ? [result[0], result[1]] : [];
}

/**
 * Render data rows for a table.
 * @param body {d3.Selection<HTMLTableSectionElement, undefined, null, undefined>}
 * @param data {TableRow[]}
 */
function renderTableData(body, data) {
	const rows = body.selectAll("tr").data(data).enter().append("tr");

  rows
    .selectAll("td")
    .data((/** @type {TableRow} */ row) => Object.entries(row))
    .enter()
    .append("td")
		.style("text-align", (/** @type {[string, string | number]} */ [_, data]) => {
			if (typeof data === "string") {
				return "left";
			} else {
				return "right";
			}
		})
    .text(([key, data]) => {
			if (key === "coverage" && (typeof data === "number")) {
				return `${(100 * data).toFixed(1)}%`;
			}
			if (typeof data === "string") {
				return data;
			} else {
				return data.toFixed(3).replace(/\.?0+$/, "");
			}
		});
}

/**
 * Create a table.
 * @param rows {TableRow[]}
 * @returns {Node[]}
 */
function table(rows) {
	/** @type {[number, TableRow[]][]} */
	const pages = d3.groups(rows, (_, index) => Math.floor(index / 25))
				.map(([key, group]) => [key + 1, group]);

	const heading = d3.create("h4").text("Details of the perturbed complexes.");
	const table = d3.create("table").attr("id", "data-table");
  const thead = table.append("thead");

  thead
    .append("tr")
    .selectAll("th")
    .data([
      "Complex ID",
      "Complex Name",
      "Coverage",
      "Perturbation Type",
      "Perturbation Score",
      "Normalized Score",
      "Subunit ID",
      "Gene Name",
      "log2(FC)",
      "Adj. p-value",
    ])
    .enter()
    .append("th")
    .attr("tabindex", (_, i) => i)
    .attr("role", "button")
    .attr("aria-label", "Sort column")
    .text(identity);

	const tbody = table.append("tbody");
	renderTableData(tbody, pages[0]?.[1] ?? []);

	const paginationContainer = d3.create("div")
				.attr("class", "table-pagination-container");
	paginationContainer
		.selectAll()
		.data(pages)
		.join("span")
		.style("border", (d, i) => i === 0 ? "1px solid grey" : "none")
	  .on("click", (event, [index, page]) => {
			tbody.node()?.replaceChildren();
			renderTableData(tbody, page);
			d3
				.selectAll("div.table-pagination-container > span")
				.style("border", "none");
			d3
				.select(`div.table-pagination-container > span:nth-child(${index})`)
				.style("border", "1px solid grey");
		})
		.text(([page, _]) => page.toString());

  const result = [heading.node(), table.node(), paginationContainer.node()];
	if (result.every((x) => x !== null)) {
		return result;
	}
  return [];
}

export { barPlot, histogramPlot, vennPlot, volcanoPlot, table };
