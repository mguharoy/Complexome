// @ts-check
/**
 * @template T
 * @param {T} x - value to return
 * @returns {T}
 */
function identity(x) {
	return x;
}

/**
 * Get the min and max from an array of numbers.
 * @overload
 * @param {number[]} data - The data to summarize.
 * @returns {[number, number]} - [minimum value, maximum value]
 */
/**
 * Get the min and max from an array using an accessor function.
 * @template T
 * @overload
 * @param {T[]} data - the data to summarize
 * @param {(x: T) => number} [accessor] - Access members to compute min and max on.
 * @return {[number, number]}  - [minimum value, maximum value]
 */
/**
 * @param {number[] | T[]} data
 * @param {((x: T) => number) | undefined} [accessor]
 */
function minmax(data, accessor) {
	if (accessor === undefined) {
		const numberData = /** @type {number[]} */ (data);
		return numberData.reduce(([min, max], datum) => [Math.min(min ?? Infinity, datum), Math.max(max ?? -Infinity, datum)], [Infinity, -Infinity]);
	} else {
		const genericData = /** @type {T[]} */ (data);
		return genericData.reduce(
			([min, max], datum) => {
				const value = accessor(datum);
				return [Math.min(min ?? Infinity, value), Math.max(max ?? -Infinity, value)];
			},
			[Number.POSITIVE_INFINITY, Number.NEGATIVE_INFINITY]
		);
	}
}

/**
 * Find root of the circle-circle intersection equation:
 *   from: https://mathworld.wolfram.com/Circle-CircleIntersection.html
 * @param {number} R         - the radius of the first circle
 * @param {number} r         - the radius of the second circle
 * @param {number} AB        - the area of the intersection
 * @param {[number, number]} bracket - an initial range for the root
 * @return {number}          - the x coordinate of the second circle
 */
function bisect(R, r, AB, bracket) {
	let left = bracket[0];
	let right = bracket[1];
	const r2 = r * r;
	const R2 = R * R;
	const area = (/** @type {number} */ d) => {
		const n = Math.max(0, Math.min(1, (d * d + r2 - R2) / (2 * d * r)));
		const o = Math.max(0, Math.min(1, (d * d + R2 - r2) / (2 * d * R)));
		const p = Math.max(0.001, (-d + r + R) * (d + r - R) * (d - r + R) * (d + r + R));
		return r2 * Math.acos(n) + R2 * Math.acos(o) - 0.5 * Math.sqrt(p);
	}
	for (let n = 0; n < 25; n++) {
		const d = (left + right) / 2;
		const f = area(d);
		if ((Math.abs(AB - f) < 0.0001) || ((right - left) / 2 < 0.001)) {
			return d;
		} else if (Math.sign(f - AB) === Math.sign(area(left) - AB)) {
			left = d;
		} else {
			right = d;
		}
	}
	throw new Error("Could not bisect function");
}


export { identity, minmax, bisect };
