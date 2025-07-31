/**
 * @param {Array<number>} data - the data to summarize
 * @return {[number, number]}  - [minimum value, maximum value]
 */
function minmax(data) {
	return data.reduce(([min, max], value) => [Math.min(min, value), Math.max(max, value)], [Infinity, -Infinity]);
}

/**
 * Find root of the circle-circle intersection equation:
 *   from: https://mathworld.wolfram.com/Circle-CircleIntersection.html
 * @param {number} R         - the radius of the first circle
 * @param {number} r         - the radius of the second circle
 * @param {number} AB        - the area of the intersection
 * @param {[number, number]} - an initial range for the root
 * @return {number}          - the x coordinate of the second circle
 */
function bisect(R, r, AB, bracket) {
	let left = bracket[0];
	let right = bracket[1];
	const r2 = r * r;
	const R2 = R * R;
	const area = (d) => {
		const n = Math.max(0, Math.min(1, (d*d + r2 - R2) / (2 * d * r)));
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


export { minmax, bisect };
