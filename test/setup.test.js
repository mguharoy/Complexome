import test from 'node:test';
import assert from 'node:assert';

import { bisect, minmax } from '../src/numeric.mjs';

test('The numeric module', (t) => {
	t.test('bisect 1', (t) => {
		assert.strictEqual(Math.abs(bisect(0.8, 0.4, 0.08359, [0.8, 1.2]) - 1) < 0.01, true);
	});

	t.test('bisect 2', (t) => {
		assert.strictEqual(Math.abs(bisect(0.5, 0.5, 0.39634, [0.001, 1]) - 0.4) < 0.01, true);
	});

	t.test('bisect 3', (t) => {
		assert.strictEqual(Math.abs(bisect(0.25, 0.75, 0.00896, [0.25, 1]) - 0.95) < 0.01, true)
	});

	t.test('minmax 1', (t) => {
		assert.deepStrictEqual(minmax([1, 2]), [1, 2]);
	});

	t.test('minmax 2', (t) => {
		assert.deepStrictEqual(minmax([2, 3, 1]), [1, 3]);
	});

	t.test('minmax 3', (t) => {
		assert.deepStrictEqual(minmax([2]), [2, 2]);
	});
})
