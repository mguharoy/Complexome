// @ts-check
/**
 * @typedef Field A tabular data entry.
 *
 * @property offset {number} The index into the buffer that this entry starts at.
 * @property length {number} The length of this field.
 */

/**
 * @typedef TSV Tabular data
 *
 * @property buffer {string} The original text.
 * @property headers {string[]} Header text values.
 * @property columns {Field[][]} Tabular data fields arranged in column-major order.
 * @property rows {number} The total number of _data_ rows in the table.
 */

/**
 * Parse CSV formatted data.
 * @param {string} text
 * @param {string} separator - Field separator
 * @returns {TSV}
 */
export function tsvParse(text, separator) {
  let finishedHeaders = false;
  let fieldIndex = 0;
  let columns = [];
  let columnIndex = 0;
  /** @type {TSV} */
  let tsv = {
    buffer: text,
    headers: [],
    columns: [],
    rows: 0,
  };
  let rows = 0;

  for (let i = 0; i < text.length; i++) {
    if (text[i] === separator) {
      if (!finishedHeaders) {
        tsv.headers.push(text.slice(fieldIndex, i));
      } else if (columns.length <= columnIndex) {
        columns.push([{ offset: fieldIndex, length: i - fieldIndex }]);
      } else {
        columns[columnIndex]?.push({
          offset: fieldIndex,
          length: i - fieldIndex,
        });
      }
      columnIndex++;
      fieldIndex = i + 1;
    } else if (text[i] === "\n") {
      if (finishedHeaders) {
        rows++;
        if (columns.length <= columnIndex) {
          columns.push([{ offset: fieldIndex, length: i - fieldIndex }]);
        } else {
          columns[columnIndex]?.push({
            offset: fieldIndex,
            length: i - fieldIndex,
          });
        }
      } else {
        tsv.headers.push(text.slice(fieldIndex, i));
      }
      finishedHeaders = true;
      columnIndex = 0;
      fieldIndex = i + 1;
    }
  }

  tsv.columns = columns;
  tsv.rows = rows;
  return tsv;
}

/**
 * Extract rows from tabular data.
 * @param {TSV} tsv - A table of data.
 * @returns {Generator<string[]>}
 */
export function* rows(tsv) {
  let row = 0;
  while (row < tsv.rows) {
    yield tsv.columns.map((column) => {
      const offset = column[row]?.offset;
      const length = column[row]?.length;
      if (offset !== undefined && length !== undefined) {
        return tsv.buffer.slice(offset, offset + length);
      } else {
        return "";
      }
    });

    row++;
  }
}

/**
 * Extract rows from tabular data matching a column filter.
 * @param {TSV} tsv - A table of data.
 * @param {string[]} columnFilter - A list of column headers to extract rows from.
 * @returns {Generator<Record<string, string>>}
 */
export function* columnFilteredRows(tsv, columnFilter) {
  let row = 0;
  while (row < tsv.rows) {
    yield columnFilter.reduce(
      (/** @type {Record<string, string>} */ obj, field) => {
        const column = tsv.headers.indexOf(field);
        if (column > -1) {
          const offset = tsv.columns[column]?.[row]?.offset;
          const length = tsv.columns[column]?.[row]?.length;
          if (offset !== undefined && length !== undefined) {
            obj[field] = tsv.buffer.slice(offset, offset + length);
          }
        }
        return obj;
      },
      {},
    );

    row++;
  }
}
