export function tsvParse(text, separator) {
	let finishedHeaders = false;
  let fieldIndex = 0;
  let columns = [];
  let columnIndex = 0;
  let tsv = {
    buffer: text,
    headers: [],
    columns: [],
		rows: 0
  };
	let rows = 0;

  for (let i = 0; i < text.length; i++) {
    if (text[i] === separator) {
      if (!finishedHeaders) {
        tsv.headers.push(text.slice(fieldIndex, i));
      }
      else if (columns.length <= columnIndex) {
        columns.push([{offset: fieldIndex, length: i - fieldIndex}]);
      } else {
        columns[columnIndex]?.push({offset: fieldIndex, length: i - fieldIndex});
      }
      columnIndex++;
      fieldIndex = i + 1;
    } else if (text[i] === '\n') {
      if (finishedHeaders){
				rows++;
        if (columns.length <= columnIndex) {
          columns.push([{offset: fieldIndex, length: i - fieldIndex}]);
        } else {
          columns[columnIndex]?.push({offset: fieldIndex, length: i - fieldIndex});
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

export function* rows(tsv) {
	let row = 0;
	while (row < tsv.rows) {
		yield tsv.columns.map((column) => {
			const offset = column[row]?.offset;
			const length = column[row]?.length;
			if (offset !== undefined && length !== undefined) {
				return tsv.buffer.slice(offset, offset + length);
			}
		});

		row++;
	}
}

export function* columnFilteredRows(tsv, columnFilter) {
	let row = 0;
	while (row < tsv.rows) {
		yield columnFilter.reduce((obj, field) => {
			const column = tsv.headers.indexOf(field);
			if (column > -1){
				const offset = tsv.columns[column][row]?.offset;
				const length = tsv.columns[column][row]?.length;
				if (offset !== undefined && length !== undefined) {
					obj[field] = tsv.buffer.slice(offset, offset + length);
				}
			}
			return obj;
		}, {});

		row++;
	}
}
