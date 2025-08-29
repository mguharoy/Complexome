import { columnFilteredRows, tsvParse, rows } from "./tsv.js";

self.cache = new Map();

function expandParticipants(participants) {
	let members = new Set();
	for (const participant of participants.split("|")) {
		const participant_id = participant.split("(")[0];
		if (participant_id?.includes("[")) {
			for (const paralog of participant_id.replaceAll("[", "").replaceAll("]", "").split(",")) {
				members.add(paralog);
			}
		} else if (!members.has(participant_id)) {
			members.add(participant_id);
		}
	}

	return members;
}

self.processComplexomeData = function(tsv) {
	const complexes = new Map();
	const complexNames = new Map();
	const ontologyTerms = new Map();
	const crossRefs = new Map();

	for (const row of columnFilteredRows(tsv, ["#Complex ac", "Recommended name", "Taxonomy identifier", "Identifiers (and stoichiometry) of molecules in complex", "Expanded participant list", "Go Annotations", "Cross references"])) {
		const complexID = row["#Complex ac"];
		const participants = row["Identifiers (and stoichiometry) of molecules in complex"];
		const complexMembers = expandParticipants(participants).union(participants.includes("CPX-") ? expandParticipants(row["Expanded participant list"]) : new Set());
		if (!complexes.has(complexID)) {
			complexes.set(complexID, complexMembers);
			complexNames.set(complexID, row["Recommended name"]);
		} else {
			throw new Error(`Complex ID already exists: ${complexID}`);
		}

		const geneOntologyTerms = row["Go Annotations"].split("|");
		if (!ontologyTerms.has(complexID)) {
			ontologyTerms.set(complexID, geneOntologyTerms);
		} else {
			throw new Error(`Complex ID already exists: ${complexID}`);
		}

		crossRefs.set(complexID, row["Cross references"]
			.split("|")
			.filter((ref) => ref.startsWith("reactome:"))
			.map((ref) => ref.slice(9)));
	}

	if ((complexes.size === complexNames.size) && (complexNames.size === ontologyTerms.size)) {
		return [complexes, complexNames, ontologyTerms, crossRefs];
	}

	throw new Error(`I found an inconsistency: num_complexes=${complexes.size}, num_names=${complexNames.size}, num_gene_ontology_terms=${ontologyTerms.size}`);
};

function processUserData(tsv) {
	const result = new Map();

	for (const row of rows(tsv)) {
		const log2fc = parseFloat(row[1]);
		const adjpval = parseFloat(row[2]);
		if (isNaN(log2fc) || isNaN(adjpval)) {
			continue;
		}
		result.set(row[0], [log2fc, adjpval]);
	}

	return result;
}

self.onmessage = async (event) => {
	const { op, data } = event.data;

	switch (op) {

		case "getcomplex": {
			const complexid = data;
			console.log(`Worker now requesting ${complexid}`);
			if (self.cache.has(complexid)) {
				self.postMessage(self.cache.get(complexid));
			} else {
				const response = await fetch(`https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/${complexid}.tsv`);
				if (response.ok) {
					const text = await response.text();
					const tsv = tsvParse(text, "\t");
					const cpxdata = self.processComplexomeData(tsv);
					self.cache.set(complexid, cpxdata);
					self.postMessage({ "id": complexid, "complex": cpxdata });
				} else {
					self.postMessage({ "error": "network error" });
				}
			}
			break;
		}

		case "csv": {
			const result = tsvParse(data, ",");
			self.postMessage({ "userdata": processUserData(result) });
			break;
		}
	}
};
