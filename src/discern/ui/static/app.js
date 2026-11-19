const slider = document.getElementById("prevalence");
if (slider) {
  const output = document.getElementById("prevalence-value");
  slider.addEventListener("input", () => {
    output.value = (Number(slider.value) * 100).toFixed(1) + "%";
  });
}

// Live CSV Quality Evaluation
const REQUIRED_PROFILE_COLS = ["sample_id", "gene", "mutation", "frs"];
const ALLOWED_GENES = new Set(["Rv0678", "atpE", "pepQ", "mmpL5", "mmpS5"]);
const SAMPLE_ID_RE = /^[A-Za-z0-9_\-\.]{1,64}$/;

function parseCSV(text) {
  const lines = [];
  let row = [];
  let inQuotes = false;
  let current = "";
  for (let i = 0; i < text.length; i++) {
    const char = text[i];
    const next = text[i + 1];
    if (char === '"') {
      if (inQuotes && next === '"') {
        current += '"';
        i++;
      } else {
        inQuotes = !inQuotes;
      }
    } else if (char === ',' && !inQuotes) {
      row.push(current.trim());
      current = "";
    } else if ((char === '\r' || char === '\n') && !inQuotes) {
      if (char === '\r' && next === '\n') i++;
      row.push(current.trim());
      current = "";
      if (row.some(c => c.length > 0)) lines.push(row);
      row = [];
    } else {
      current += char;
    }
  }
  if (current.length > 0 || row.length > 0) {
    row.push(current.trim());
    if (row.some(c => c.length > 0)) lines.push(row);
  }
  return lines;
}

function isValidMutationSyntax(mut) {
  if (!mut) return false;
  const suf = mut.split("@").pop();
  if (!suf) return false;
  // Synonymous: same first & last letter, length >= 2, no _, no !
  if (suf.length > 1 && suf[0] === suf[suf.length - 1] && !suf.includes("_") && !suf.endsWith("!")) return true;
  // Null call: ends with Z
  if (suf.endsWith("Z")) return true;
  // Stop codon: ends with !
  if (suf.endsWith("!")) return true;
  // Indels
  if (suf.includes("_ins_") || suf.includes("_del_") || suf.includes("indel") || suf.startsWith("del_") || suf.startsWith("-") || suf.includes("ins") || suf.includes("del")) return true;
  // Promoter: e.g. c-11a, t-8c
  if (/^[acgt]-\d+[acgt]?$/i.test(suf) || /^[acgt]-\d+$/i.test(suf)) return true;
  // Missense: e.g. R50Q, G41V
  if (/^[A-Z]\d+[A-Z]$/.test(suf)) return true;
  return false;
}

function evaluateCSVQuality(text) {
  const stripped = text.trim();
  if (!stripped) return null;

  let table;
  try {
    table = parseCSV(stripped);
  } catch (err) {
    return {
      completeness: 0,
      totalRows: 0,
      validRows: 0,
      missingCount: 0,
      duplicateCount: 0,
      formatCount: 1,
      issues: [{ line: 1, type: "format", field: "csv", message: "Malformed CSV structure: " + err.message }],
      isValid: false
    };
  }

  if (!table || table.length === 0) return null;

  const rawHeader = table[0];
  const headers = rawHeader[rawHeader.length - 1] === "banner" ? rawHeader.slice(0, -1) : rawHeader;

  const isRoster = headers.length === 1 && headers[0] === "sample_id";
  const requiredCols = isRoster ? ["sample_id"] : REQUIRED_PROFILE_COLS;

  const issues = [];
  const rowsWithIssues = new Set();

  // Check missing required columns in header
  for (const req of requiredCols) {
    if (!headers.includes(req)) {
      issues.push({ line: 1, type: "format", field: req, message: `Header line 1 is missing required column '${req}'` });
      rowsWithIssues.add(1);
    }
  }

  const dataRows = table.slice(1);
  const totalRows = dataRows.length;
  if (totalRows === 0) {
    issues.push({ line: 1, type: "format", field: "data", message: "No data rows found in the CSV" });
    return {
      completeness: 0,
      totalRows: 0,
      validRows: 0,
      missingCount: 0,
      duplicateCount: 0,
      formatCount: issues.length,
      issues,
      isValid: false
    };
  }

  const totalRequiredCells = totalRows * requiredCols.length;
  let filledRequiredCells = 0;

  const seenExact = new Map();
  const seenVariants = new Map();

  dataRows.forEach((rowVals, idx) => {
    const lineNum = idx + 2; // 1-based, header is line 1

    if (rowVals.length < headers.length) {
      issues.push({ line: lineNum, type: "format", field: "columns", message: `Line ${lineNum}: row has ${rowVals.length} columns, expected ${headers.length}` });
      rowsWithIssues.add(lineNum);
    } else if (rowVals.length > headers.length) {
      issues.push({ line: lineNum, type: "format", field: "columns", message: `Line ${lineNum}: row has extra ${rowVals.length - headers.length} columns` });
      rowsWithIssues.add(lineNum);
    }

    const rowObj = {};
    headers.forEach((h, colIdx) => {
      rowObj[h] = rowVals[colIdx] !== undefined ? rowVals[colIdx].trim() : "";
    });

    // Check completeness of required columns
    requiredCols.forEach(req => {
      const val = rowObj[req] || "";
      if (val !== "") {
        filledRequiredCells++;
      } else {
        issues.push({ line: lineNum, type: "missing", field: req, message: `Line ${lineNum}: missing required value in '${req}'` });
        rowsWithIssues.add(lineNum);
      }
    });

    // Check duplicate rows
    const exactKey = headers.map(h => rowObj[h] || "").join("||");
    if (seenExact.has(exactKey)) {
      const prior = seenExact.get(exactKey);
      issues.push({ line: lineNum, type: "duplicate", field: "row", message: `Line ${lineNum}: exact duplicate of line ${prior}` });
      rowsWithIssues.add(lineNum);
    } else {
      seenExact.set(exactKey, lineNum);
    }

    // Check duplicate variants for same sample
    if (rowObj.sample_id && rowObj.gene && rowObj.mutation) {
      const varKey = `${rowObj.sample_id}@@${rowObj.gene}@@${rowObj.mutation}`;
      if (seenVariants.has(varKey) && !seenExact.has(exactKey)) {
        const prior = seenVariants.get(varKey);
        issues.push({ line: lineNum, type: "duplicate", field: "variant", message: `Line ${lineNum}: duplicate variant ${rowObj.gene}@${rowObj.mutation} for sample '${rowObj.sample_id}' (earlier at line ${prior})` });
        rowsWithIssues.add(lineNum);
      } else if (!seenVariants.has(varKey)) {
        seenVariants.set(varKey, lineNum);
      }
    }

    // Check format of sample_id
    if (rowObj.sample_id) {
      if (!SAMPLE_ID_RE.test(rowObj.sample_id)) {
        issues.push({ line: lineNum, type: "format", field: "sample_id", message: `Line ${lineNum}: invalid sample_id formatting` });
        rowsWithIssues.add(lineNum);
      }
    }

    // Check format of gene
    if (rowObj.gene) {
      if (!ALLOWED_GENES.has(rowObj.gene)) {
        issues.push({ line: lineNum, type: "format", field: "gene", message: `Line ${lineNum}: invalid gene '${rowObj.gene}' (expected: Rv0678, atpE, pepQ, mmpL5, mmpS5)` });
        rowsWithIssues.add(lineNum);
      }
    }

    // Check format of mutation
    if (rowObj.mutation) {
      if (!isValidMutationSyntax(rowObj.mutation)) {
        issues.push({ line: lineNum, type: "format", field: "mutation", message: `Line ${lineNum}: invalid mutation syntax '${rowObj.mutation}'` });
        rowsWithIssues.add(lineNum);
      }
    }

    // Check format of frs
    if (rowObj.frs !== undefined && rowObj.frs !== "") {
      const num = Number(rowObj.frs);
      if (isNaN(num) || num < 0 || num > 1) {
        issues.push({ line: lineNum, type: "format", field: "frs", message: `Line ${lineNum}: frs must be a valid number between 0.0 and 1.0` });
        rowsWithIssues.add(lineNum);
      }
    }
  });

  const completeness = totalRequiredCells > 0 ? (filledRequiredCells / totalRequiredCells * 100) : 0;
  const missingCount = issues.filter(i => i.type === "missing").length;
  const duplicateCount = issues.filter(i => i.type === "duplicate").length;
  const formatCount = issues.filter(i => i.type === "format").length;
  const validRows = Math.max(0, totalRows - rowsWithIssues.size);

  return {
    completeness: Math.round(completeness * 10) / 10,
    totalRows,
    validRows,
    missingCount,
    duplicateCount,
    formatCount,
    issues,
    isValid: issues.length === 0
  };
}

function renderUploadQuality(text) {
  const card = document.getElementById("upload-quality-card");
  if (!card) return;

  const result = evaluateCSVQuality(text);
  if (!result) {
    card.hidden = true;
    return;
  }

  card.hidden = false;
  const scoreVal = document.getElementById("quality-score-val");
  const scoreBadge = document.getElementById("quality-score-badge");
  const totalRowsEl = document.getElementById("quality-total-rows");
  const validRowsEl = document.getElementById("quality-valid-rows");
  const missingCountEl = document.getElementById("quality-missing-count");
  const duplicateCountEl = document.getElementById("quality-duplicate-count");
  const formatCountEl = document.getElementById("quality-formatting-count");
  const issuesContainer = document.getElementById("quality-issues-container");
  const issuesBody = document.getElementById("quality-issues-body");
  const issuesSummary = document.getElementById("quality-issues-summary");
  const cleanMsg = document.getElementById("quality-clean-msg");

  if (scoreVal) scoreVal.textContent = result.completeness.toFixed(1) + "%";
  if (totalRowsEl) totalRowsEl.textContent = result.totalRows;
  if (validRowsEl) validRowsEl.textContent = `${result.validRows} / ${result.totalRows}`;
  if (missingCountEl) missingCountEl.textContent = result.missingCount;
  if (duplicateCountEl) duplicateCountEl.textContent = result.duplicateCount;
  if (formatCountEl) formatCountEl.textContent = result.formatCount;

  if (result.isValid && result.completeness === 100) {
    if (scoreBadge) {
      scoreBadge.className = "quality-score-badge quality-pass";
    }
    if (cleanMsg) cleanMsg.hidden = false;
    if (issuesContainer) issuesContainer.hidden = true;
  } else {
    if (scoreBadge) {
      scoreBadge.className = "quality-score-badge quality-warn";
    }
    if (cleanMsg) cleanMsg.hidden = true;
    if (issuesContainer) {
      issuesContainer.hidden = false;
      if (issuesSummary) {
        issuesSummary.textContent = `${result.issues.length} issue(s) detected across ${result.totalRows - result.validRows} of ${result.totalRows} row(s):`;
      }
      if (issuesBody) {
        issuesBody.innerHTML = "";
        result.issues.slice(0, 30).forEach(iss => {
          const tr = document.createElement("tr");
          tr.innerHTML = `<td>Line ${iss.line}</td><td><span class="issue-pill issue-${iss.type}">${iss.type}</span></td><td><code>${iss.field}</code></td><td>${iss.message}</td>`;
          issuesBody.appendChild(tr);
        });
      }
    }
  }
}

const fileInput = document.getElementById("csv-file");
if (fileInput) {
  const status = document.getElementById("file-status");
  const contents = document.getElementById("csv");
  const form = fileInput.form;
  const submit = form.querySelector('button[type="submit"]');
  const action = submit ? submit.textContent : "Run model scoring";
  const zone = fileInput.closest("[data-dropzone]");
  const maxBytes = Number(fileInput.dataset.maxBytes);
  let selection = 0;

  function setStatus(message, error = false) {
    status.textContent = message;
    status.dataset.state = error ? "error" : "ready";
    if (error) fileInput.setAttribute("aria-invalid", "true");
    else fileInput.removeAttribute("aria-invalid");
  }

  async function loadFile(file) {
    if (!file) return;
    const current = ++selection;
    contents.value = "";
    const serverError = form.querySelector("[data-upload-error]");
    if (serverError) serverError.hidden = true;
    fileInput.setAttribute("aria-describedby", zone ? "upload-help file-status" : "file-status");
    if (submit) submit.disabled = false;
    if (!/\.csv$/i.test(file.name) || !file.size || file.size > maxBytes) {
      setStatus("Choose a non-empty .csv file of up to " + maxBytes / 1024 + " KiB.", true);
      renderUploadQuality("");
      return;
    }
    if (submit) submit.disabled = true;
    setStatus("Reading file…");
    try {
      const text = await file.text();
      if (current !== selection) return;
      if (!text.trim()) {
        setStatus("This CSV is empty. Choose a file that contains a recorded sample.", true);
        renderUploadQuality("");
        return;
      }
      contents.value = text;
      setStatus(file.name + " loaded. Ready for scoring.");
      renderUploadQuality(text);
    } catch {
      if (current === selection) setStatus("The file could not be read. You can paste the CSV contents.", true);
    } finally {
      if (current === selection && submit) submit.disabled = false;
    }
  }

  fileInput.addEventListener("change", () => loadFile(fileInput.files[0]));

  if (contents) {
    contents.addEventListener("input", () => {
      renderUploadQuality(contents.value);
    });
  }

  if (zone && "DataTransfer" in window) {
    zone.dataset.enabled = "";
    zone.addEventListener("dragover", event => {
      if (!event.dataTransfer.types.includes("Files")) return;
      event.preventDefault();
      event.dataTransfer.dropEffect = "copy";
      zone.classList.add("dragging");
    });
    zone.addEventListener("dragleave", event => {
      if (!zone.contains(event.relatedTarget)) zone.classList.remove("dragging");
    });
    zone.addEventListener("drop", event => {
      event.preventDefault();
      zone.classList.remove("dragging");
      if (event.dataTransfer.files.length !== 1) {
        ++selection;
        contents.value = "";
        fileInput.value = "";
        if (submit) submit.disabled = false;
        setStatus("Drop one CSV file at a time.", true);
        renderUploadQuality("");
        return;
      }
      fileInput.value = "";
      loadFile(event.dataTransfer.files[0]);
    });
  }

  form.addEventListener("submit", event => {
    if (!contents.value.trim()) {
      event.preventDefault();
      setStatus("Choose a CSV file or paste its contents to continue.", true);
      fileInput.focus();
      return;
    }
    if (submit) {
      submit.textContent = "Scoring…";
      submit.disabled = true;
    }
  });

  window.addEventListener("pageshow", () => {
    if (submit) {
      submit.textContent = action;
      submit.disabled = false;
    }
  });

  const uploadError = form.querySelector("[data-upload-error]");
  if (uploadError) uploadError.focus();
}
