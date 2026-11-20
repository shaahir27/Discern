const slider = document.getElementById("prevalence");
if (slider) {
  const output = document.getElementById("prevalence-value");
  slider.addEventListener("input", () => {
    output.value = (Number(slider.value) * 100).toFixed(1) + "%";
  });
}

const studyNavigation = document.querySelector(".study-navigation");
if (studyNavigation && !studyNavigation.querySelector("[aria-current]") &&
    matchMedia("(max-width: 950px)").matches) {
  studyNavigation.open = false;
}

const queueFilter = document.getElementById("queue-filter");
if (queueFilter) {
  const rows = [...document.querySelectorAll("[data-sample-id]")];
  const count = document.querySelector("[data-queue-count]");
  const empty = document.querySelector("[data-queue-empty]");
  queueFilter.closest("[data-queue-controls]").hidden = !rows.length;
  queueFilter.addEventListener("input", () => {
    const query = queueFilter.value.trim().toLowerCase();
    let shown = 0;
    for (const row of rows) {
      row.hidden = !row.dataset.sampleId.toLowerCase().includes(query);
      if (!row.hidden) shown++;
    }
    count.textContent = shown === rows.length ? `${shown} samples` : `${shown} of ${rows.length} samples`;
    empty.hidden = shown > 0;
  });
}

const panelFilter = document.getElementById("panel-filter");
if (panelFilter) {
  const assays = Array.from(document.querySelectorAll("[data-assay]"));
  panelFilter.closest("[data-panel-controls]").hidden = false;
  panelFilter.addEventListener("change", () => {
    for (const row of assays) {
      row.hidden = Boolean(panelFilter.value) && row.dataset.assay !== panelFilter.value;
    }
  });
}

const cohortComparison = document.querySelector("[data-cohort-comparison]");
if (cohortComparison) {
  const options = cohortComparison.querySelector(".comparison-options");
  const rows = cohortComparison.querySelectorAll("[data-compare-metric]");
  options.hidden = false;
  options.addEventListener("change", event => {
    for (const row of rows) {
      const value = Number(event.target.dataset[row.dataset.compareMetric]);
      row.querySelector("meter").value = value;
      const output = row.querySelector("data");
      output.value = value;
      output.textContent = `${(Math.trunc(value * 1000) / 10).toFixed(1)}%`;
      if (!matchMedia("(prefers-reduced-motion: reduce)").matches) {
        output.getAnimations().forEach(animation => animation.cancel());
        output.animate([{opacity: .5}, {opacity: 1}], {duration: 160, easing: "ease-out"});
      }
    }
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
    const isExactDup = seenExact.has(exactKey);
    if (isExactDup) {
      const prior = seenExact.get(exactKey);
      issues.push({ line: lineNum, type: "duplicate", field: "row", message: `Line ${lineNum}: exact duplicate of line ${prior}` });
      rowsWithIssues.add(lineNum);
    } else {
      seenExact.set(exactKey, lineNum);
    }

    // Check duplicate variants for same sample
    if (rowObj.sample_id && rowObj.gene && rowObj.mutation) {
      const varKey = `${rowObj.sample_id}@@${rowObj.gene}@@${rowObj.mutation}`;
      if (seenVariants.has(varKey) && !isExactDup) {
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
        result.issues.slice(0, 30).forEach(issue => {
          const row = issuesBody.insertRow();
          for (const content of [`Line ${issue.line}`, issue.type, issue.field, issue.message]) {
            row.insertCell().textContent = content;
          }
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
  const assay = form.elements.assay;
  const formatStatus = form.querySelector("[data-format-status]");
  const filename = document.createElement("input");
  filename.type = "hidden";
  filename.name = "filename";
  form.append(filename);
  let selection = 0;
  let manualAssay = false;
  assay?.addEventListener("change", () => { manualAssay = true; });

  function inspectUpload(text) {
    const header = text.trimStart().split(/\r?\n/, 1)[0];
    const columns = header.toLowerCase().split(/[,\t]/).map(value => value.trim().replace(/^"|"$/g, ""));
    const isDiscern = columns[0] === "sample_id" &&
      (columns[1] === "gene" && columns[2] === "mutation" && columns[3] === "frs" ||
       columns.length === 1 || columns[1] === "banner");
    let label = text.trim() ? "Discern CSV" : "";
    let inferred = assay?.options[0]?.value ?? "wgs";
    if (header.startsWith("{")) {
      label = "TB-Profiler JSON";
    } else if (!isDiscern && (/freq.*%|variant frequency|gene name/i.test(header))) {
      label = "Deeplex CSV";
      inferred = "deeplex_myc_tb";
    } else if (!isDiscern && columns.some(value => ["vaf", "af", "freq"].includes(value))) {
      label = "AmPORE-TB CSV/TSV";
      inferred = "ampore_tb";
    } else if (isDiscern && columns.includes("assay")) {
      const assays = new Set(parseCSV(text).slice(1)
        .map(row => row[columns.indexOf("assay")]?.trim()).filter(Boolean));
      if (assays.size === 1 && assay) {
        const candidate = assays.values().next().value;
        if (Array.from(assay.options).some(option => option.value === candidate)) inferred = candidate;
      }
    }
    if (assay && !manualAssay) assay.value = inferred;
    if (formatStatus) {
      formatStatus.textContent = label ? `Detected format: ${label}` : "";
      formatStatus.hidden = !label;
    }
    renderUploadQuality(isDiscern ? text : "");
  }

  function setStatus(message, error = false) {
    status.textContent = message;
    status.dataset.state = error ? "error" : "ready";
    if (zone) zone.dataset.state = error ? "error" : contents.value ? "ready" : "idle";
    if (error) fileInput.setAttribute("aria-invalid", "true");
    else fileInput.removeAttribute("aria-invalid");
  }

  async function loadFile(file) {
    if (!file) return;
    const current = ++selection;
    contents.value = "";
    filename.value = "";
    const serverError = form.querySelector("[data-upload-error]");
    if (serverError) serverError.hidden = true;
    fileInput.setAttribute("aria-describedby", zone ? "upload-help file-status" : "file-status");
    if (submit) submit.disabled = false;
    inspectUpload("");
    if (!/\.(csv|tsv|json)$/i.test(file.name) || !file.size || file.size > maxBytes) {
      setStatus("Choose a non-empty CSV, TSV, or JSON file of up to " + maxBytes / 1024 + " KiB.", true);
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
      filename.value = file.name;
      setStatus(file.name + " is ready to check.");
      inspectUpload(text);
    } catch {
      if (current === selection) setStatus("The file could not be read. You can paste the CSV contents.", true);
    } finally {
      if (current === selection && submit) submit.disabled = false;
    }
  }

  fileInput.addEventListener("change", () => loadFile(fileInput.files[0]));

  if (contents) {
    contents.addEventListener("input", () => {
      ++selection;
      filename.value = "";
      if (submit) submit.disabled = false;
      inspectUpload(contents.value);
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
        filename.value = "";
        fileInput.value = "";
        if (submit) submit.disabled = false;
        setStatus("Drop one file at a time.", true);
        inspectUpload("");
        return;
      }
      fileInput.value = "";
      loadFile(event.dataTransfer.files[0]);
    });
  }

  form.addEventListener("submit", event => {
    if (!contents.value.trim()) {
      event.preventDefault();
      setStatus("Choose a file or paste its contents to continue.", true);
      fileInput.focus();
      return;
    }
    if (submit) {
      submit.textContent = "Checking…";
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

// ---------------------------------------------------------------------------
// Longitudinal Serial-Isolate Tracking & Clonal-Ribbon Timeline
// ---------------------------------------------------------------------------
function initLongitudinalTimeline() {
  const dataScript = document.getElementById("timeline-data");
  if (!dataScript) return;

  let timeline;
  try {
    timeline = JSON.parse(dataScript.textContent);
  } catch (e) {
    return;
  }

  const visits = timeline.timepoints;
  const n = visits.length;
  if (n < 2) return;

  const variants = timeline.variants;
  const trajectories = timeline.trajectories;

  // Series colors distinguish variants independently of the catalogue call.
  const colorPalette = ["#126580", "#957515", "#567d76", "#6a516d"];

  const variantColors = {};
  variants.forEach((v, idx) => {
    variantColors[v] = colorPalette[idx % colorPalette.length];
  });

  // Color legend boxes
  document.querySelectorAll(".legend-color-box").forEach(el => {
    const v = el.getAttribute("data-variant");
    if (variantColors[v]) el.style.backgroundColor = variantColors[v];
  });

  const svg = document.getElementById("muller-svg");
  const ribbonsGroup = document.getElementById("ribbons-group");
  const markersGroup = document.getElementById("visit-markers-group");
  const cursor = document.getElementById("scrubber-cursor");

  const scrubber = document.getElementById("visit-scrubber");
  const activeLabel = document.getElementById("active-visit-label");
  const btnPrev = document.getElementById("btn-prev-visit");
  const btnNext = document.getElementById("btn-next-visit");
  const btnPlay = document.getElementById("btn-play-pause");
  const visitButtons = Array.from(document.querySelectorAll("[data-visit-index]"));

  const leftX = 60;
  const rightX = 720;
  const topY = 20;
  const botY = 240;
  const usableW = rightX - leftX;
  const usableH = botY - topY;

  const firstDay = visits[0].sort_key;
  const spanDays = visits[n - 1].sort_key - firstDay;
  const getX = idx => leftX + (spanDays ? (visits[idx].sort_key - firstDay) / spanDays : idx / (n - 1)) * usableW;
  const getY = frs => botY - Math.min(1.0, Math.max(0.0, frs)) * usableH;

  // Calculate stacked ribbon intervals for each visit
  // order: variants then wildtype
  const allLayers = [...timeline.core_variants, "__WT__"];
  const stackedTops = {};
  const stackedBots = {};

  allLayers.forEach(layer => {
    stackedTops[layer] = new Array(n).fill(0);
    stackedBots[layer] = new Array(n).fill(0);
  });

  for (let i = 0; i < n; i++) {
    let currentBot = 0.0;
    timeline.core_variants.forEach(v => {
      const frs = (trajectories[v] && trajectories[v][i]) || 0.0;
      const top = currentBot + frs;
      stackedBots[v][i] = currentBot;
      stackedTops[v][i] = Math.min(1.0, top);
      currentBot = top;
    });
    // WT remainder fills up to 1.0
    stackedBots["__WT__"][i] = 1.0 - timeline.wt_series[i];
    stackedTops["__WT__"][i] = 1.0;
  }

  // Draw ribbons
  allLayers.forEach(layer => {
    let pathD = "";
    // Top boundary left-to-right
    for (let i = 0; i < n; i++) {
      const x = getX(i);
      const y = getY(stackedTops[layer][i]);
      if (i === 0) pathD += `M ${x} ${y}`;
      else {
        const prevX = getX(i - 1);
        const prevY = getY(stackedTops[layer][i - 1]);
        const midX = (prevX + x) / 2;
        pathD += ` C ${midX} ${prevY}, ${midX} ${y}, ${x} ${y}`;
      }
    }
    // Bottom boundary right-to-left
    for (let i = n - 1; i >= 0; i--) {
      const x = getX(i);
      const y = getY(stackedBots[layer][i]);
      if (i === n - 1) pathD += ` L ${x} ${y}`;
      else {
        const nextX = getX(i + 1);
        const nextY = getY(stackedBots[layer][i + 1]);
        const midX = (nextX + x) / 2;
        pathD += ` C ${midX} ${nextY}, ${midX} ${y}, ${x} ${y}`;
      }
    }
    pathD += " Z";

    const path = document.createElementNS("http://www.w3.org/2000/svg", "path");
    path.setAttribute("d", pathD);
    path.setAttribute("fill", layer === "__WT__" ? "#afc3ca" : variantColors[layer]);
    path.setAttribute("opacity", layer === "__WT__" ? "0.35" : "0.85");
    path.setAttribute("stroke", "#ffffff");
    path.setAttribute("stroke-width", "1");
    ribbonsGroup.appendChild(path);
  });

  // Draw visit markers & date labels
  visits.forEach((v, idx) => {
    const x = getX(idx);
    const line = document.createElementNS("http://www.w3.org/2000/svg", "line");
    line.setAttribute("x1", x);
    line.setAttribute("y1", topY);
    line.setAttribute("x2", x);
    line.setAttribute("y2", botY);
    line.setAttribute("stroke", "#94a3b8");
    line.setAttribute("stroke-width", "1.5");
    line.setAttribute("stroke-dasharray", "2,2");
    markersGroup.appendChild(line);

    const circle = document.createElementNS("http://www.w3.org/2000/svg", "circle");
    circle.setAttribute("cx", x);
    circle.setAttribute("cy", botY);
    circle.setAttribute("r", "4");
    circle.setAttribute("fill", "#0f172a");
    markersGroup.appendChild(circle);

    const text = document.createElementNS("http://www.w3.org/2000/svg", "text");
    text.setAttribute("x", x);
    text.setAttribute("y", botY + 18);
    text.setAttribute("font-size", "10");
    text.setAttribute("fill", "#475569");
    text.setAttribute("text-anchor", "middle");
    text.setAttribute("font-weight", "600");
    text.textContent = v.date_display;
    markersGroup.appendChild(text);
  });

  // Active visit update function
  let cursorAnimation;
  function updateActiveVisit(idx, play = false) {
    if (idx < 0) idx = 0;
    if (idx >= n) idx = n - 1;
    const current = visits[idx];
    const x = getX(idx);

    const previousX = Number(cursor.getAttribute("x1"));
    cursorAnimation?.cancel();
    cursor.setAttribute("x1", x);
    cursor.setAttribute("x2", x);
    if (play && previousX !== x && !matchMedia("(prefers-reduced-motion: reduce)").matches) {
      cursorAnimation = cursor.animate(
        [{ transform: `translateX(${previousX - x}px)` }, { transform: "translateX(0px)" }],
        { duration: 200, easing: "cubic-bezier(.23, 1, .32, 1)" }
      );
    }
    scrubber.value = idx;
    activeLabel.textContent = `${idx + 1} of ${n}`;
    btnPrev.disabled = idx === 0;
    btnNext.disabled = idx === n - 1;
    visitButtons.forEach((button, index) => button.setAttribute("aria-pressed", String(index === idx)));

    // Update active visit summary card
    document.getElementById("visit-date-title").textContent = current.date_display;
    document.getElementById("visit-sample-id").textContent = current.sample_id;
    const link = document.getElementById("visit-sample-link");
    const batch = dataScript.dataset.batch;
    link.href = "/result?" + new URLSearchParams({ sample_id: current.sample_id, batch });

    const callBadge = document.getElementById("visit-call-badge");
    callBadge.textContent = current.call_label || current.call;
    callBadge.className = `status status-${current.call}`;

    const probDisp = document.getElementById("visit-prob-display");
    if (typeof current.probability === "number") {
      probDisp.textContent = `${(current.probability * 100).toFixed(1)}%`;
    } else {
      probDisp.textContent = "—";
    }

    // Populate variant bars
    const barsContainer = document.getElementById("visit-variant-bars");
    barsContainer.replaceChildren();
    if (current.variants && current.variants.length > 0) {
      current.variants.forEach(v => {
        const row = document.createElement("div");
        row.className = "frequency-row";

        const label = document.createElement("code");
        label.textContent = `${v.gene}@${v.mutation}`;

        const barTrack = document.createElement("div");
        barTrack.className = "frequency-track";

        const barFill = document.createElement("div");
        barFill.className = "frequency-fill";
        barFill.style.width = `${Math.min(100, Math.max(0, v.frs * 100))}%`;
        const mutKey = `${v.gene}@${v.mutation}`;
        barFill.style.background = variantColors[mutKey] || colorPalette[0];

        barTrack.appendChild(barFill);

        const val = document.createElement("strong");
        val.textContent = `${(v.frs * 100).toFixed(1)}%`;

        row.appendChild(label);
        row.appendChild(barTrack);
        row.appendChild(val);
        barsContainer.appendChild(row);
      });
    } else {
      const note = document.createElement("p");
      note.className = "fine";
      note.textContent = "No retained variants at this visit.";
      barsContainer.appendChild(note);
    }
    const alerts = document.getElementById("active-visit-alerts");
    alerts.replaceChildren(...timeline.alerts.filter(alert => alert.visit_idx === idx).map(alert => {
      const item = document.createElement("li");
      item.textContent = alert.message;
      return item;
    }));
  }

  // Interactive scrubber events
  scrubber.addEventListener("input", () => {
    stopPlayback();
    updateActiveVisit(parseInt(scrubber.value, 10));
  });

  btnPrev.addEventListener("click", () => {
    stopPlayback();
    updateActiveVisit(parseInt(scrubber.value, 10) - 1);
  });

  btnNext.addEventListener("click", () => {
    stopPlayback();
    updateActiveVisit(parseInt(scrubber.value, 10) + 1);
  });

  for (const button of visitButtons) {
    button.addEventListener("click", () => {
      stopPlayback();
      updateActiveVisit(Number(button.dataset.visitIndex));
    });
  }

  svg.addEventListener("click", event => {
    stopPlayback();
    const start = svg.createSVGPoint();
    start.x = event.clientX;
    start.y = event.clientY;
    const x = start.matrixTransform(svg.getScreenCTM().inverse()).x;
    const visit = visits.reduce((nearest, _, index) =>
      Math.abs(getX(index) - x) < Math.abs(getX(nearest) - x) ? index : nearest, 0);
    updateActiveVisit(visit);
  });

  // Play / Pause timer
  let playTimer = null;
  function stopPlayback() {
    clearInterval(playTimer);
    playTimer = null;
    btnPlay.textContent = "Play";
    btnPlay.setAttribute("aria-pressed", "false");
  }
  btnPlay.addEventListener("click", () => {
    if (playTimer) {
      stopPlayback();
    } else {
      if (Number(scrubber.value) === n - 1) updateActiveVisit(0);
      btnPlay.textContent = "Pause";
      btnPlay.setAttribute("aria-pressed", "true");
      playTimer = setInterval(() => {
        const next = Number(scrubber.value) + 1;
        updateActiveVisit(next, true);
        if (next >= n - 1) stopPlayback();
      }, 1200);
    }
  });

  // Initialize on first visit
  updateActiveVisit(0);
  document.querySelector(".visit-picker").hidden = false;
}

document.addEventListener("DOMContentLoaded", () => {
  initLongitudinalTimeline();
});
