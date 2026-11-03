/* Static OS + RFS prediction calculator for M3/M4 Elastic-Net Cox models.
   Two-level selection: endpoint (OS / RFS) -> model (M3 / M4).
   All computation is client-side; no data leaves the browser. */
(function () {
  "use strict";
  const A = ASSETS;                 // { levels, rawmap, endpoints:{OS,RFS} }
  const MCOL = { M3: "#5B5FA2", M4: "#C98AA6" };
  let curEP = "OS";                 // default endpoint
  let curM  = "M3";                 // default model

  const nice = {
    LVI:"Lymphovascular invasion", PNI:"Perineural invasion", Diff:"Differentiation",
    T_stage:"T stage", N_stage:"N stage", Age:"Age", Sex:"Sex", BMI:"BMI",
    CA199:"CA19-9", CEA:"CEA",
    IHC_Treg:"IHC Treg", IHC_epiN:"IHC epithelial-nuclear", IHC_epiC:"IHC epithelial-cytoplasmic",
    RNA_Treg:"Stromal FOXP3 (SightRNA)", RNA_epiN:"Nuclear E-FOXP3 (SightRNA)", RNA_epiC:"Cytoplasmic E-FOXP3 (SightRNA)"
  };
  // endpoint-specific wording
  function epShort() { return A.endpoints[curEP].meta.endpoint_short; }      // "OS" / "RFS"
  function riskWord() { return curEP === "OS" ? "death risk" : "recurrence/progression risk"; }
  function endpointFull() { return A.endpoints[curEP].meta.endpoint; }

  function mdef() { return A.endpoints[curEP].models[curM]; }
  function tpm() { return A.endpoints[curEP].timepoints_months; }
  function tpd() { return mdef().timepoints_days; }

  // ---- helpers ----
  function s0at(day) {
    const g = mdef().baseline.t_days, s = mdef().baseline.S0;
    if (day <= g[0]) return s[0];
    if (day >= g[g.length - 1]) return s[s.length - 1];
    let lo = 0, hi = g.length - 1;
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (g[mid] <= day) lo = mid; else hi = mid; }
    const frac = (day - g[lo]) / (g[hi] - g[lo]);
    return s[lo] + frac * (s[hi] - s[lo]);
  }
  function rawValue(varName, optIdx) {
    if (A.rawmap[varName]) return A.rawmap[varName][optIdx];
    return optIdx;
  }
  function contributions() {
    const coef = mdef().coefficients, used = mdef().used_vars, out = [];
    used.forEach(function (v) {
      const sel = document.getElementById("f_" + v);
      const x = rawValue(v, parseInt(sel.value, 10));
      out.push({ v: v, contrib: coef[v] * x });
    });
    return out;
  }
  function linpredFrom(cs) { return cs.reduce((a, c) => a + c.contrib, 0); }
  function strataLabel(lp) {
    const t = mdef().lp_tertiles;
    if (lp < t[0]) return ["Low risk", "low"];
    if (lp < t[1]) return ["Intermediate risk", "med"];
    return ["High risk", "high"];
  }

  // ---- build input fields (grouped) ----
  function buildInputs() {
    const used = mdef().used_vars;
    const molec = new Set(["IHC_Treg","IHC_epiN","IHC_epiC","RNA_Treg","RNA_epiN","RNA_epiC"]);
    const clin = used.filter(v => !molec.has(v));
    const mol  = used.filter(v => molec.has(v));
    function fieldHTML(v) {
      const opts = A.levels[v].map((lab, i) => `<option value="${i}">${lab}</option>`).join("");
      return `<div class="field"><label for="f_${v}">${nice[v] || v}</label>
              <select id="f_${v}">${opts}</select></div>`;
    }
    const grpStyle = "font-size:12px;font-weight:700;color:#8B91A1;text-transform:uppercase;letter-spacing:.4px;margin:12px 0 2px;";
    let html = `<div style="${grpStyle};margin-top:0">Clinical &amp; histopathology</div>${clin.map(fieldHTML).join("")}`;
    const molTitle = mol.some(v => v.startsWith("IHC")) ? "Molecular (IHC + SightRNA)" : "Molecular (SightRNA)";
    if (mol.length) html += `<div style="${grpStyle}">${molTitle}</div>${mol.map(fieldHTML).join("")}`;
    document.getElementById("inputPanels").innerHTML = html;
    used.forEach(v => document.getElementById("f_" + v).addEventListener("change", compute));
  }

  // ---- update info cards / header text on endpoint change ----
  function refreshChrome() {
    const ep = A.endpoints[curEP];
    document.getElementById("icEndpoint").textContent = curEP === "OS" ? "Overall survival" : "Recurrence-free survival";
    document.getElementById("icEndpointX").textContent = curEP === "OS" ? "Death event" : "Progression or death";
    document.getElementById("icTP").textContent = tpm().join(" / ") + " mo";
    document.getElementById("cindexVal").textContent = ep.cindex[curM];
    document.getElementById("cindexNote").textContent = "Model " + curM + ", internal (OOB)";
    document.getElementById("subEndpoint").textContent =
      curEP === "OS" ? "overall-survival prediction (research tool)"
                     : "recurrence-free-survival prediction (research tool)";
    document.getElementById("curveTitle").textContent = "Predicted individual " + epShort() + " curve";
    document.getElementById("curveCap").innerHTML =
      "Line = predicted " + epShort() + " probability; markers = " + tpm().join(" / ") + "-month points.";
  }

  // ---- compute + render ----
  function compute() {
    const cs = contributions();
    const lp = linpredFrom(cs);
    const eLp = Math.exp(lp);
    const td = tpd(), tm = tpm();
    const cls = ["t6","t12","t18"];   // three color slots (reused for both endpoints)
    let cards = "";
    tm.forEach(function (mo, i) {
      const S = Math.pow(s0at(td[i]), eLp);
      const surv = S * 100, risk = (1 - S) * 100;
      cards += `<div class="rc ${cls[i]}"><div class="tp">${mo}-month ${epShort()}</div>
        <div class="risk">${surv.toFixed(1)}%</div>
        <div class="os">${riskWord()} ${risk.toFixed(1)}%</div></div>`;
    });
    document.getElementById("results").innerHTML = cards;

    const [slab, scls] = strataLabel(lp);
    document.getElementById("strata").innerHTML =
      `Risk stratum: <b class="tag ${scls}">${slab}</b>
       <span class="muted">&nbsp;(lp = ${lp.toFixed(3)}; training-cohort tertiles — reference only)</span>`;

    drawContrib(cs);
    drawCurve(eLp);
  }

  function drawContrib(cs) {
    const sorted = cs.slice().sort((a, b) => Math.abs(b.contrib) - Math.abs(a.contrib));
    const maxAbs = Math.max.apply(null, sorted.map(c => Math.abs(c.contrib)).concat([1e-9]));
    let html = "";
    sorted.forEach(function (c) {
      const a = Math.abs(c.contrib);
      const zero = a < 1e-9;
      const w = zero ? 0 : Math.max(3, Math.round(a / maxAbs * 100));
      const col = c.contrib > 0 ? "#C98AA6" : (c.contrib < 0 ? "#5B5FA2" : "#D9DDE3");
      const sign = c.contrib > 0 ? "+" : (c.contrib < 0 ? "\u2212" : "");
      const vtxt = zero ? "0" : (sign + a.toFixed(2));
      html += `<div class="crow"><div class="cn" title="${nice[c.v]||c.v}">${nice[c.v]||c.v}</div>
        <div class="cbar"><span style="width:${w}%;background:${col}"></span></div>
        <div class="cv">${vtxt}</div></div>`;
    });
    document.getElementById("contrib").innerHTML = html;
  }

  function drawCurve(eLp) {
    const g = mdef().baseline.t_days, s0 = mdef().baseline.S0;
    const xs = g.map(d => d / 30.4375);
    const ys = s0.map(v => Math.pow(v, eLp) * 100);
    const col = MCOL[curM];
    const tm = tpm(), td = tpd();
    const mx = tm, my = td.map(d => Math.pow(s0at(d), eLp) * 100);
    const lbl = "Predicted " + epShort();
    const line = { x: xs, y: ys, mode: "lines", line: { color: col, width: 3 },
                   name: lbl, hovertemplate: "%{x:.1f} mo: " + epShort() + " %{y:.1f}%<extra></extra>" };
    const pts = { x: mx, y: my, mode: "markers+text", marker: { color: col, size: 9, line:{color:"#fff",width:1.5} },
                  text: my.map((v,i)=>`  ${mx[i]}mo`), textposition: "middle right",
                  textfont: { size: 10, color: "#8B91A1" },
                  hovertemplate: "%{x} mo: " + epShort() + " %{y:.1f}%<extra></extra>", showlegend:false };
    Plotly.newPlot("curve", [line, pts], {
      margin: { l: 48, r: 14, t: 8, b: 40 },
      xaxis: { title: "Months since surgery", range: [0, Math.max.apply(null, xs) * 1.02], zeroline:false, gridcolor:"#EEF0F3" },
      yaxis: { title: epShort() + " probability (%)", range: [0, 100], zeroline:false, gridcolor:"#EEF0F3" },
      font: { family: "Inter, Source Sans 3, Liberation Sans, Arial, sans-serif", size: 11, color:"#2E3340" },
      paper_bgcolor: "#fff", plot_bgcolor: "#fff", showlegend: false
    }, { displayModeBar: false, responsive: true });
  }

  // ---- endpoint toggle ----
  document.querySelectorAll(".epbtn").forEach(function (btn) {
    btn.addEventListener("click", function () {
      document.querySelectorAll(".epbtn").forEach(b => b.classList.remove("active"));
      btn.classList.add("active");
      curEP = btn.getAttribute("data-ep");
      refreshChrome(); buildInputs(); compute();
    });
  });
  // ---- model toggle ----
  document.querySelectorAll(".mbtn").forEach(function (btn) {
    btn.addEventListener("click", function () {
      document.querySelectorAll(".mbtn").forEach(b => b.classList.remove("active"));
      btn.classList.add("active");
      curM = btn.getAttribute("data-m");
      refreshChrome(); buildInputs(); compute();
    });
  });

  // ---- init ----
  refreshChrome();
  buildInputs();
  compute();
})();
