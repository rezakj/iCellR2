
#' Teal interactive HTML for cluster.plot (static ggplot unchanged)
#' @keywords internal
.ic2_cluster_teal_pal <- function(n) {
  base <- c(
    "#0E7490", "#EA580C", "#15803D", "#BE123C", "#0891B2",
    "#A16207", "#DB2777", "#475569", "#0F766E"
  )
  if (n <= length(base)) base[seq_len(n)] else rep(base, length.out = n)
}

#' @keywords internal
.ic2_cluster_plot_html <- function(
    emb2,
    emb3 = NULL,
    meta,
    plot.type,
    dims2,
    dims3 = NULL,
    col.by = "clusters",
    out.name = NULL,
    anno.clust = TRUE
) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for interactive cluster.plot().")
  }
  stopifnot(nrow(emb2) >= 1L)
  cells <- rownames(emb2)
  meta <- meta[cells, , drop = FALSE]

  cl <- if ("cluster_id" %in% colnames(meta)) {
    as.character(.ic2_factor_sorted(meta$cluster_id))
  } else {
    rep("0", length(cells))
  }
  cd <- if ("condition_id" %in% colnames(meta)) {
    as.character(meta$condition_id)
  } else {
    rep("all", length(cells))
  }
  cd[is.na(cd) | !nzchar(cd)] <- "all"

  has_mito <- "mito.percent" %in% colnames(meta)
  has_g <- "nGenes" %in% colnames(meta)
  has_u <- "UMIs" %in% colnames(meta)
  mito <- if (has_mito) as.numeric(meta$mito.percent) else rep(NA_real_, length(cells))
  ngenes <- if (has_g) as.numeric(meta$nGenes) else rep(NA_real_, length(cells))
  umis <- if (has_u) as.numeric(meta$UMIs) else rep(NA_real_, length(cells))

  has3d <- !is.null(emb3) && ncol(emb3) >= 3L && nrow(emb3) == nrow(emb2)
  points <- vector("list", length(cells))
  for (i in seq_along(cells)) {
    pt <- list(
      x = unname(emb2[i, 1L]),
      y = unname(emb2[i, 2L]),
      cl = cl[[i]],
      cd = cd[[i]],
      m = if (is.finite(mito[[i]])) unname(mito[[i]]) else NULL,
      g = if (is.finite(ngenes[[i]])) unname(ngenes[[i]]) else NULL,
      u = if (is.finite(umis[[i]])) unname(umis[[i]]) else NULL
    )
    if (has3d) pt$z <- unname(emb3[i, 3L])
    points[[i]] <- pt
  }

  qrng <- function(v) {
    v <- as.numeric(v)
    v <- v[is.finite(v)]
    if (!length(v)) return(list(min = 0, max = 1))
    qs <- stats::quantile(v, probs = c(0.01, 0.99), names = FALSE, na.rm = TRUE)
    if (!is.finite(qs[1]) || !is.finite(qs[2]) || qs[2] <= qs[1]) {
      qs <- range(v, na.rm = TRUE)
    }
    list(min = unname(qs[1]), max = unname(qs[2]))
  }

  default_color <- "cluster"
  if (identical(col.by, "conditions")) default_color <- "condition"
  else if (identical(col.by, "mito.percent") || identical(col.by, "mito")) default_color <- "mito"
  else if (identical(col.by, "nGenes")) default_color <- "ngenes"
  else if (identical(col.by, "UMIs")) default_color <- "umis"
  else if (!identical(col.by, "clusters") && !identical(col.by, "monochrome") &&
           !identical(col.by, "cc") && col.by %in% colnames(meta)) {
    if (is.numeric(meta[[col.by]])) {
      ## continuous custom column not wired in chips  -  fall back to cluster
      default_color <- "cluster"
    }
  }

  cl_levels <- levels(.ic2_factor_sorted(cl))
  cd_levels <- sort(unique(cd))

  payload <- list(
    n = length(cells),
    reduction = plot.type,
    has3d = has3d,
    has_mito = has_mito,
    has_ngenes = has_g,
    has_umis = has_u,
    default_color = default_color,
    anno = isTRUE(anno.clust),
    axis2d = paste0(toupper(plot.type), " ", as.integer(dims2)),
    axis3d = if (has3d) paste0(toupper(plot.type), " ", as.integer(dims3)) else list(),
    clusters = as.list(cl_levels),
    conditions = as.list(cd_levels),
    scales = list(
      mito = qrng(mito),
      ngenes = qrng(ngenes),
      umis = qrng(umis)
    ),
    points = points
  )

  json <- jsonlite::toJSON(payload, dataframe = "rows", auto_unbox = TRUE, digits = 6, null = "null")
  title <- toupper(plot.type)
  html <- .ic2_cluster_plot_html_template(json, title = title)

  if (is.null(out.name) || !nzchar(as.character(out.name)[1L])) {
    out.name <- paste0("cluster_plot_", plot.type)
  }
  out.name <- path.expand(as.character(out.name)[1L])
  out.name <- sub("\\.html$", "", out.name, ignore.case = TRUE)
  file <- paste0(out.name, ".html")
  parent <- dirname(file)
  if (!dir.exists(parent) && !identical(parent, ".")) {
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  }
  writeLines(html, file, useBytes = TRUE)
  file
}

#' @keywords internal
.ic2_cluster_plot_html_template <- function(json_payload, title = "UMAP") {
  pre <- '<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>@@TITLE@@</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
*{box-sizing:border-box}
body{margin:0;background:#fff;color:#14242A;font-family:"Avenir Next","Segoe UI","Helvetica Neue",Arial,sans-serif}
.wrap{max-width:1100px;margin:0 auto;padding:22px 18px 48px}
h1{margin:0 0 4px;font-size:22px;letter-spacing:-.02em}
.meta{color:#5A6B72;font-size:13px;margin:0 0 14px;line-height:1.45}
.card{border:1px solid #E5E7EB;border-radius:14px;padding:10px 10px 4px;margin-bottom:16px;background:#fff}
.card-head{display:flex;align-items:center;justify-content:space-between;gap:8px;margin:4px 8px 6px}
.card-head h2{margin:0;font-size:14px;font-weight:700}
.dl-group{display:flex;gap:6px}
.dl{appearance:none;border:1px solid #D1D5DB;background:#F8FAFC;color:#14242A;border-radius:8px;padding:5px 9px;font:inherit;font-size:11.5px;cursor:pointer;font-weight:600}
.dl:hover{border-color:#0F766E;color:#0F766E}
.plot{width:100%;height:560px}
.plot3{width:100%;height:620px}
.row{display:flex;gap:8px;flex-wrap:wrap;margin:0 0 10px;align-items:center}
.lab{font-size:11px;color:#5A6B72;text-transform:uppercase;letter-spacing:.06em;margin-right:2px}
.chip{appearance:none;border:1px solid #D1D5DB;background:#fff;border-radius:999px;padding:5px 11px;font:inherit;font-size:12px;font-weight:600;cursor:pointer;color:#14242A}
.chip.active{background:#0F766E;border-color:#0F766E;color:#ECFEFF}
.hidden{display:none !important}
</style>
</head>
<body>
<div class="wrap">
  <h1 id="h1">@@TITLE@@</h1>
  <p class="meta" id="meta"></p>
  <div class="row"><span class="lab">Colour</span><span id="color-by" class="row" style="margin:0"></span></div>
  <div class="row" id="scale-row"><span class="lab">Scale</span>
    <button class="chip active" type="button" data-scale="clip">Clip 1-99%</button>
    <button class="chip" type="button" data-scale="full">Full range</button>
  </div>
  <div class="card">
    <div class="card-head">
      <h2 id="h2-2d">2D</h2>
      <div class="dl-group">
        <button class="dl" type="button" data-dl="p2" data-fmt="png">PNG</button>
        <button class="dl" type="button" data-dl="p2" data-fmt="svg">SVG</button>
      </div>
    </div>
    <div class="plot" id="p2"></div>
  </div>
  <div class="card" id="card-3d">
    <div class="card-head">
      <h2 id="h2-3d">3D</h2>
      <div class="dl-group">
        <button class="dl" type="button" data-dl="p3" data-fmt="png">PNG</button>
        <button class="dl" type="button" data-dl="p3" data-fmt="svg">SVG</button>
      </div>
    </div>
    <div class="plot3" id="p3"></div>
  </div>
</div>
<script>
window.IC2_RED = '
  post <- ';
</script>
<script>
(function(){
  const D=window.IC2_RED;
  const PAL=["#0E7490","#EA580C","#15803D","#BE123C","#0891B2","#A16207","#DB2777","#475569","#0F766E"];
  const INK="#14242A", MUTE="#5A6B72";
  const CONT={
    mito:{key:"m", title:"mito%", scale:"mito", colorscale:[[0,"#E2E8F0"],[0.5,"#FB923C"],[1,"#BE123C"]]},
    ngenes:{key:"g", title:"nGenes", scale:"ngenes", colorscale:[[0,"#E8F0F3"],[0.5,"#0E7490"],[1,"#134E4A"]]},
    umis:{key:"u", title:"UMIs", scale:"umis", colorscale:[[0,"#F3EFE8"],[0.5,"#EA580C"],[1,"#9A3412"]]}
  };
  const has3d=!!D.has3d;
  const red=D.reduction||"UMAP";
  const ax2lab=D.axis2d||["Dim1","Dim2"];
  const ax3lab=D.axis3d||["Dim1","Dim2","Dim3"];
  let colorBy=D.default_color||"cluster", scaleMode="clip";

  document.getElementById("h1").textContent=red.toUpperCase();
  document.getElementById("h2-2d").textContent="2D  |  "+ax2lab[0]+"-"+ax2lab[1].replace(/^[A-Za-z]+/, function(s){return s;});
  document.getElementById("h2-2d").textContent="2D  |  "+ax2lab.join("-");
  document.getElementById("h2-3d").textContent="3D  |  "+ax3lab.join("-");
  if(!has3d) document.getElementById("card-3d").classList.add("hidden");

  function setMeta(){
    document.getElementById("meta").textContent=
      D.n.toLocaleString()+" cells  |  "+D.clusters.length+" clusters  |  "+D.conditions.length+" samples  |  continuous colours use "+(scaleMode==="clip"?"1-99% quantile clip":"full range");
  }
  setMeta();

  const colorEl=document.getElementById("color-by");
  const opts=[["cluster","Cluster"],["condition","Condition"]];
  if(D.has_mito) opts.push(["mito","Mito %"]);
  if(D.has_ngenes) opts.push(["ngenes","nGenes"]);
  if(D.has_umis) opts.push(["umis","UMIs"]);
  opts.forEach(([id,lab])=>{
    const b=document.createElement("button"); b.type="button";
    b.className="chip"+(id===colorBy?" active":""); b.dataset.id=id; b.textContent=lab;
    b.onclick=()=>{
      colorBy=id;
      [...colorEl.querySelectorAll(".chip")].forEach(x=>x.classList.toggle("active",x.dataset.id===colorBy));
      document.getElementById("scale-row").style.display = CONT[colorBy] ? "" : "none";
      draw();
    };
    colorEl.appendChild(b);
  });
  document.getElementById("scale-row").style.display = CONT[colorBy] ? "" : "none";
  document.querySelectorAll("[data-scale]").forEach(b=>{
    b.onclick=()=>{
      scaleMode=b.dataset.scale;
      document.querySelectorAll("[data-scale]").forEach(x=>x.classList.toggle("active",x===b));
      setMeta(); draw();
    };
  });

  function levels(){return colorBy==="condition"?D.conditions:D.clusters;}
  function key(p){return colorBy==="condition"?p.cd:p.cl;}
  function centroids2(){const m={};D.points.forEach(p=>{const k=key(p);if(!m[k])m[k]={x:0,y:0,n:0};m[k].x+=p.x;m[k].y+=p.y;m[k].n++;});return Object.entries(m).map(([lab,v])=>({lab,x:v.x/v.n,y:v.y/v.n}));}
  function centroids3(){const m={};D.points.forEach(p=>{const k=key(p);if(!m[k])m[k]={x:0,y:0,z:0,n:0};m[k].x+=p.x;m[k].y+=p.y;m[k].z+=p.z;m[k].n++;});return Object.entries(m).map(([lab,v])=>({lab,x:v.x/v.n,y:v.y/v.n,z:v.z/v.n}));}
  function ax2(title){return{title:{text:title,font:{size:12,color:MUTE}},tickfont:{color:MUTE,size:10},showgrid:false,zeroline:false,showline:true,linecolor:INK,linewidth:1.2,ticks:"outside",mirror:false,color:INK};}
  function ax3(title){return{title:title,titlefont:{size:12,color:MUTE},tickfont:{color:MUTE,size:10},showgrid:true,gridcolor:"rgba(20,36,42,0.08)",zeroline:false,showline:true,linecolor:INK,backgroundcolor:"#fff",showspikes:false};}
  function contRange(cfg){
    if(scaleMode==="clip" && D.scales && D.scales[cfg.scale]) return [D.scales[cfg.scale].min, D.scales[cfg.scale].max];
    const vals=D.points.map(p=>p[cfg.key]).filter(v=>v!=null && isFinite(v));
    return [Math.min(...vals), Math.max(...vals)];
  }
  function discreteTraces(dims3){
    return levels().map((lab,i)=>{
      const pts=D.points.filter(p=>key(p)===lab);
      const t={type:dims3?"scatter3d":"scattergl", mode:"markers", name:String(lab),
        x:pts.map(p=>p.x), y:pts.map(p=>p.y),
        marker:{size:dims3?2.8:5.5, opacity:0.9, color:PAL[i%PAL.length], line:{width:0}},
        hovertemplate:String(lab)+"<br>(%{x:.3f}, %{y:.3f}"+(dims3?", %{z:.3f}":"")+")<extra></extra>"};
      if(dims3) t.z=pts.map(p=>p.z);
      return t;
    });
  }
  function contTrace(dims3){
    const cfg=CONT[colorBy], [cmin,cmax]=contRange(cfg);
    const t={
      type:dims3?"scatter3d":"scattergl", mode:"markers", name:cfg.title,
      x:D.points.map(p=>p.x), y:D.points.map(p=>p.y),
      marker:{size:dims3?2.8:5.5, opacity:0.9, color:D.points.map(p=>p[cfg.key]),
        colorscale:cfg.colorscale, cmin:cmin, cmax:cmax,
        colorbar:{title:cfg.title+(scaleMode==="clip"?"  |  1-99%":""), thickness:12, len:0.5, tickfont:{color:MUTE}},
        line:{width:0}},
      hovertemplate:cfg.title+" %{marker.color:.4g}<extra></extra>"
    };
    if(dims3) t.z=D.points.map(p=>p.z);
    return [t];
  }

  function draw(){
    const cont=!!CONT[colorBy];
    const t2=cont?contTrace(false):discreteTraces(false);
    const ann=(!cont && colorBy==="cluster")?centroids2().map(c=>({x:c.x,y:c.y,text:"<b>"+c.lab+"</b>",showarrow:false,font:{size:12,color:INK},bgcolor:"rgba(255,255,255,0.88)",borderpad:2})):[];
    Plotly.react("p2", t2, {
      paper_bgcolor:"#fff", plot_bgcolor:"#fff",
      font:{family:"Avenir Next, Segoe UI, Helvetica Neue, Arial, sans-serif", color:INK},
      margin:{l:56,r:18,t:36,b:52},
      title:{text:red.toUpperCase(), font:{size:16,color:INK}, x:0, xanchor:"left"},
      showlegend:!cont,
      legend:{orientation:"v", x:1.02, y:1, bgcolor:"rgba(0,0,0,0)", font:{size:11,color:MUTE}, itemsizing:"constant"},
      xaxis:Object.assign({},ax2(ax2lab[0])),
      yaxis:Object.assign({},ax2(ax2lab[1]),{scaleanchor:"x",scaleratio:1}),
      annotations:ann
    },{responsive:true, displaylogo:false, modeBarButtonsToRemove:["lasso2d","select2d"],
       toImageButtonOptions:{format:"png", filename:red+"_2d_"+colorBy, height:900, width:1100, scale:2}});

    if(!has3d) return;
    let t3=cont?contTrace(true):discreteTraces(true);
    if(!cont && colorBy==="cluster"){
      const cents=centroids3();
      t3=t3.concat([{type:"scatter3d", mode:"text", name:"labels",
        x:cents.map(c=>c.x), y:cents.map(c=>c.y), z:cents.map(c=>c.z),
        text:cents.map(c=>c.lab), textfont:{size:13,color:INK}, hoverinfo:"skip", showlegend:false}]);
    }
    Plotly.react("p3", t3, {
      paper_bgcolor:"#fff",
      font:{family:"Avenir Next, Segoe UI, Helvetica Neue, Arial, sans-serif", color:INK},
      margin:{l:0,r:0,t:36,b:0},
      title:{text:red.toUpperCase()+" 3D", font:{size:16,color:INK}, x:0.02, xanchor:"left"},
      showlegend:!cont,
      legend:{orientation:"v", x:0.98, y:0.98, xanchor:"right", bgcolor:"rgba(255,255,255,0.75)", font:{size:11,color:MUTE}},
      scene:{bgcolor:"#fff", xaxis:ax3(ax3lab[0]), yaxis:ax3(ax3lab[1]), zaxis:ax3(ax3lab[2]),
        camera:{eye:{x:1.45,y:1.45,z:1.05}}, aspectmode:"data"}
    },{responsive:true, displaylogo:false,
       toImageButtonOptions:{format:"png", filename:red+"_3d_"+colorBy, height:900, width:1100, scale:2}});
  }

  document.querySelectorAll("[data-dl]").forEach(btn=>{
    btn.addEventListener("click",()=>{
      const id=btn.getAttribute("data-dl");
      if(id==="p3" && !has3d) return;
      Plotly.downloadImage(id,{
        format: btn.getAttribute("data-fmt")||"png",
        width:1200, height:900, scale:2,
        filename: red+(id==="p3"?"_3d_":"_2d_")+colorBy
      });
    });
  });

  draw();
})();
</script>
</body>
</html>'
  pre <- gsub("@@TITLE@@", title, pre, fixed = TRUE)
  paste0(pre, json_payload, post)
}
