var vp_width = 1000;
var vp_height = 1000;
var vp_padding = 20;

var stat_data = null;

var embed_x = function(v){ return v % stat_data.g_emb.width; };
var embed_y = function(v){ return Math.floor(v / stat_data.g_emb.width); };

var dx = [-1, 0, 1, -1, 1, -1, 0, 1];
var dy = [-1, -1, -1, 0, 0, 1, 1, 1];

var make_xscale = function(){
  return d3.scaleLinear()
    .domain([0, stat_data.g_emb.width - 1])
    .range([vp_padding, vp_width - vp_padding]);
};
var make_yscale = function(){
  return d3.scaleLinear()
    .domain([0, stat_data.g_emb.height - 1])
    .range([vp_padding, vp_height - vp_padding]);
};

var render_neighboring_groups = function(gid){
  var w = stat_data.g_emb.width, h = stat_data.g_emb.height;
  var n = stat_data.g.num_vertices, n_emb = w * h;
  var group = stat_data.solution[gid];
  var mapping = new Array(n_emb);
  for(var i = 0; i < n_emb; ++i){ mapping[i] = -1; }
  for(var i = 0; i < n; ++i){
    stat_data.solution[i].forEach(function(x){ mapping[x] = i; });
  }
  var neighboring = new Array(n);
  for(var i = 0; i < n; ++i){ neighboring[i] = false; }
  stat_data.g.edges.forEach(function(e){
    if(e[0] == gid){ neighboring[e[1]] = true; }
    if(e[1] == gid){ neighboring[e[0]] = true; }
  });
  var neighboring_emb = new Array(n);
  for(var i = 0; i < n; ++i){ neighboring_emb[i] = false; }
  group.forEach(function(u){
    var uy = embed_y(u), ux = embed_x(u);
    for(var d = 0; d < 8; ++d){
      var vy = uy + dy[d], vx = ux + dx[d];
      if(vy < 0 || h <= vy || vx < 0 || w <= vx){ continue; }
      var v = vy * w + vx;
      if(mapping[v] >= 0 && mapping[v] != gid){
        neighboring_emb[mapping[v]] = true;
      }
    }
  });

  var valid_edges = [];
  for(var uy = 0; uy < h; ++uy){
    for(var ux = 0; ux < w; ++ux){
      var u = uy * w + ux;
      if(mapping[u] < 0 || !neighboring[mapping[u]]){ continue; }
      for(var d = 0; d < 8; ++d){
        var vy = uy + dy[d], vx = ux + dx[d];
        if(vy < 0 || h <= vy || vx < 0 || w <= vx){ continue; }
        var v = vy * w + vx;
        if(u < v && mapping[u] == mapping[v]){
          valid_edges.push({
            x1: ux, y1: uy, x2: vx, y2: vy, c: neighboring_emb[mapping[u]]
          })
        }
      }
    }
  }
  var xscale = make_xscale(), yscale = make_yscale();
  d3.select('#neighboring-edges').selectAll('line').remove();
  d3.select('#neighboring-edges')
    .selectAll('line')
    .data(valid_edges)
    .enter()
    .append('line')
    .attr('x1', function(d){ return xscale(d.x1); })
    .attr('y1', function(d){ return xscale(d.y1); })
    .attr('x2', function(d){ return xscale(d.x2); })
    .attr('y2', function(d){ return xscale(d.y2); })
    .attr('stroke-width', 6.0)
    .attr('stroke', function(d){ return d.c ? '#c0ffc0' : '#ffc0c0'; });

  var vertices = [];
  stat_data.solution.forEach(function(g, index){
    if(neighboring[index]){
      g.forEach(function(v){
        var x = embed_x(v), y = embed_y(v);
        vertices.push({ x: x, y: y, c: neighboring_emb[index] });
      });
    }
  });
  d3.select('#neighboring-vertices').selectAll('circle').remove();
  d3.select('#neighboring-vertices')
    .selectAll('line')
    .data(vertices)
    .enter()
    .append('circle')
    .attr('cx', function(d){ return xscale(d.x); })
    .attr('cy', function(d){ return yscale(d.y); })
    .attr('r', 4.0)
    .attr('stroke', function(d){ return d.c ? '#00c000' : '#c00000'; })
    .attr('stroke-width', 3.0)
    .attr('fill', 'white');
};
var clear_neighboring_groups = function(){
  d3.select('#neighboring-edges').selectAll('line').remove();
  d3.select('#neighboring-vertices').selectAll('circle').remove();
};

var render_group_edges = function(){
  var w = stat_data.g_emb.width, h = stat_data.g_emb.height;
  var n = stat_data.g.num_vertices, n_emb = w * h;
  var mapping = new Array(n_emb);
  for(var i = 0; i < n_emb; ++i){ mapping[i] = -1; }
  for(var i = 0; i < n; ++i){
    stat_data.solution[i].forEach(function(x){ mapping[x] = i; });
  }
  var valid_edges = [];
  for(var uy = 0; uy < h; ++uy){
    for(var ux = 0; ux < w; ++ux){
      var u = uy * w + ux;
      for(var d = 0; d < 8; ++d){
        var vy = uy + dy[d], vx = ux + dx[d];
        if(vy < 0 || h <= vy || vx < 0 || w <= vx){ continue; }
        var v = vy * w + vx;
        if(u < v && mapping[u] >= 0 && mapping[u] == mapping[v]){
          valid_edges.push({ x1: ux, y1: uy, x2: vx, y2: vy, g: mapping[u] })
        }
      }
    }
  }
  var xscale = make_xscale(), yscale = make_yscale();
  d3.select('#embed-edges').selectAll('line').remove();
  d3.select('#embed-edges')
    .selectAll('line')
    .data(valid_edges)
    .enter()
    .append('line')
    .attr('x1', function(d){ return xscale(d.x1); })
    .attr('y1', function(d){ return xscale(d.y1); })
    .attr('x2', function(d){ return xscale(d.x2); })
    .attr('y2', function(d){ return xscale(d.y2); })
    .attr('stroke-width', 6.0)
    .attr('stroke', 'black')
    .attr('opacity', 0.2)
    .on('mouseover', function(d){
      render_current_group(d.g);
      render_neighboring_groups(d.g);
    })
    .on('mouseout', function(){
      clear_neighboring_groups();
      clear_current_group();
    });
};

var render_current_group = function(gid){
  var w = stat_data.g_emb.width, h = stat_data.g_emb.height;
  var g = stat_data.solution[gid], gset = new Set(g);
  var valid_edges = [];
  g.forEach(function(u){
    var uy = embed_y(u), ux = embed_x(u);
    for(var d = 0; d < 8; ++d){
      var vy = uy + dy[d], vx = ux + dx[d];
      if(vy < 0 || h <= vy || vx < 0 || w <= vx){ continue; }
      if(gset.has(vy * w + vx)){
        valid_edges.push({ x1: ux, y1: uy, x2: vx, y2: vy })
      }
    }
  });
  var xscale = make_xscale(), yscale = make_yscale();
  d3.select('#current-edges').selectAll('line').remove();
  d3.select('#current-edges')
    .selectAll('line')
    .data(valid_edges)
    .enter()
    .append('line')
    .attr('x1', function(d){ return xscale(d.x1); })
    .attr('y1', function(d){ return xscale(d.y1); })
    .attr('x2', function(d){ return xscale(d.x2); })
    .attr('y2', function(d){ return xscale(d.y2); })
    .attr('stroke-width', 6.0)
    .attr('stroke', 'black');
};
var clear_current_group = function(){
  d3.select('#current-edges').selectAll('line').remove();
};

var render_scored_edges = function(){
  var w = stat_data.g_emb.width, h = stat_data.g_emb.height;
  var n = stat_data.g.num_vertices, n_emb = w * h;
  var mapping = new Array(n_emb);
  for(var i = 0; i < n_emb; ++i){ mapping[i] = -1; }
  for(var i = 0; i < n; ++i){
    stat_data.solution[i].forEach(function(x){ mapping[x] = i; });
  }
  var edges = new Array(n);
  for(var i = 0; i < n; ++i){ edges[i] = new Set(); }
  stat_data.g.edges.forEach(function(e){
    edges[e[0]].add(e[1]);
    edges[e[1]].add(e[0]);
  });
  var valid_edges = [];
  for(var uy = 0; uy < h; ++uy){
    for(var ux = 0; ux < w; ++ux){
      var u = uy * w + ux;
      if(mapping[u] < 0){ continue; }
      for(var d = 0; d < 8; ++d){
        var vy = uy + dy[d], vx = ux + dx[d];
        if(vy < 0 || h <= vy || vx < 0 || w <= vx){ continue; }
        var v = vy * w + vx;
        if(u > v || mapping[v] < 0 || mapping[u] == mapping[v]){ continue; }
        if(edges[mapping[u]].has(mapping[v])){
          valid_edges.push({ x1: ux, y1: uy, x2: vx, y2: vy });
          edges[mapping[u]].delete(mapping[v]);
          edges[mapping[v]].delete(mapping[u]);
        }
      }
    }
  }
  var xscale = make_xscale(), yscale = make_yscale();
  d3.select('#scored-edges').selectAll('line').remove();
  d3.select('#scored-edges')
    .selectAll('line')
    .data(valid_edges)
    .enter()
    .append('line')
    .attr('x1', function(d){ return xscale(d.x1); })
    .attr('y1', function(d){ return xscale(d.y1); })
    .attr('x2', function(d){ return xscale(d.x2); })
    .attr('y2', function(d){ return xscale(d.y2); })
    .attr('stroke-width', 1.5)
    .attr('stroke', '#c0c0ff');
};
/*
var render_embed_edges = function(){
  var valid_edges = []
  stat_data.g.edges.forEach(function(e){
    var transform = stat_data.solution
    var u = transform[e.u - 1], v = transform[e.v - 1];
    var ux = embed_x(u), uy = embed_y(u);
    var vx = embed_x(v), vy = embed_y(v);
    if(Math.abs(ux - vx) <= 1 && Math.abs(uy - vy) <= 1){
      valid_edges.push({ x1: ux, y1: uy, x2: vx, y2: vy, weight: e.w })
    }
  });
  var xscale = make_xscale(), yscale = make_yscale();
  d3.select('#embed-edges').selectAll('line').remove();
  d3.select('#embed-edges')
    .selectAll('line')
    .data(valid_edges)
    .enter()
    .append('line')
    .attr('x1', function(d){ return xscale(d.x1); })
    .attr('y1', function(d){ return xscale(d.y1); })
    .attr('x2', function(d){ return xscale(d.x2); })
    .attr('y2', function(d){ return xscale(d.y2); })
    .attr('stroke-width', 1.5)
    .attr('stroke', function(d){ return palette[d.weight]; });
};

var render_original_edges = function(c){
  var valid_edges = []
  stat_data.g.edges.forEach(function(e){
    if(e.u != c && e.v != c){ return; }
    var transform = stat_data.solution;
    var u = transform[e.u - 1], v = transform[e.v - 1];
    var ux = embed_x(u), uy = embed_y(u);
    var vx = embed_x(v), vy = embed_y(v);
    valid_edges.push({ x1: ux, y1: uy, x2: vx, y2: vy, weight: e.w })
  });
  var xscale = make_xscale(), yscale = make_yscale();
  d3.select('#original-edges').selectAll('line').remove();
  d3.select('#original-edges')
    .selectAll('line')
    .data(valid_edges)
    .enter()
    .append('line')
    .attr('x1', function(d){ return xscale(d.x1); })
    .attr('y1', function(d){ return xscale(d.y1); })
    .attr('x2', function(d){ return xscale(d.x2); })
    .attr('y2', function(d){ return xscale(d.y2); })
    .attr('stroke-width', 1.5)
    .attr('stroke', function(d){ return palette[d.weight]; });
  d3.select('#embed-edges')
    .attr('opacity', 0.2);
};

var clear_original_edges = function(){
  d3.select('#original-edges').selectAll('line').remove();
  d3.select('#embed-edges')
    .attr('opacity', 1.0);
};
*/
var render_embed_vertices = function(){
  var vertices = [];
  stat_data.solution.forEach(function(g, index){
    g.forEach(function(v){
      var x = embed_x(v), y = embed_y(v);
      vertices.push({ x: x, y: y, g: index });
    });
  });
  var xscale = make_xscale(), yscale = make_yscale();
  d3.select('#embed-vertices').selectAll('circle').remove();
  d3.select('#embed-vertices')
    .selectAll('line')
    .data(vertices)
    .enter()
    .append('circle')
    .attr('cx', function(d){ return xscale(d.x); })
    .attr('cy', function(d){ return yscale(d.y); })
    .attr('r', 4.0)
    .attr('stroke', 'black')
    .attr('stroke-width', 1.0)
    .attr('fill', 'white')
    .on('mouseover', function(d){
      render_current_group(d.g);
      render_neighboring_groups(d.g);
    })
    .on('mouseout', function(){
      clear_neighboring_groups();
      clear_current_group();
    });
};

$(document).on('dragenter', function(e){
  e.stopPropagation();
  e.preventDefault();
});

$(document).on('dragover', function(e){
  e.stopPropagation();
  e.preventDefault();
});

$(document).on('drop', function(e){
  e.stopPropagation();
  e.preventDefault();
  e = e.originalEvent;
  var files = e.dataTransfer.files;
  if(files.length == 0){ return; }
  var reader = new FileReader();
  reader.onload = (function(f){
    return function(e){
      stat_data = JSON.parse(e.target.result);
      render_group_edges();
      render_scored_edges();
      render_embed_vertices();
      $("#info-score").text(stat_data["score"]);
      $("#info-total").text(stat_data["total"]);
      $("#info-time").text(stat_data["time"]);
    }
  })(files[0]);
  reader.readAsText(files[0]);
});
