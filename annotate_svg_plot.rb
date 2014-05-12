#!/usr/bin/env ruby


require 'rexml/document'
include REXML

if ARGV.size < 2
  $stderr.puts "USAGE: #{$0} input_svg ratio_info [options]"
  exit -1
end
infn = ARGV.shift
ratiofn = ARGV.shift
options = ARGV
# parse plotting options
fetfn = nil
doAnnotation = false
doSequence = false
doZoom = false
doMarkers = false
marker_n = -1
doEnds = false
doId = false
sequenceId = ""
#color_scheme = "blue-red"
color_scheme=""
hmin = false
hmax = false
fade = false
doLegend = true
annotationtype = "unknown"
discrete_colors = Hash.new

while (opt = options.shift) do
	if (opt == "-a")
		doAnnotation = true
		annotationtype = options.shift
	elsif (opt == "-s")
		doSequence = true
	elsif (opt == "-z")
		doZoom = true
	elsif (opt == "-nolegend")
		doLegend = false
	elsif (opt == "-m")
		doMarkers = true
		opt = options.shift
		marker_n = opt.to_i
		if !(marker_n > 0)
			$stderr.puts "marker_n #{marker_n} must be a positive integer"
			exit -1
		end
	elsif (opt == "-e")
		doEnds = true
	elsif (opt == "-i")
		doId = true
		sequenceID = options.shift
	elsif (opt == "-colorscheme")
		color_scheme = options.shift
		if !( (color_scheme == "blue-red") || (color_scheme == "orange-red") || (color_scheme == "green-red") || (color_scheme == "discrete"))
			$stderr.puts "invalid color_scheme #{color_scheme}"
			exit -1
		end
		# store discrete colors if necessary
		if (color_scheme == "discrete")
			discrete_color_str = options.shift # should be of the format [comma-separated VALUES];[comma-separated COLORS in hex]
			vstr, cstr = discrete_color_str.split(";")
			values = vstr.split(",")
			colors = cstr.split(",")
			if values.length != colors.length
				$stderr.puts "Number of values (#{values.length}) does not match number of colors (#{colors.length}) using DISCRETE color scheme"
				exit -1
			end
			0.upto(values.length-1) do |i|
				discrete_colors[values[i].to_f] = colors[i]
			end
		end			
	elsif (opt == "-hmin")
		hmin = options.shift.to_f
	elsif (opt == "-hmax")
		hmax = options.shift.to_f
        elsif (opt == "-fade")
                fade = true
	else
		$stderr.puts "Invalid option #{opt}"
		exit -1
	end
end

svg = Document.new(File.new(infn))

svg.root.attributes['width'] = "100%"
svg.root.attributes['height'] = "100%"

js_txt = <<ENDSCRIPT

var currScale=1;
var currX=0.0;
var currY=0.0;
var currRot=0.0;
var readCountThickness=1;
var startPanX = 0.0;
var startPanY = 0.0;
var showSeq = 1;
grabbed = 0;

function handle_wheel(delta) {
  if (delta > 0) 
    zoom(1.1);
  else
    zoom(0.9);
}

function wheel(event) {
        var delta = 0;
        if (!event) /* For IE. */
                event = window.event;
        if (event.wheelDelta) { /* IE/Opera. */
                delta = event.wheelDelta/120;
                /** In Opera 9, delta differs in sign as compared to IE.
                 */
                if (window.opera)
                        delta = -delta;
        } else if (event.detail) { /** Mozilla case. */
                /** In Mozilla, sign of delta is different than in IE.
                 * Also, delta is multiple of 3.
                 */
                delta = -event.detail/3;
        }
        /** If delta is nonzero, handle it.
         * Basically, delta is now positive if wheel was scrolled up,
         * and negative, if wheel was scrolled down.
         */
        if (delta)
                handle_wheel(delta);
        /** Prevent default actions caused by mouse wheel.
         * That might be ugly, but we handle scrolls somehow
         * anyway, so don't bother here..
         */
        if (event.preventDefault)
                event.preventDefault();
	event.returnValue = false;
}


function init() {
  var o = document.getElementById("outline");
<!--  o.setAttribute('style', 'stroke:#999999; fill:none; stroke-width:1.5') -->
  o = document.getElementById("pairs");
  o.setAttribute('style', 'stroke:#999999; stroke-width:1')
  grabbed = 0;

  if (window.addEventListener)
    /** DOMMouseScroll is for mozilla. */
    window.addEventListener('DOMMouseScroll', wheel, false);
  /** IE/Opera. */
  window.onmousewheel = document.onmousewheel = wheel;


  reset();
}


function reset() {
 currScale=1.0;
 currX=0.0;
 currY=0.0;
 currRot=0.0;
 readCountThickness=4;
 startPanX = 0.0;
 startPanY = 0.0;
 showSeq = 1;
 ungrab();

 var d = document.getElementById("draghelp");
 d.removeAttribute('display');

 updateView();
}

function zoom(k)
{
 currScale *= k;
 updateView();
}

function rotate(a)
{
 currRot += a;
 updateView();
}

function grab(evt){
  var g=document.getElementById("grabber"); 
  g.setAttribute("onmousemove", "slide(evt)");
  startPanX = evt.clientX;
  startPanY = evt.clientY;

  var s=document.getElementById("seq");
  s.setAttribute('display', 'none');

  var d = document.getElementById("draghelp");
  d.setAttribute('display', 'none');
 
  grabbed = 1;

  evt.preventDefault();
}

function ungrab(){
  var g=document.getElementById("grabber"); 
  g.setAttribute("onmousemove", "");
  var s=document.getElementById("seq");
  if(showSeq) {
    s.removeAttribute('display')
  }

  grabbed = 0;
}

function slide(evt){
       currX += (evt.clientX - startPanX) / currScale;
       currY += (evt.clientY - startPanY) / currScale;
       startPanX = evt.clientX;
       startPanY = evt.clientY;
       updateView();
}

function updateView() {
   var r = document.getElementById("annotation");
   r.setAttribute("stroke-width", readCountThickness)

   var o=document.getElementById("main");    
   o.setAttribute("transform", "scale(" + currScale + ") translate(" + currX + ", " + currY + ") rotate(" + currRot + ")");

   o = document.getElementById("seq");
   if (showSeq == 0) {
      o.setAttribute('display', 'none');
   } else {
      if (!grabbed) {
        o.removeAttribute('display');
      }
   }
  var child = o.firstChild;
   while (child) {
   	if (child.nodeName == "text") {
   		var ox = child.getAttributeNS(null, "x");
   		var oy = child.getAttributeNS(null, "y");
		 	child.setAttributeNS(null, "transform", "rotate(" + currRot*-1.0 + " " + ox + " " + oy + ")");
		}
		child = child.nextSibling;
	}
}

function highlight(o) {
  o.setAttribute("opacity", "1.0");
}

function unhighlight(o) {
  o.setAttribute("opacity", "0.4");
}

function changeReadCountThickness(d) {
  readCountThickness += d;
  if (readCountThickness < 1)
    readCountThickness = 1;
  updateView();
}

function toggleShowSeq(c) {
  if (c) {
     showSeq = 1;
  } else {
     showSeq = 0;
  }
  updateView();
    
}

ENDSCRIPT

if (doZoom)
	svg.root.add_attribute('onload', 'init()')
	svg.root.elements[1] = Element.new('script');
	svg.root.elements[1].add_attribute('type', 'text/ecmascript')
	svg.root.elements[1].text = "\n"
	svg.root.elements[1].add(CData.new(js_txt))
end

t = svg.root.elements[3].attribute('transform').to_s
curscale = /scale\((.+?),(.+?)\)/.match(t).captures.map { |s| s.to_f }

### view area

viewframe = svg.root.add_element('svg', {'x' => '0', 'y' => '0', 'width' => '100%', 'height'=>"100%"})

# help text
#txt = viewframe.add_element('text', {'x'=>'25', 'y'=>'150', 'text-anchor'=>'left', 'fill'=>'#777777', 'id'=>'draghelp'})
#txt.text = "Drag to move view"

# move structure contents to 'main' element
main = viewframe.add_element('g', {'id'=>'main', 'transform'=>'scale(1.0), translate(0.0, 0.0)'})

img = svg.root.elements[3].deep_clone
svg.root.elements.delete(svg.root.elements[3])

main.elements.add(img)

# move the structure model so it doesnt overlap w/ the legend area
imgattrib = img.attribute('transform').to_s
if imgattrib =~ /scale\((.*),(.*)\)\s+translate\((.*),(.*)\)/
	xsc = $1.to_f
	ysc = $2.to_f
	xtr = $3.to_f*1.2
	ytr = $4.to_f*1.2
	img.delete_attribute('transform')
	img.add_attribute('transform', "scale(#{xsc},#{ysc}) translate(#{xtr},#{ytr})")
else
	$stderr.puts "ERROR: failed to parse img transform #{imgattrib}"
	exit -1
end

# mouse grab rect
grabber = viewframe.add_element('rect', {'id'=>'grabber', 'x'=>0, 'y'=>0, 'width'=>"100%", 'height'=>'100%', 'style'=>"fill:white; opacity:0;", 'onmousedown'=>'grab(evt)', 'onmouseup'=>'ungrab()'})


### legend area
if color_scheme == "discrete"
	boxheight = 15*(discrete_colors.length+4)
else
	boxheight = 95
end
navframe = svg.root.add_element('svg', {'x' => '0', 'y' => '0', 'width' => "130", 'height'=>boxheight})
navframe.add_element('rect', {'x'=>0, 'y'=>0, 'width'=>'100%', 'height'=>'100%', 'style'=>"fill:none; stroke:none"})
idframe = svg.root.add_element('svg', {'x' => '150', 'y' => '0', 'width' => "500", 'height'=>60})
idframe.add_element('rect', {'x'=>0, 'y'=>0, 'width'=>'500', 'height'=>'100%', 'style'=>"fill:none; stroke:none"})

# extract structure coordinates
backbone = img.elements[1]
outline = backbone.attribute('points').to_s.split(/\n/).reject { |s| s !~ /,/ }.map { |s| s.split(/,/).map { |p| p.to_f } }
outlineN = outline.length-1

# remove the translate on the sequence; replace it vertical translate and horizontal centered text
seq = img.elements[3]
seq.delete_attribute('transform')
seq.add_attribute('transform', 'translate(0,5.6)')
seq.elements.each { |ele|
	ele.add_attribute('text-anchor', 'middle')
	#ele.add_attribute('dominant-baseline', 'central')
}

### compute norms of structure backbone (vectors perpendicular to backbone)
def norm(v)
  len = Math.sqrt(v[0]**2 + v[1]**2)
  [v[0] / len, v[1] / len] 
end
def perp(u)
  [-u[1], u[0]]
end
def midpoint(p0, p1)
  [(p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2]
end
def add(p0,p1)
  [p0[0] + p1[0], p0[1] + p1[1]]
end
def scalar_mult(u, a)
  [u[0]*a, u[1]*a]
end
norms = []
1.upto(outline.length-1) do |i| 
  p0 = outline[i-1]
  p1 = outline[i]
  n = norm(perp(add(p1,scalar_mult(p0,-1.0))))
  m = midpoint(p0,p1)
  norms << n
end
scale = 10

if doEnds
	# add 5' and 3' indicators
	fivep_back = outline.first
	threep_back = outline.last	              
	fivep_coord = [outline[0][0] - (outline[1][0] - outline[0][0]), 
		             outline[0][1] - (outline[1][1] - outline[0][1]) ]
	threep_coord = [outline[outlineN][0] + (outline[outlineN][0] - outline[outlineN-1][0]),
		              outline[outlineN][1] + (outline[outlineN][1] - outline[outlineN-1][1])]
		              
	#seq.delete_attribute('style')
	#seq.add_attribute('style', 'font-family: SansSerif; font-weight: bold')
	ends = img.add_element('g', {'transform'=>seq.attributes['transform'], 'id'=>'ends', 'style'=>'font-family: SansSerif'})
	fivep = ends.add_element('text', {'id'=>'fivep', 'x'=>"#{fivep_coord[0]}", 'y'=>"#{fivep_coord[1]}", 'font-size'=>'10'})
	fivep.text = "5'"
	threep = ends.add_element('text', {'id'=>'threep', 'x'=>"#{threep_coord[0]}", 'y'=>"#{threep_coord[1]}", 'font-size'=>'10'})
	threep.text = "3'"
end

if doMarkers
	i = marker_n
	markers = img.add_element('g', {'transform'=>'translate(1.3,-0.25)', 'id'=>'markers', 'style'=>'font-family: SansSerif'})
	while (i-1 <= outlineN) do
    n = if (i-1) <= 0 
        norms.first
      elsif (i-1) >= outlineN
        norms.last
      else  # average neighboring norms
        scalar_mult(add(norms[i-1], norms[i-2]), 0.5)
      end
		marker_coord = add(outline[i-1], scalar_mult(n, scale*3))
		marker_line_start_coord = add(marker_coord, scalar_mult(n, scale*-1))
		marker_line_end_coord = add(outline[i-1], scalar_mult(n, scale*1))
#		puts marker_coord, marker_line_start_coord, marker_line_end_coord
#		exit -1
		
		marker = markers.add_element('text', {'id'=>"marker#{i}", 'x'=>"#{marker_coord[0]}", 'y'=>"#{marker_coord[1]}", 'font-size'=>'8', 'style'=>"text-anchor: middle"})
		marker.text = "#{i}"
		markers.add_element('line', {"x1"=>"#{marker_line_start_coord[0]}", "y1"=>"#{marker_line_start_coord[1]}", "x2"=>"#{marker_line_end_coord[0]}", "y2"=>"#{marker_line_end_coord[1]}", "style"=>"stroke:rgb(0,0,0);stroke-width:0.5"})
		i += marker_n
	end
end

if doId
#	seqid = svg.root.add_element('sequenceID')
	g = idframe.add_element('svg', {'width'=>'100%', 'height'=>'100%', 'viewBox'=>"0 0 850 95", 'preserveAspectRatio' => 'none'})
	txt = g.add_element('text', {"x"=>"5", "y"=>45, "font-size"=>"16pt", "font-weight"=>"bold"})
	if (sequenceID.length > 50)
		txt.text = "#{sequenceID[0..49]}..."
	else
		txt.text = sequenceID
	end
end

if doAnnotation
	Read = Struct.new(:level, :end_pos, :count, :read_id)
	ratio_ind = Array.new
	i = 0
	seq.each { |x|
		line = x.to_s;
		if line =~ />([ATUGCNWGRatugcnwgr-])</
			if $1 == "-"
				# do nothing
			else
				ratio_ind << i
			end
			i+=1
		end
	}
	ratio = Array.new
	File.open(ratiofn).each_line do |line|
		next if line.chomp == ""
	#  puts "reading in #{line.chomp}"
		ratio << line.chomp.to_f
	end

	if ratio.length != ratio_ind.length
		$stderr.puts "ratio length #{ratio.length} does not match index length #{ratio_ind.length}!"
		exit -1
	end
  
#         def compute_opacity(a, b, x, hmin, hmax, fade)
#           sx = if hmin && hmax
#                  if x > hmax
#                    sx = (b-hmax <= 0) ? 0.5 : Float((x-hmax)/(b-hmax))
#                  elsif x < hmin
#                    sx = (hmin-a <= 0) ? 1 : Float((x-a)/(hmin-a))
#                  else
#                    sx = 0.25
#                  end
#                else
#                  sx = (b-a == 0) ? 0.5 : Float((x-a)/(b-a))
#                end
#           if fade
#             0.1+sx*0.8  # 0.1-0.9
#           else
#             if hmin && hmax
#               [0.25+sx*0.5, 0.75].max
#             else
#               0.75
#             end
#           end
#         end


	def heat_scale(a, b, x, color_scheme, hmin, hmax, fade)
		# x     R G B
		# 0     0 0 1
		# 0.25  0 1 1
		# 0.5   0 1 0
		# 0.75  1 1 0
		# 1     1 0 0
	
		# blue colorscale for negative x, red colorscale for positive x
		# grey if within [hmin, hmax] (no information)
		if color_scheme == "blue-red" && hmin && hmax
			if x >= hmin && x <= hmax
				return "NA" # "#EEEEEE-0.25"
			elsif x > hmax
				sx = (b-hmax <= 0) ? 0.5 : Float((x-hmax)/(b-hmax))
        opacity = if fade
            0.1+sx*0.8  # 0.1-0.9
          else
            [0.25+sx*0.5, 0.75].max
          end
				red, green, blue = [0,0,0]
				red = 1
				# [0,1] -> [1,0]
				green = 1-sx
				red *= 255
				green *= 255
				blue *= 255
				heatstr = sprintf("#%02x%02x%02x", red, green, blue)
				return "#{heatstr}-#{opacity}"
			else
				sx = (hmin-a <= 0) ? 1 : Float((x-a)/(hmin-a))
        opacity = if fade
            0.1+sx*0.8  # 0.1-0.9
          else
            [0.25+sx*0.5, 0.75].max
          end
				red, green, blue = [0,0,0]
				blue = 1
				# [0,1]
				green = sx
				red *= 255
				green *= 255
				blue *= 255
				heatstr = sprintf("#%02x%02x%02x", red, green, blue)
				return "#{heatstr}-#{opacity}"
			end
		# continuous blue-red colorscale
		elsif color_scheme == "blue-red"
			sx = (b-a == 0) ? 0.5 : Float((x-a)/(b-a))
      opacity = if fade
		      0.1+sx*0.8  # 0.1-0.9
		    else
		      0.75
		    end
			red, green, blue = [0,0,0]
			red = if sx < (1.0/3.0)
					0
				elsif sx > (2.0/3.0)
					1
				else
					# [1/3,2/3] -> [0,1]
					(sx*3)-1
				end
			green = if sx < (1.0/3.0)
					# [0,1/3] -> [0,1]
					sx*3
				elsif sx > (2.0/3.0)
					# [2/3,1] -> [1,0]
					3-(sx*3)
				else
					1
				end
			blue = if sx < (1.0/3.0)
					1
				elsif sx > (2.0/3.0)
					0
				else
					# [1/3,2/3] -> [1,0]
					2-(sx*3)
				end
			red *= 255
			green *= 255
			blue *= 255
			heatstr = sprintf("#%02x%02x%02x", red, green, blue)
			return "#{heatstr}-#{opacity}"
		# orange-red colorscheme with grey area
		elsif color_scheme == "orange-red" && hmin && hmax
			if x >= hmin && x <= hmax
				return "NA"
			elsif x > hmax
				sx = (b-hmax <= 0) ? 0.5 : Float((x-hmax)/(b-hmax))
        opacity = if fade
            0.1+sx*0.8  # 0.1-0.9
          else
            [0.25+sx*0.5, 0.75].max
           end
				red, green, blue = [0,0,0]
				red = 1
				# [0,1] -> [0.5,0]
				green = (1-sx)/2
				red *= 255
				green *= 255
				blue *= 255
				heatstr = sprintf("#%02x%02x%02x", red, green, blue)
				return "#{heatstr}-#{opacity}"
			else
				sx = (hmin-a <= 0) ? 0.5 : Float((x-a)/(hmin-a))
        opacity = if fade
            0.1+sx*0.8  # 0.1-0.9
          else
            [0.25+sx*0.5, 0.75].max
           end
				red, green, blue = [0,0,0]
				red = 1
				# [0,1] -> [1,0.5]
				green = 1-(sx/2)
				red *= 255
				green *= 255
				blue *= 255
				heatstr = sprintf("#%02x%02x%02x", red, green, blue)
				return "#{heatstr}-#{opacity}"
			end
		# continuous orange-red colorscheme
		elsif color_scheme == "orange-red"
#			sx = Float(x/b)
			sx = (b-a == 0) ? 0.5 : Float((x-a)/(b-a))
		  opacity = if fade
          0.1+sx*0.8  # 0.1-0.9
        else
          0.75
        end
			red, green, blue = [0,0,0]
			red = 1
			# [0,1] -> [1,0]
			green = 1-sx
			red *= 255
			green *= 255
			blue *= 255
			heatstr = sprintf("#%02x%02x%02x", red, green, blue)
			return "#{heatstr}-#{opacity}"
		# green colorscale for negative x, red colorscale for positive x
		# grey if within [hmin, hmax] (no information)
		elsif color_scheme == "green-red" && hmin && hmax
			if x >= hmin && x <= hmax
				return "NA"
			elsif x > hmax
				sx = (b-hmax <= 0) ? 0.5 : Float((x-hmax)/(b-hmax))
				opacity = if fade
						0.1+sx*0.8  # 0.1-0.9
					else
					 [0.25+sx*0.5, 0.75].max
					end
				red, green, blue = [0,0,0]
				red = 1
				# [0,1] -> [.5,0]
				green = 0.5-(sx*0.5)
				red *= 255
				green *= 255
				blue *= 255
				heatstr = sprintf("#%02x%02x%02x", red, green, blue)
				return "#{heatstr}-#{opacity}"
			else
				sx = (hmin-a <= 0) ? 0.5 : Float((x-a)/(hmin-a))
        opacity = if fade
            0.1+sx*0.8  # 0.1-0.9
          else
           [0.25+sx*0.5, 0.75].max
          end
				red, green, blue = [0,0,0]
				# [0,1] -> [0,.5]
				red = sx*0.5
				green = 1
				red *= 255
				green *= 255
				blue *= 255
				heatstr = sprintf("#%02x%02x%02x", red, green, blue)
				return "#{heatstr}-#{opacity}"
			end
		# continuous green-red color scheme
		elsif color_scheme == "green-red"
			sx = (b-a==0) ? 0.5 : Float((x-a)/(b-a))
      opacity = if fade
            0.1+sx*0.8  # 0.1-0.9
          else
            0.75
          end
			red, green, blue = [0,0,0]
			red = if sx < 0.5
					# [0,0.5] -> [0,1]
					sx*2
				else
					1
				end
			green = if sx < 0.5
					1
				else
					# [0.5,1] -> [1,0]
					2-(sx*2)
				end
			red *= 255
			green *= 255
			blue *= 255
			heatstr = sprintf("#%02x%02x%02x", red, green, blue)
			return "#{heatstr}-#{opacity}"
		else
			$stderr.puts "unsupported color scheme #{color_scheme}"
			exit -1
		end
		
	end

	### annotation
	g = img.add_element('g', {'id'=>'annotation', 'transform'=>'translate(0,0)'})
	0.upto(ratio.size-1) do |i|
		c = ratio[i]
		if color_scheme == "discrete"
			if discrete_colors.has_key?(c)
				stroke = "\##{discrete_colors[c]}"
				opacity = 0.75
			else
				next
			end
		else
			heatstr = heat_scale(ratio.min, ratio.max, c, color_scheme, hmin, hmax, fade)
			next if heatstr == "NA" # for non-informative positions
			stroke, opacity = heatstr.split("-")
		end
		next if c.nil?
		j = ratio_ind[i]
		g.add_element('circle', {'style'=>"stroke: none; fill: #{stroke}", 'r'=>"7", 'cx'=>"#{outline[j][0]}", 'cy'=>"#{outline[j][1]}", 'fill-opacity'=>"#{opacity}"})
		# display it
		width = [c, 4].min
	end
	### read scale
	if doLegend
		defs = svg.root.add_element('defs')
		if color_scheme == "discrete"
			g = navframe.add_element('svg', {'width'=>'100%', 'height'=>boxheight, 'viewBox'=> sprintf("0 0 120 %d",boxheight), 'preserveAspectRatio'=>'none'})
			ypos = 30
			discrete_colors.keys.sort.each { |value|
				color = discrete_colors[value]
				g.add_element('rect', {"x"=>"10", "y"=>ypos-12, "width"=>"10", "height"=>"15", 'style'=>"fill: \##{color}"})
				txt = g.add_element('text', {"x"=>"25", "y"=>ypos, "font-size"=>"8pt"})
				str = sprintf("%g", value)
				txt.text = str
				ypos += 15
			}
		else
			grad = defs.add_element('linearGradient', {"id"=>"heat", "x1"=>"0%", "y1"=>"0%", "x2"=>"100%", "y2"=>"0%"})
			if color_scheme == "blue-red"
				grad.add_element('stop', {"offset"=>"0%",  "style"=>"stop-color:rgb(0,0,255);stop-opacity:1"})
				grad.add_element('stop', {"offset"=>"33%", "style"=>"stop-color:rgb(0,255,255);stop-opacity:1"})
				if hmin && hmax
					grad.add_element('stop', {"offset"=>"33%", "style"=>"stop-color:rgb(223,223,223);stop-opacity:1"})
					grad.add_element('stop', {"offset"=>"66%", "style"=>"stop-color:rgb(223,223,223);stop-opacity:1"})
				end
				grad.add_element('stop', {"offset"=>"66%", "style"=>"stop-color:rgb(255,255,0);stop-opacity:1"})
				grad.add_element('stop', {"offset"=>"100%","style"=>"stop-color:rgb(255,0,0);stop-opacity:1"})
			elsif color_scheme == "orange-red" && hmin && hmax
				grad.add_element('stop', {"offset"=>"0%",  "style"=>"stop-color:rgb(255,255,0);stop-opacity:1"})
				grad.add_element('stop', {"offset"=>"33%", "style"=>"stop-color:rgb(255,127,0);stop-opacity:1"})
				grad.add_element('stop', {"offset"=>"33%", "style"=>"stop-color:rgb(223,223,223);stop-opacity:1"})
				grad.add_element('stop', {"offset"=>"66%", "style"=>"stop-color:rgb(223,223,223);stop-opacity:1"})
				grad.add_element('stop', {"offset"=>"66%", "style"=>"stop-color:rgb(255,127,0);stop-opacity:1"})
				grad.add_element('stop', {"offset"=>"100%","style"=>"stop-color:rgb(255,0,0);stop-opacity:1"})
			elsif color_scheme == "orange-red"
				grad.add_element('stop', {"offset"=>"0%",  "style"=>"stop-color:rgb(255,255,0);stop-opacity:1"})
				grad.add_element('stop', {"offset"=>"100%","style"=>"stop-color:rgb(255,0,0);stop-opacity:1"})
			elsif color_scheme == "green-red" && hmin && hmax
				grad.add_element('stop', {"offset"=>"0%",  "style"=>"stop-color:rgb(0,255,0);stop-opacity:1"})
				grad.add_element('stop', {"offset"=>"33%", "style"=>"stop-color:rgb(127,255,0);stop-opacity:1"})
				grad.add_element('stop', {"offset"=>"33%", "style"=>"stop-color:rgb(223,223,223);stop-opacity:1"})
				grad.add_element('stop', {"offset"=>"66%", "style"=>"stop-color:rgb(223,223,223);stop-opacity:1"})
				grad.add_element('stop', {"offset"=>"66%", "style"=>"stop-color:rgb(255,127,0);stop-opacity:1"})
				grad.add_element('stop', {"offset"=>"100%","style"=>"stop-color:rgb(255,0,0);stop-opacity:1"})
			elsif color_scheme == "green-red"
				grad.add_element('stop', {"offset"=>"0%",  "style"=>"stop-color:rgb(0,255,0);stop-opacity:1"})
				grad.add_element('stop', {"offset"=>"50%",  "style"=>"stop-color:rgb(255,255,0);stop-opacity:1"})
				grad.add_element('stop', {"offset"=>"100%","style"=>"stop-color:rgb(255,0,0);stop-opacity:1"})
			end

			g = navframe.add_element('svg', {'width'=>'100%', 'height'=>'100', 'viewBox'=>"0 0 120 100", 'preserveAspectRatio'=>'none'})
			txt = g.add_element('text', {"x"=>"10", "y"=>"30", "font-size"=>"8pt"})
			str = sprintf("%.2f", ratio.min);
			txt.text = str
			if (hmin && hmax)
				txt = g.add_element('text', {"x"=>"30", "y"=>"90", "font-size"=>"8pt"})
				str = sprintf("%.1f", hmin);
				txt.text = str
				txt = g.add_element('text', {"x"=>"65", "y"=>"90", "font-size"=>"8pt"})
				str = sprintf("%.1f", hmax);
				txt.text = str
			end
			txt = g.add_element('text', {"x"=>"120", "y"=>"30", "font-size"=>"8pt", "text-anchor"=>"end"})
			str = sprintf("%.2f", ratio.max);
			txt.text = str
			txt = g.add_element('text', {"x"=>"15", "y"=>"15", "font-weight"=>"bold", "font-size"=>"8pt"})
			txt.text = annotationtype.dup
			g.add_element('line', {"x1"=>"11", "y1"=>"35", "x2"=>"11", "y2"=>"50", "style"=>"stroke:rgb(0,0,0);stroke-width:2"})
			g.add_element('line', {"x1"=>"109", "y1"=>"35", "x2"=>"109", "y2"=>"50", "style"=>"stroke:rgb(0,0,0);stroke-width:2"})
			if (hmin && hmax)
				g.add_element('line', {"x1"=>"44", "y1"=>"70", "x2"=>"44", "y2"=>"80", "style"=>"stroke:rgb(0,0,0);stroke-width:2"})
				g.add_element('line', {"x1"=>"75", "y1"=>"70", "x2"=>"75", "y2"=>"80", "style"=>"stroke:rgb(0,0,0);stroke-width:2"})
			end
			g.add_element('rect', {"x"=>"10", "y"=>"50", "width"=>"100", "height"=>"20", "style"=>"fill:url(#heat)"})
		end
	end
end

if doSequence
	img.elements.delete(img.elements[3])
	img.elements.add(seq)
else
	img.elements.delete(img.elements[3])
end

if false
	### subfeature legend
	g = navframe.add_element('svg', {'width'=>'100%', 'height'=>'100%', 'x'=>'0', 'y'=>'100'})
	currY = 10
	subfeature_legend.values.each do |sfl|
		g.add_element('line', {'x1'=>'5', 'x2'=>'15', 'y1'=>currY.to_s, 'y2'=>currY.to_s, 'style'=>"stroke-width:2; stroke:#{sfl.color}"})
		txt = g.add_element('text', {'x'=>'20', 'y'=>currY.to_s, 'dy'=>'0.5ex', 'font-size'=>'10'})
		txt.text = sfl.label

		currY += 20
	end
	navframe.attributes['height'] = (navframe.attributes['height'].to_i + currY).to_s
end

### finish
svg.root.elements.delete(svg.root.elements[2])
puts svg


