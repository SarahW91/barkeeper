jQuery(function() {
    $('#new_primer_read').fileupload(
        {
            dataType: "script",
            add: function(e, data) {
                data.context = $(tmpl("template-upload", data.files[0]));
                $('#new_primer_read').append(data.context);
                return data.submit();
            },
            progress: function(e, data) {
                var progress;
                if (data.context) {
                    progress = parseInt(data.loaded / data.total * 100, 10);
                    return data.context.find('.progress-bar.progress-bar-striped.active').css('width', progress + '%');
                }
            }
        }
    );
    $('#primer_reads').DataTable({
        bProcessing: true,
        bServerSide: true,
        sAjaxSource: $('#primer_reads').data('source'),
        "columnDefs": [
            { "orderable": false, "targets": 2 },
            { "orderable": false, "targets": 4 }
        ],
        "order": [ 3, 'desc' ]
    });
    $('#primer_read_contig_name').autocomplete( {
        source: $('#primer_read_contig_name').data('autocomplete-source')
    });

//    draw chromatogram
    var chromatogram1 = $('#chromatogram').data('url');

    if (chromatogram1){
        draw_chromatogram(chromatogram1);
    }

    //$('#chromatogram_container').scrollLeft(100);

});


function draw_chromatogram(chromatogram1){

    var change_base_primer_read_url='/primer_reads/'+chromatogram1.id+'/change_base';

    var ymax=250;

    var scale=4;

    var lineFunction = d3.svg.line()
        .x(function(d,i) { return i; })
        .y(function(d) { return ymax-d/scale; })
        .interpolate("linear");

    var svg=d3.select('#chromatogram')
        .append('svg')
        .attr('width', chromatogram1.atrace.length)
        .attr('height', 250);

    //draw clipped areas
    if (chromatogram1.trimmedReadStart){

        svg.append('rect')
            .attr("x", 0)
            .attr("y", 0)
            .attr("width", chromatogram1.peak_indices[chromatogram1.trimmedReadStart-1]-5)
            .attr("height", ymax)
            .attr("fill", "#d3d3d3");

        svg.append('rect')
            .attr("x", chromatogram1.peak_indices[chromatogram1.trimmedReadEnd-1]+5)
            .attr("y", 0)
            .attr("width", chromatogram1.atrace.length-chromatogram1.peak_indices[chromatogram1.trimmedReadEnd-1]+5)
            .attr("height", ymax)
            .attr("fill", "#d3d3d3");
    }


    //draw traces

    svg.append("path")
        .attr("d", lineFunction(chromatogram1.atrace))
        .attr("stroke", "green")
        .attr("stroke-width", 1)
        .attr("fill", "none");
    svg.append("path")
        .attr("d", lineFunction(chromatogram1.ctrace))
        .attr("stroke", "blue")
        .attr("stroke-width", 1)
        .attr("fill", "none");
    svg.append("path")
        .attr("d", lineFunction(chromatogram1.gtrace))
        .attr("stroke", "black")
        .attr("stroke-width", 1)
        .attr("fill", "none");
    svg.append("path")
        .attr("d", lineFunction(chromatogram1.ttrace))
        .attr("stroke", "red")
        .attr("stroke-width", 1)
        .attr("fill", "none");


    //draw base calls

    for(var i = 0; i < chromatogram1.peak_indices.length; i++){
        var pos = chromatogram1.peak_indices[i];
        var ch = chromatogram1.sequence[i];

        var color='gray';
        var ta='middle';

        // position indicator
        var disp = i+1;

        if(disp % 10 == 0){

            svg.append("text")
                .attr("x", pos)
                .attr("y", 10)
                .text(disp)
                .attr("font-family", "sans-serif")
                .attr("font-size", "7px")
                .attr("fill", color)
                .attr("text-anchor", ta);
            svg.append("text")
                .attr("x", pos)
                .attr("y", 17)
                .text('.')
                .attr("font-family", "sans-serif")
                .attr("font-size", "7px")
                .attr("fill", color)
                .attr("text-anchor", ta);
        }

        //base calls
        if (ch == 'A') {
            color = 'green';
        } else if (ch == 'C') {
            color = 'blue';
        } else if (ch == 'G') {
            color = 'black';
        } else if (ch == 'T') {
            color = 'red';
        } else {
            color = 'gray';
        }

        svg.append("text")
            .attr("x", pos)
            .attr("y", 30)
            .text(ch)
            .attr("font-family", "sans-serif")
            .attr("font-size", "10px")
            .attr("fill", color)
            .attr("text-anchor", ta)
            .attr("id", i)

            .on('mouseover', function(){
                d3.select(this)
                    .style('font-size','14px')
                    .style('font-weight', 'bold')
            })
            .on('mouseout', function(){
                d3.select(this)
                    .style('font-size','10px')
                    .style('font-weight', 'normal')
            })
            .on('click', function(){
                var p = this.parentNode;


                var selected_base = d3.select(this);
                var p_el = d3.select(p);


                var current_x=selected_base.attr("x");
                var current_y=selected_base.attr("y");
                var current_char=selected_base.text();
                var base_index = selected_base.attr("id");

                var frm = p_el.append("foreignObject");

                var inp = frm
                    .attr({
                        'x': current_x-5,
                        'y': 12,
                        'width': 20,
                        'height': 20
                    })
                    .append("xhtml:form")
                    .append('xhtml:input')
                    .attr("value", current_char)
                    .on("keypress", function() {

                        if (d3.event.keyCode===13) {

                            event.preventDefault(); // cancel default behavior

                            var newBase = inp.node().value;

                            if (newBase==" " || newBase=="" || newBase=="_") {
                                newBase = "-";
                            }

                            change_base(base_index, newBase, change_base_primer_read_url);

                            selected_base.text(newBase);

                            if (newBase=="A"){
                                selected_base.attr("fill", 'green');
                            } else if (newBase=="C"){
                                selected_base.attr("fill", 'blue');
                            } else if (newBase=="G"){
                                selected_base.attr("fill", 'black');
                            } else if (newBase=="T"){
                                selected_base.attr("fill", 'red');
                            } else {
                                selected_base.attr("fill", 'grey');
                            }


                            frm.remove();
                        }
                    });

            });

        color='gray';

        //quality scores
        var q=chromatogram1.qualities[i];
        //ignore manually entered bases with fake qualities "-10"
        if (q > -10) {
            svg.append("text")
                .attr("x", pos)
                .attr("y", 40)
                .text(q)
                .attr("font-family", "sans-serif")
                .attr("font-size", "7px")
                .attr("fill", color)
                .attr("text-anchor", ta);
        }
    }
}

function change_base(base_index, base, change_base_primer_read_url) {
    //console.log(pos, base, change_base_primer_read_url);
    $.ajax({
        data: {
            'position': base_index,
            'base': base
        },
        type: 'POST',
        url: change_base_primer_read_url,
        success: function () {
        },
        error: function (response) {
            // we had an error
            alert('Not authorized? Could not change base at index '+base_index+' to '+base);
        }
    });
    return 0;
}

//var resize = d3.behavior.drag()
//    .origin(function() {
//        var current = d3.select(this);
//        return {x: current.attr("x"), y: current.attr("y") };
//    })
//    .on("drag", dragResize);
//
//function dragResize(){
//    var dragx = Math.max(dx + (16/2), Math.min(w, dx + width + d3.event.dx));
//
//
//    //console.log("resize x:"+x+" y:"+y);
//    console.log("d3.event.x:"+d3.event.dx);
//
//    var dragTarget = d3.select(this);
//    var dragObject = d3.select(this.parentNode);
//
//    var o = dragObject.select("rect.box");
//    var o1 = dragObject.select("rect.titleBox");
//
//
//    var oldx = dx;
//    var oldy = dy;
//
//    dx = Math.max(0, Math.min(dx + width - (16 / 2), d3.event.x));
//    dy = Math.max(0, Math.min(dy + height - (16 ), d3.event.y));
//    w = w - (oldx - dx);
//    h = h - (oldy - dy);
//
//
//    dragTarget
//        .attr("x", function(d) { return dragx - (16/2) })
//        .attr("y", function(d) { return dragy - (16) })
//
//    o.attr("width", w)
//        .attr("height", h);
//
//    o1.attr("width", w);
//};