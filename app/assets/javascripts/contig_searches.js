jQuery(function() {
    $('#contig_searches').dataTable({
        bProcessing: true,
        bServerSide: true,
        sAjaxSource: $('#contig_searches').data('source'),
        "columnDefs": [{
            "targets": 3,
            "orderable": false
        }],
        "order": [ 2, 'desc' ]
    });

    $('#contig_search_results').dataTable({
        bProcessing: true,
        bServerSide: true,
        sAjaxSource: $('#contig_search_results').data('source'),
        "columnDefs": [{
            "targets": [1, 2, 5],
            "orderable": false
        }],
        "order": [ 0, 'asc' ]
    });

    $('#contig_search_min_age').datepicker({
        dateFormat: 'yy-mm-dd'
    });

    $('#contig_search_max_age').datepicker({
        dateFormat: 'yy-mm-dd'
    });

    $('#contig_search_min_update').datepicker({
        dateFormat: 'yy-mm-dd'
    });

    $('#contig_search_max_update').datepicker({
        dateFormat: 'yy-mm-dd'
    });

    $('#contig_search_name').autocomplete({
        source: $('#contig_search_name').data('autocomplete-source'),
        minLength: 2
    });

    $('#contig_search_marker').autocomplete({
        source: $('#contig_search_marker').data('autocomplete-source')
    });

    $('#contig_search_specimen').autocomplete({
        source: $('#contig_search_specimen').data('autocomplete-source'),
        minLength: 2
    });

    $('#contig_search_taxon').autocomplete({
        source: $('#contig_search_taxon').data('autocomplete-source')
    });

    $('#contig_search_verified_by').autocomplete({
        source: $('#contig_search_verified_by').data('autocomplete-source')
    });
});
