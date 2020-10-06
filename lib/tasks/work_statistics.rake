# frozen_string_literal: true

namespace :data do
  desc 'Do work statistics for <year> and Lab with <labcode> for higher order taxa'
  task :work_statistics_yearly, %i[year labcode] => [:environment] do |_t, args|
    lab_prefixes = { nees: 'gbol', bgbm: 'db' }

    range = Time.new(args[:year].to_i, 1, 1).all_year
    prefix = lab_prefixes[args[:labcode].downcase.to_sym]

    isolates = Isolate.where('lab_isolation_nr ilike ?', "#{prefix}%")
                      .where(created_at: range)
                      .select(:id, :lab_isolation_nr)
    contigs = Contig.where('contigs.name ilike ?', "#{prefix}%")
                    .where(verified_at: range)
                    .select(:id, :name)
    reads = PrimerRead.where('primer_reads.name ilike ?', "#{prefix}%")
                      .where(created_at: range)
                      .select(:id, :name)
    marker_sequences = MarkerSequence.where('marker_sequences.name ilike ?', "#{prefix}%")
                                     .where(created_at: range)
                                     .select(:id, :sequence, :name)

    puts "Calculating activities for #{range.first.year} in lab #{args[:labcode]}...\n"

    puts "Number of new isolates: #{isolates.size}"
    puts "Number of new contigs: #{contigs.size}"
    puts "Number of new primer reads: #{reads.size}"
    puts "Number of new marker sequences: #{marker_sequences.size}"
    puts ''

    # Split output for groups 4 (Marchantiophytina), 9 (Bryophytina), 5 (Filicophytina) and 1 + 10 (Magnoliopsida + Coniferopsida / Spermatophytina)
    group_ids = [1, 4, 5, 9, 10]
    taxa = HigherOrderTaxon.where(id: group_ids).select(:id, :name)

    taxa.each do |taxon|
      isolates_in_taxon = isolates.joins(individual: [species: [family: [order: :higher_order_taxon]]])
                                  .where(individual: { species: { family: { orders: { higher_order_taxon_id: taxon.id } } } })
      puts "Newly created isolates for #{taxon.name}: #{isolates_in_taxon.distinct.size}"

      contigs_in_taxon = contigs.joins(isolate: [individual: [species: [family: [order: :higher_order_taxon]]]])
                                .where(isolate: { individual: { species: { family: { orders: { higher_order_taxon_id: taxon.id } } } } })
      puts "Newly verified contigs for #{taxon.name}: #{contigs_in_taxon.distinct.size}"

      primer_reads_in_taxon = reads.joins(contig: [isolate: [individual: [species: [family: [order: :higher_order_taxon]]]]])
                                   .where(contig: { isolate: { individual: { species: { family: { orders: { higher_order_taxon_id: taxon.id } } } } } })
      puts "Newly uploaded primer reads for #{taxon.name}: #{primer_reads_in_taxon.distinct.size}"

      marker_sequences_in_taxon = marker_sequences.joins(contigs: [isolate: [individual: [species: [family: [order: :higher_order_taxon]]]]])
                                                  .where(contigs: { isolate: { individual: { species: { family: { orders: { higher_order_taxon_id: taxon.id } } } } } })
      base_count = marker_sequences_in_taxon.sum('length(sequence)')
      puts "Number of newly generated barcode sequences for #{taxon.name}: #{marker_sequences_in_taxon.distinct.size}"
      puts "Sequenced base pairs for #{taxon.name}: #{base_count}"
      puts ''
    end
  end

  desc 'Do work statistics'
  task work_statistics_total: [:environment] do |_t, _args|
    gbol_marker = Marker.gbol_marker

    puts 'Calculating activities in total...'

    all_isolates_count = Isolate.all.size
    isolates_bonn_cnt = Isolate.where('lab_isolation_nr ilike ?', '%gbol%').size
    isolates_berlin_cnt = Isolate.where('lab_isolation_nr ilike ?', '%db%').size
    other_isolates_count = all_isolates_count - isolates_bonn_cnt - isolates_berlin_cnt
    puts "Number of isolates: #{all_isolates_count} (total), #{isolates_bonn_cnt} (Bonn), #{isolates_berlin_cnt} (Berlin), #{other_isolates_count} (not assigned to a lab)"

    sequence_cnt_per_marker = {}
    verified_count_per_marker = {}

    gbol_marker.each do |marker|
      sequence_cnt_per_marker[marker.name] = MarkerSequence.where(marker_id: marker.id).size
      verified_count_per_marker[marker.name] = Isolate.joins(contigs: :marker).where(contigs: { verified: true, marker_id: marker.id }).count('contigs.id')
    end

    puts "Number of barcode sequences per marker: #{sequence_cnt_per_marker}"
    puts "Isolates with verified contigs per marker: #{verified_count_per_marker}"

    # isolates_larger_three = Isolate.joins(:marker_sequences).group('isolates.id').having('count(marker_sequences) > 3')
    # isolates_three = Isolate.joins(:marker_sequences).group('isolates.id').having('count(marker_sequences.where(marker.gbol_marker)) = 3')
    # isolates_two = Isolate.joins(:marker_sequences).group('isolates.id').having('count(marker_sequences) = 2')
    # isolates_one = Isolate.joins(:marker_sequences).group('isolates.id').having('count(marker_sequences) = 1')
    # isolates_zero = Isolate.joins(:marker_sequences).group('isolates.id').having('count(marker_sequences) = 0').where(:marker_sequences.first.marker.gbol_marker)
  end

  task :primer_usage, [:marker_name] => [:environment] do |_t, args|
    marker = Marker.find_by_name(args[:marker_name])

    primers = PrimerRead.group(:primer_id).count
    primers_marker = primers.select { |k, _v| Primer.find_by_id(k)&.marker_id == marker.id unless k.nil? }
    primers_marker.keys.each { |k| primers_marker[Primer.find(k).name] = primers_marker[k]; primers_marker.delete(k) }

    puts Hash[primers_marker.sort_by { |_k, v| v }.reverse]
  end

  desc 'Do work statistics for <year> and Lab with <labcode> for all barcode regions'
  task :work_statistics_regions_yearly, %i[year labcode] => [:environment] do |_t, args|
    labcode = args[:labcode].downcase
    year = args[:year].to_i

    if labcode == 'nees' || labcode == 'bgbm'
      lab_prefixes = { nees: 'gbol', bgbm: 'db' }

      range = Time.new(year, 1, 1).all_year
      prefix = lab_prefixes[labcode.to_sym]

      isolates = Isolate.where('lab_isolation_nr ilike ?', "#{prefix}%")
                     .where(created_at: range)
                     .select(:id, :lab_isolation_nr)
      contigs = Contig.where('contigs.name ilike ?', "#{prefix}%")
                    .where(verified_at: range)
                    .select(:id, :name)
      reads = PrimerRead.where('primer_reads.name ilike ?', "#{prefix}%")
                  .where(created_at: range)
                  .select(:id, :name)
      marker_sequences = MarkerSequence.where('marker_sequences.name ilike ?', "#{prefix}%")
                             .where(created_at: range)
                             .select(:id, :sequence, :name)

      puts "Calculating activities for #{range.first.year} in lab #{args[:labcode]}...\n"

      puts "Number of new isolates: #{isolates.size}"
      puts "Number of new contigs: #{contigs.size}"
      puts "Number of new primer reads: #{reads.size}"
      puts "Number of new marker sequences: #{marker_sequences.size}"
      puts ''

      # Split output for markers
      marker_ids = [4, 5, 6, 7]
      markers = Marker.where(id: marker_ids).select(:id, :name)

      markers.each do |marker|
        isolates_for_marker = isolates.joins(marker_sequences: :marker).where(marker_sequences: { marker_id: marker.id })
        puts "Newly created isolates for #{marker.name}: #{isolates_for_marker.distinct.size}"

        contigs_for_marker = contigs.joins(marker_sequence: :marker).where(marker_sequences: { marker_id: marker.id })
        puts "Newly verified contigs for #{marker.name}: #{contigs_for_marker.distinct.size}"

        primer_reads_for_marker = reads.joins(contig: [ marker_sequence: :marker ])
                                    .where(contigs: { marker_sequences: { marker_id: marker.id }})
        puts "Newly uploaded primer reads for #{marker.name}: #{primer_reads_for_marker.distinct.size}"

        marker_sequences_for_marker = marker_sequences.where(marker_id: marker.id)
        base_count = marker_sequences_for_marker.sum('length(sequence)')
        puts "Number of newly generated barcode sequences for #{marker.name}: #{marker_sequences_for_marker.distinct.size}"
        puts "Sequenced base pairs for #{marker.name}: #{base_count}"
        puts ''
      end
    else
      puts 'Labcode must be one of BGBM or NEES.'
    end
  end
end
