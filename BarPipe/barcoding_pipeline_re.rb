# frozen_string_literal: true

require 'bio'
require 'tre-ruby'
require 'fuzzystringmatch'
require 'nokogiri'
require 'parallel'
require 'csv'
require 'optparse'
require 'ostruct'
require 'pp'

require_relative 'barcoding_pipeline_re_methods'

################
# USER ARGUMENTS
################

options = OptparseExample.parse(ARGV)
options[:output] = ARGV[0].gsub(%r{/$}, '')

#############
# VALIDATIONS
#############

if options[:map].length > 1
  if options['tpm_splitter'].nil?
    raise InputError, 'Please provide tags for assigning the according tagPrimerMap.'
  end

  multiple_plates = true
else
  multiple_plates = false
end

if !options[:fa] && !options[:webdav]
  raise InputError, 'You did not provide a fastq/fasta file or a webdav address. Please do so when running the script again.'
end

if options[:fa] && options[:webdav]
  raise InputError, 'Please provide EITHER a webdav address or a fastq/fasta file.'
end

###########
# VARIABLES
###########

job_title = File.basename(options[:output])

intermediate_out_dir = "#{options[:output]}/intermediate_out"
results_dir = "#{options[:output]}/results"
data_folder = "#{options[:output]}/data"
fastq_folder = "#{data_folder}/fastq"
fasta_folder = "#{data_folder}/fasta"
webdav_folder = "#{data_folder}/webdav"

homopolymers = 0
too_long_too_short = 0
ambiguous_seqs_counter = 0
no_second_primer_counter = 0
barcode_counter = {}
no_barcode_counter = 0
more_than_one_hit = 0

#############
# DIR PRESETS
#############

if Dir.exist?(intermediate_out_dir)
  `rm -r #{intermediate_out_dir}/*parsed 2> /dev/null`
end

`rm -r #{results_dir} 2> /dev/null`
`mkdir -p #{intermediate_out_dir} #{results_dir} #{fasta_folder} #{fastq_folder}`

######
# PIPE
######

begin

  ################
  # TPM PROCESSING

  tag_primer_maps = {}
  sample_id_to_description = {}
  plate_names = []

  options[:map].map! do |tpm|
    plate_names << tpm.scan(/(?<!\.|\D)[^._\d]+/)[-1]
    plate = tpm.scan(/[0-9]+[0-9A-Za-z]*/)[-1]
    tpm_processing = TpmProcessing.new

    # TPM CORRECTION
    altered_tpm = tpm_processing.alter_tpm(tpm)

    # TPM READING
    processed_tpm = tpm_processing.read_tpm(options[:primer_mismatches])
    tag_primer_maps[plate] = processed_tpm[0]
    sample_id_to_description = sample_id_to_description.merge(processed_tpm[1])

    # save altered tpm in options[:map]
    altered_tpm
  end

  ####################
  # REGION DIR PRESETS

  regions_array = tag_primer_maps.map { |_plate, plate_info| plate_info[:region] }.uniq
  regions_array.each do |region|
    `mkdir '#{intermediate_out_dir}/#{region}_parsed' 2> /dev/null`
    `mkdir '#{intermediate_out_dir}/#{region}_clusters' 2> /dev/null`
    `mkdir '#{intermediate_out_dir}/#{region}_blast' 2> /dev/null`
  end

  ########################
  # SEQUEL DATA PROCESSING

  if options[:webdav] && Dir["#{fasta_folder}/*"].empty?
    sequel_processing = SequelDataProcessing.new(data_folder, job_title)

    # download from webdav and unzip
    sequel_processing.retrieve_webdav(options[:webdav])

    bam_exists = !`find #{webdav_folder} -name [^.]*.subreads.bam`.chomp.empty?
    if bam_exists
      # run ccs, lima, bam2fq
      sequel_processing.bam2fastq_wrapper(bam, multiple_plates)
    end
  end

  ###################
  # QUALITY FILTERING

  if options[:mode] == 1
    raise FinishedPrematurely, 'Finished sequel data processing. Job finished successfully.'
  end

  # TODO: think about approach again
  # TODO: we already have a fasta and fastq dir, i.e. the file types separated
  is_fastq = !`head -n 4 '#{options[:fa]}'`.scan(/^@.*\n.*\n\+/).empty?
  if options[:fa] && is_fastq
    `seqtk seq -a -q #{options[:quality_stat]} '#{options[:fa]}' > '#{fasta_folder}/seqtk_filtered.fasta'`
  elsif !Dir["#{fastq_folder}/*"].empty?
    Dir["#{fastq_folder}/*"].each do |fastq_file|
      `seqtk seq -A #{fastq_file} > #{fasta_folder}/#{job_title}.fasta`
    end
  end

  options[:filtered_fasta] = Dir["#{fasta_folder}/*.fasta"]

  seqs_pre_filtering = 'NA' if options[:webdav] # TODO: why???
  seqs_post_filtering_seqtk = `grep -c "^>" #{options[:filtered_fasta].join(' ')}`
                              .split("\n").map { |c| c.rpartition(':')[-1].to_i }

  ######################################
  # QIIME REPLACEMENT AKA DEMULTIPLEXING

  if options[:mode] == 2
    raise FinishedPrematurely, 'Finished quality filtering. Job finished successfully.'
  end

  no_reverse = 0
  one_primer_only = {}
  seq_to_id = {}
  multiple_barcodes = {}

  DemuxTime.out_dir(intermediate_out_dir)
  options[:filtered_fasta].each do |fasta|
    plate = fasta.scan(/[0-9]+[0-9A-Za-z]*/)[-1]

    Bio::FlatFile.open(Bio::FastaFormat, fasta).each do |entry|
      sequence = entry.seq
      demultiplexing = DemuxTime.new(entry,
                                     plate,
                                     options[:tag_mismatches],
                                     tag_primer_maps[plate])

      # quality check
      # TODO: doesn't work that way bc quality check will be called several times inside DemuxTime
      length_check = sequence.length.between?(150, 7500)
      homopolymer_check = sequence.scan(/(\w)\1{31,}/).empty?

      homopolymers += 1 unless homopolymer_check
      too_long_too_short += 1 unless length_check

      if length_check && homopolymer_check
        demultiplexing.demultiplex
      end
    end
  end

  ###################
  # USEARCH AND BLAST

  if options[:mode] == 3
    raise FinishedPrematurely, 'Finished demultiplexing. Job finished successfully.'
  end

  seq_name_to_lineage = {}
  cluster_info = {}
  clusters_per_sample = {}
  cluster_nums = {}
  e_values = {}

  regions_array.each do |region|
    usearch_output_dir = "#{intermediate_out_dir}/#{region}_clusters"
    blast_output_dir = "#{intermediate_out_dir}/#{region}_blast"
    region_dir = Dir["#{intermediate_out_dir}/#{region}_parsed/*"]

    region_dir.each do |parsed_file|
      sample_id = parsed_file.split(%r{[/.]})[-2]

      unless File.exist?("#{usearch_output_dir}/log_#{sample_id}")
        `'/data/data2/lara/Barcoding/usearch11.0.667_i86linux32' -id #{options[:identity]} -cluster_fast '#{parsed_file}' -quiet -centroids '#{usearch_output_dir}/#{sample_id}_centroid' -uc '#{usearch_output_dir}/log_#{sample_id}' -consout '#{usearch_output_dir}/#{sample_id}_consensus'`
      end
      unless File.exist?("#{blast_output_dir}/#{sample_id}")
        threads = options[:threads] || ServerJobManager.calculate_num_threads
        `blastn -num_threads #{threads} -db nt -max_target_seqs 1 -outfmt '6 qseqid sskingdoms sscinames staxids evalue length pident nident mismatch gaps sacc sseqid stitle' -evalue 1e-10 -query '#{usearch_output_dir}/#{sample_id}_centroid' -out '#{blast_output_dir}/#{sample_id}'`
      end

      # how many clusters?
      clusters_per_sample["#{region}#{sample_id}"] = `grep -c "^C" '#{usearch_output_dir}/log_#{sample_id}'`.to_i

      # retrieve ids
      seq_name = `awk -F"\t" '{print $1}' '#{blast_output_dir}/#{sample_id}'`.split("\n")
      tax_ids = `awk -F"\t" '{print $4}' '#{blast_output_dir}/#{sample_id}'`.split("\n")
      e_value = `awk -F"\t" '{print $5}' '#{blast_output_dir}/#{sample_id}'`.split("\n")
      seq_name_to_id = seq_name.zip(tax_ids).to_h

      times = *(0...e_value.length)
      times.each do |i|
        e_values[seq_name[i]] = e_value[i]
      end

      seq_name_to_id.each do |seqname, id|
        lineage = MyNokogiri.lineage(id)
        seq_name_to_lineage[seqname] = lineage
      end
    end

    `cat '#{usearch_output_dir}/log_'* > '#{usearch_output_dir}/clusterlog_#{region}.txt'`

    # retrieve usearch log file information
    CSV.foreach("#{usearch_output_dir}/clusterlog_#{region}.txt'", col_sep: "\t") do |row|
      record_type = row[0]
      query_seq = row[8]
      target_seq = row[9]
      cluster_number = row[1]

      if record_type == 'C'
        sample_id = query_seq.split('_')[0]
        new_id = "#{sample_id}_#{region}_#{cluster_number}"
        cluster_info[new_id] = [row[2], query_seq, seq_name_to_lineage[query_seq]]
      end

      cluster_nums[query_seq] = cluster_number
      cluster_nums[target_seq] = cluster_number
    end
  end

  ########################
  # O1: tagPrimerMap (ext)

  options[:map].each do |tpm|
    plate = tpm.scan(/[0-9]+[0-9A-Za-z]*/)[-1]
    this_tpm = tag_primer_maps[plate]
    CSV.open("#{results_dir}/tagPrimerMap_#{plate}_expanded.txt", 'w+', col_sep: "\t") do |out|
      csv = CSV.read(map, headers: true, col_sep: "\t")
      out << csv.headers + %w[HighQualSeqs LinkerPrimerOnly Clusters TotalSequences]

      csv.each do |row|
        region = row['Region']
        description = row['Description']
        id = row['#SampleID'].scan(/(\D+\d+|\D+\z)/)
        gbol = id[0][0]
        id_plus_marker = id.join('_')

        if File.exist?("#{intermediate_out_dir}/#{region}_parsed/#{gbol}.fasta")
          high_qual = `grep -c "^>" '#{intermediate_out_dir}/#{region}_parsed/#{gbol}.fasta'`.chomp.to_i
          # evtl anders berechnen, wenn gbol IDs auch plattenübergreifend auftauchen können
        else
          high_qual = 0
        end
        linker_primer_only = one_primer_only[description]
        linker_primer_only = 0 if linker_primer_only == nil
        clusters = clusters_per_sample["#{region}#{gbol}"]
        out << [id_plus_marker] + row[1..-1] + [high_qual, linker_primer_only, clusters, high_qual + linker_primer_only]
      end
    end
  end

  ##################
  # O2: cluster info

  # für jedes cluster drei Spalten: no.seqs, centroid seq, taxFilterResult
  File.open("#{results_dir}/#{job_title}_cluster_log.txt", 'w+') do |out|
    out.puts "cluster\tno_seqs\tcentroid_seq\tTaxFilterResult\tfamily\tgenus"
    cluster_info.each do |number, value|
      out.puts "#{number}\t#{value[0]}\t#{value[1]}\t#{value[2]}"
    end
  end

  ##################
  # O3: barcode info

  File.open("#{results_dir}/barcode_log.txt", 'w+') do |out|
    out.puts "barcode\ttimes_found"
    barcode_counter.sort_by{ |k, v| -v }.each do |barcode, count|
      out.puts "#{barcode}\t#{count}"
    end
  end

  #######################
  # O4: multiple barcodes

  File.open("#{results_dir}/multiple_barcodes.txt", 'w+') do |out|
    out.puts "sequence_id\tno_of_barcodes\tfound_barcodes"
    multiple_barcodes.sort_by { |_k, v| -v.length }.each do |seq_id, barcodes|
      count = barcodes.length
      str = ''
      barcodes.each do |bar|
        str = "#{str}#{bar},"
      end
      out.puts "#{seq_id}\t#{count}\t#{str[0...-1]}"
    end
  end

  #######################
  # O5: Multifasta/region

  info_for_consensus = {}
  last_descr = nil

  regions_array.each do |region|
    input_fasta_dir = Dir["#{intermediate_out_dir}/#{region}_parsed/*"]
    File.open("#{results_dir}/#{region}.fasta", 'w+') do |out|
      input_fasta_dir.each do |input_fasta_file|

        Bio::FlatFile.open(Bio::FastaFormat, input_fasta_file).each do |entry|
          seqname = entry.definition
          splitted_seqname = seqname.split('|')
          gbol_descr = splitted_seqname[0]
          sequence_id = splitted_seqname[1]
          rc_info = splitted_seqname[-1]
          cluster_num = cluster_nums[entry.definition]
          gbol = gbol_descr.split('_')[0]
          this_cluster_info = cluster_info["#{gbol}_#{region}_#{cluster_num}"]

          begin
            lineage = this_cluster_info[2]
            centroid_seqid = this_cluster_info[1]
            cluster_size = this_cluster_info[0]
          rescue
            lineage = 'NA'
            centroid_seqid = 'NA'
            cluster_size = 'NA'
            File.open("#{results_dir}/#{job_title}_stderr.txt", 'a+') do |out|
              out.puts "===\nNo cluster information for the following entry:"
              out.puts entry.definition
              out.puts "#{entry.seq}\n"
            end
          end

          lineage_to_print = if lineage && lineage != 'NA'
                               if lineage.include?(options[:taxon])
                                 lineage.split.join(',')
                               else
                                 lineage.split("\t")[0]
                               end
                             else
                               'NA'
                             end

          if gbol_descr == last_descr
            num += 1
          else
            num = 0
          end

          last_descr = gbol_descr
          gbol_descr += "_#{num}"

          # obtain cluster size
          if cluster_num
            cluster_num = "#{cluster_num}*centroid" if seqname == centroid_seqid
          end

          evalue = e_values[centroid_seqid]

          definition = ">#{gbol_descr}|#{cluster_num}(#{cluster_size})|#{lineage_to_print}|#{evalue}|#{sequence_id}|#{rc_info}"
          out.puts definition
          out.puts entry.seq

          if definition.include?(options[:taxon])
            File.open("#{results_dir}/#{region}_#{options[:taxon]}.fasta", 'a+') do |e_out|
              e_out.puts definition
              e_out.puts entry.seq
            end
          else
            File.open("#{results_dir}/#{region}_non-#{options[:taxon]}.fasta", 'a+') do |no_e_out|
              no_e_out.puts definition
              no_e_out.puts entry.seq
            end
          end

          begin
            info_for_consensus["#{gbol}#{region}#{cluster_num.split('*')[0]}"] = [cluster_size, evalue, lineage_to_print]
          rescue
            info_for_consensus["#{gbol}#{region}#{cluster_num}"] = [cluster_size, evalue, lineage_to_print]
          end
        end
      end
    end
  end

  ####################
  # O6: consensus seqs

  regions_array.each do |region|
    consensus_files = Dir["#{intermediate_out_dir}/#{region}_clusters/*_consensus"]
    consensus_files.each do |consensus_file|
      gbol = consensus_file.scan(%r{(?<=/)[^/]+(?=_)})[-1]

      Bio::FlatFile.open(Bio::FastaFormat, consensus_file).each do |entry|
        cluster_num = entry.definition.scan(/\d+/)[0]
        gbol_definition = sample_id_to_description["#{gbol}#{region}"]
        info = info_for_consensus["#{gbol}#{region}#{cluster_num}"]
        cluster_size = info[0]
        taxonomy = info[2]
        evalue = info[1]

        if taxonomy.include?(options[:taxon])
          File.open("#{results_dir}/#{region}_consensus_#{options[:taxon]}.fasta", 'a+') do |e_out|
            e_out.puts ">#{gbol_definition}|#{cluster_num}*consensus(#{cluster_size})|#{taxonomy}|#{evalue}"
            e_out.puts entry.seq
          end
        else
          File.open("#{results_dir}/#{region}_consensus_non-#{options[:taxon]}.fasta", 'a+') do |no_e_out|
            no_e_out.puts ">#{gbol_definition}|#{cluster_num}*consensus(#{cluster_size})|#{taxonomy}|#{evalue}"
            no_e_out.puts entry.seq
          end
        end
      end
    end
  end

  ################
  # O7: plate info

  # TODO: if FASTQ! Should be able to handle FASTA as well!
  seqs_pre_filtering = `wc -l '#{options[:fa]}'`.to_i / 4 unless seqs_pre_filtering == 'NA'
  seqs_post_filtering_qiime = 0
  regions_array.each do |region|
    seqs_post_filtering_qiime += `grep -c "^>" '#{results_dir}/#{region}.fasta'`.to_i
  end

  File.open("#{results_dir}/#{job_title}_log.txt", 'w+') do |out|
    out.puts "Seqs pre-filtering: #{seqs_pre_filtering}"
    out.puts "Seqs post-seqtk-filtering: #{seqs_post_filtering_seqtk}"
    out.puts "Seqs post-qiimelike-filtering (highQual): #{seqs_post_filtering_qiime}"
    out.puts "Seqs with more than one hit: #{more_than_one_hit}"
    out.puts "Seqs w/o barcode: #{no_barcode_counter}"
    out.puts "Seqs w/o 2nd primer: #{no_second_primer_counter}"
    out.puts "Ambiguous seqs: #{ambiguous_seqs_counter}"
    out.puts "Seqs too short/too long: #{too_long_too_short}"
    out.puts "Homopolymers: #{homopolymers}"
  end

  ##########
  # GBOL APP

  # zip all output files
  Dir.chdir(options[:output]){
    %x[zip -r '#{job_title}_out.zip' '#{job_title}_out']
  }

  `/data/data2/lara/send_telegram.sh "Barcoding pipe finished on $(date)"`
  if options[:id]
    `curl -X POST "https://gbol5.de/ngs_runs/#{id}/import?results_path=#{options[:output]}/#{job_title}_out.zip" -u "external_user:K25#+95=@qmVgA5K"`
  end

rescue => e
  File.open("#{results_dir}/stderr.txt", 'a+') do |out|
    out.puts "===\nReason for abort:"
    out.puts e.message
    out.puts e.backtrace
  end
  `/data/data2/lara/send_telegram.sh "An error occurred."`
end